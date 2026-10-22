/* Checks of the C client (api/c/science.h), built by api/Makefile as a C99
 * program linked with libscience.a, and run by api/tests/run-tests.mjs. Offline,
 * the addresses refused and the JSON readers (numbers with exponents, escapes,
 * surrogate pairs, nested paths, texts that are not JSON). Against a running
 * service, its answers held to the Python client's for the same requests, read
 * on standard input as check_cpp.py writes them: samples, a series and grids
 * exactly, calls of the domains exactly; then a source defined, a reading
 * posted and received on a stream (in a second thread) within a second, a
 * reading refused with its reason, and the error codes. One JSON object a line,
 * as check_client.py reports.
 *
 *     check_c --offline
 *     check_c <url> <sensor source> [token] < reference.json
 */

#define _POSIX_C_SOURCE 200809L

#include "science.h"

#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ---- reporting ------------------------------------------------------------------------------------ */

static void quoted(const char *s) {
    putchar('"');
    for (; s && *s; s++) {
        unsigned char c = (unsigned char)*s;
        if (c == '"' || c == '\\') printf("\\%c", c);
        else if (c == '\n') fputs("\\n", stdout);
        else if (c < 0x20) printf("\\u%04x", c);
        else putchar(c);
    }
    putchar('"');
}

static void affirm(const char *name, int ok, const char *detail) {
    fputs("{\"kind\":\"assert\",\"name\":\"c client: ", stdout);
    fputs(name, stdout);
    printf("\",\"ok\":%s,\"detail\":", ok ? "true" : "false");
    quoted(detail ? detail : "");
    puts("}");
    fflush(stdout);
}

static void check(const char *name, double value, double limit, const char *unit) {
    printf("{\"kind\":\"check\",\"name\":\"c client: %s\",\"value\":", name);
    if (isfinite(value)) printf("%.17g", value);
    else fputs("null", stdout);
    printf(",\"limit\":%.17g,\"unit\":\"%s\"}\n", limit, unit);
    fflush(stdout);
}

static void skip(const char *what) {
    printf("{\"kind\":\"skip\",\"name\":\"c client: %s\"}\n", what);
    fflush(stdout);
}

/* The largest difference of two lists, NaN matching NaN (else infinite). */
static double difference(const double *a, const double *b, size_t n) {
    double worst = 0;
    for (size_t i = 0; i < n; i++) {
        if (isnan(a[i]) || isnan(b[i])) {
            if (isnan(a[i]) != isnan(b[i])) return INFINITY;
            continue;
        }
        if (fabs(a[i] - b[i]) > worst) worst = fabs(a[i] - b[i]);
    }
    return worst;
}

static double seconds_now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + t.tv_nsec * 1e-9;
}

/* An array of numbers at a path of the reference, allocated; NULL when there is none. */
static double *numbers(const char *json, const char *path, size_t *count) {
    *count = 0;
    if (science_json_numbers(json, path, NULL, 0, count) != SCIENCE_OK) return NULL;
    double *out = malloc((*count ? *count : 1) * sizeof *out);
    if (out && science_json_numbers(json, path, out, *count, count) != SCIENCE_OK) {
        free(out);
        return NULL;
    }
    return out;
}

/* A text at a path, compared with what it should be. */
static int text_is(const char *json, const char *path, const char *want) {
    char *got = NULL;
    int same = science_json_text(json, path, &got) == SCIENCE_OK && strcmp(got, want) == 0;
    science_free(got);
    return same;
}

/* ---- offline ---------------------------------------------------------------------------------------- */

static void offline(void) {
    science_client *c = science_connect("https://example.org", NULL, 5);
    affirm("an https address is refused, saying why", !c && strstr(science_last_error(), "TLS proxy") != NULL, science_last_error());
    c = science_connect("ftp://example.org", NULL, 5);
    affirm("an address that is not http:// is refused", !c, science_last_error());
    c = science_connect("http://example.org:8790", "token", 0);
    affirm("an http address makes a client, without a request", c != NULL && science_last_error()[0] == '\0', science_last_error());
    science_close(c);

    const char *doc = "{\"values\": [1e5, -2.5E-3, null, 6.02214076e23, 0], \"grid\": {\"nlon\": 360, \"rows\": [[1, 2], [3, {\"x\": -0.5}]]},"
                      " \"s\": \"\\ud83d\\ude00 \\u00e9\\n\\\"q\\\"\", \"lone\": \"\\ud83d!\", \"t\": true}";
    double v[8];
    size_t n = 0;
    int ok = science_json_numbers(doc, "values", v, 8, &n) == SCIENCE_OK && n == 5 && v[0] == 1e5 && v[1] == -2.5e-3 && isnan(v[2]) &&
             v[3] == 6.02214076e23 && v[4] == 0;
    affirm("json: an array of numbers with exponents, null as NaN", ok, science_last_error());
    double x = 0, y = 0;
    ok = science_json_number(doc, "grid.nlon", &x) == SCIENCE_OK && x == 360 && science_json_number(doc, "grid.rows.1.1.x", &y) == SCIENCE_OK && y == -0.5;
    affirm("json: numbers at nested paths of keys and indices", ok, science_last_error());
    affirm("json: escapes, and a surrogate pair as one character in UTF-8", text_is(doc, "s", "\xF0\x9F\x98\x80 \xC3\xA9\n\"q\""), NULL);
    affirm("json: a lone surrogate as U+FFFD", text_is(doc, "lone", "\xEF\xBF\xBD!"), NULL);
    affirm("json: any other value as its JSON text", text_is(doc, "grid.rows.1", "[3,{\"x\":-0.5}]") && text_is(doc, "t", "true"), NULL);
    affirm("json: nothing at a path is an argument error", science_json_number(doc, "grid.none", &x) == SCIENCE_ERROR_ARGUMENT &&
                                                                science_json_number(doc, "values.9", &x) == SCIENCE_ERROR_ARGUMENT &&
                                                                science_json_number(doc, "s", &x) == SCIENCE_ERROR_ARGUMENT,
           science_last_error());
    affirm("json: a text that is not JSON is an argument error, saying where",
           science_json_number("{\"a\": [1, 2,, 3]}", "a.0", &x) == SCIENCE_ERROR_ARGUMENT && strstr(science_last_error(), "byte 12") != NULL,
           science_last_error());
}

/* ---- against a service ------------------------------------------------------------------------------- */

struct listening {
    science_client *api;
    const char *source;
    volatile int connected;
    double deadline, received;
    char station[64], id[64];
    int code;
};

static int on_event(const char *id, const char *source, const char *readings, void *user) {
    struct listening *l = user;
    (void)source;
    if (!readings) { /* a heartbeat */
        l->connected = 1;
        return seconds_now() < l->deadline;
    }
    l->received = seconds_now();
    char *station = NULL;
    if (science_json_text(readings, "0.station", &station) == SCIENCE_OK) snprintf(l->station, sizeof l->station, "%s", station);
    science_free(station);
    snprintf(l->id, sizeof l->id, "%s", id);
    return 0;
}

static void *listen_thread(void *user) {
    struct listening *l = user;
    l->code = science_stream(l->api, l->source, NULL, on_event, l);
    return NULL;
}

static void fields(science_client *api, const char *ref) {
    size_t n = 0, m = 0, k = 0;
    double *lon = numbers(ref, "sample.lon", &n), *lat = numbers(ref, "sample.lat", &m), *want = numbers(ref, "sample.values", &k);
    if (!lon || !lat || !want) {
        skip("fields: gfs-forecast is not in the cache");
        free(lon), free(lat), free(want);
        return;
    }
    char *time = NULL;
    science_json_text(ref, "sample.time", &time);
    double *got = malloc(n * sizeof *got);
    int code = science_sample(api, "gfs-forecast", "tmp2m", time, lon, lat, n, got);
    if (code == SCIENCE_OK) check("fields: samples equal the Python client's (ten places)", difference(got, want, n), 0, "K");
    else affirm("fields: samples", 0, science_last_error());
    free(lon), free(lat), free(want), free(got);

    lon = numbers(ref, "many.lon", &n), lat = numbers(ref, "many.lat", &m), want = numbers(ref, "many.values", &k);
    got = malloc(n * sizeof *got);
    code = science_sample(api, "gfs-forecast", "rh2m", time, lon, lat, n, got);
    if (code == SCIENCE_OK) check("fields: 600 places equal the Python client's", difference(got, want, n), 0, "%");
    else affirm("fields: 600 places", 0, science_last_error());
    free(lon), free(lat), free(want), free(got);
    science_free(time);

    size_t count = 0;
    double *times = NULL, *values = NULL, *want_times = numbers(ref, "series.times_unix", &k), *want_values = numbers(ref, "series.values", &m);
    code = science_series(api, "gfs-forecast", "tmp2m", 139.69, 35.69, &count, &times, &values);
    if (code == SCIENCE_OK && count == m) {
        check("fields: a series equals the Python client's", difference(values, want_values, count), 0, "K");
        check("fields: ... at the same times (Unix seconds)", difference(times, want_times, count), 0, "s");
    } else {
        affirm("fields: a series", 0, science_last_error());
    }
    science_free(times), science_free(values), free(want_times), free(want_values);

    const char *which[] = {"grid", "blend"};
    for (int w = 0; w < 2; w++) {
        char path[64], *at = NULL;
        double bbox[4], stride = 1;
        size_t four = 0;
        snprintf(path, sizeof path, "%s.bbox", which[w]);
        science_json_numbers(ref, path, bbox, 4, &four);
        snprintf(path, sizeof path, "%s.stride", which[w]);
        science_json_number(ref, path, &stride);
        snprintf(path, sizeof path, "%s.time", which[w]);
        science_json_text(ref, path, &at);
        snprintf(path, sizeof path, "%s.values", which[w]);
        want = numbers(ref, path, &k);
        science_grid g = {0};
        code = science_grid_read(api, "gfs-forecast", "tmp2m", at, bbox, (int)stride, &g);
        if (code == SCIENCE_OK && g.nlon * g.nlat == k) {
            double *as_double = malloc(k * sizeof *as_double);
            for (size_t i = 0; i < k; i++) as_double[i] = g.values[i];
            check(w ? "fields: a grid between stored times equals the Python client's, as float32"
                    : "fields: a grid across the antimeridian equals the Python client's, as float32",
                  difference(as_double, want, k), 0, "K");
            free(as_double);
        } else {
            affirm(w ? "fields: a grid between stored times" : "fields: a grid across the antimeridian", 0, science_last_error());
        }
        science_grid_free(&g);
        science_free(at);
        free(want);
    }

    char *answer = NULL;
    code = science_get(api, "/fields/no-such-field", &answer);
    affirm("errors: a missing field answers SCIENCE_ERROR_HTTP with status 404", code == SCIENCE_ERROR_HTTP && science_last_status() == 404, science_last_error());
    science_free(answer);
}

static void compute(science_client *api, const char *ref) {
    double want = 0, got = 0;
    if (science_json_number(ref, "compute.0.result.ra", &want) != SCIENCE_OK) {
        skip("compute: the Python client could not call the domains");
        return;
    }
    char *result = NULL;
    int code = science_compute(api, "astronomy", "sunPosition", "[2460000.5]", &result);
    double worst = code == SCIENCE_OK ? 0 : INFINITY;
    const char *keys[] = {"ra", "dec", "longitude", "distanceAu"};
    for (int k = 0; k < 4 && code == SCIENCE_OK; k++) {
        char path[64];
        snprintf(path, sizeof path, "compute.0.result.%s", keys[k]);
        science_json_number(ref, path, &want);
        if (science_json_number(result, keys[k], &got) != SCIENCE_OK) worst = INFINITY;
        else if (fabs(got - want) > worst) worst = fabs(got - want);
    }
    check("compute: the Sun's position equals the Python client's", worst, 0, "degree");
    science_free(result);
    code = science_compute(api, "meteorology", "windFromDirection", "[5, 0]", &result);
    affirm("compute: a wind from the west comes from 270 degrees", code == SCIENCE_OK && strcmp(result, "270") == 0, result ? result : science_last_error());
    science_free(result);
    code = science_compute(api, "astronomy", "sunPosition", "[2460000.5,", &result);
    affirm("compute: arguments that are not JSON are an argument error", code == SCIENCE_ERROR_ARGUMENT, science_last_error());
    code = science_compute(api, "astronomy", "noSuchFunction", "[]", &result);
    affirm("compute: a function that is not there answers 404", code == SCIENCE_ERROR_HTTP && science_last_status() == 404, science_last_error());
}

static void observations(science_client *api, const char *source) {
    const char *variables[] = {"air_temperature", "relative_humidity"};
    int code = science_define(api, source, "API check (C)", variables, 2);
    affirm("observations: a source is defined", code == SCIENCE_OK, science_last_error());

    struct listening l = {api, source, 0, seconds_now() + 8, 0, "", "", -1};
    pthread_t thread;
    pthread_create(&thread, NULL, listen_thread, &l);
    while (!l.connected && seconds_now() < l.deadline) {
        struct timespec pause = {0, 20000000};
        nanosleep(&pause, NULL);
    }
    struct timespec settle = {0, 200000000};
    nanosleep(&settle, NULL);
    const char *names[] = {"air_temperature", "relative_humidity"};
    double values[] = {289.5, 58};
    double t0 = seconds_now();
    code = science_post_reading(api, source, "c1", 139.76, 35.68, NULL, names, values, 2);
    affirm("observations: a reading is posted", code == SCIENCE_OK, science_last_error());
    pthread_join(thread, NULL);
    if (l.received > 0) {
        check("observations: a posted reading reaches the stream", l.received - t0, 1.0, "s");
        affirm("observations: the stream carries the reading and its id", strcmp(l.station, "c1") == 0 && l.id[0] != '\0', l.station);
    } else {
        affirm("observations: a posted reading reaches the stream", 0, l.code ? science_last_error() : "nothing within 8 s");
    }

    char path[160], *answer = NULL;
    snprintf(path, sizeof path, "/observations/%s?limit=1", source);
    code = science_get(api, path, &answer);
    affirm("observations: the latest reading reads back", code == SCIENCE_OK && text_is(answer, "readings.0.station", "c1"), answer);
    science_free(answer);

    double hot = 1000;
    code = science_post_reading(api, source, "c1", 139.76, 35.68, NULL, names, &hot, 1);
    affirm("observations: an implausible reading is SCIENCE_ERROR_REFUSED, with its reason",
           code == SCIENCE_ERROR_REFUSED && strstr(science_last_error(), "air_temperature") != NULL, science_last_error());
    code = science_post(api, source, "{\"readings\": [{\"station\": \"c2\", \"time\": \"2026-10-23T06:00:00Z\", \"values\": {\"air_temperature\": 290}}]",
                        &answer);
    affirm("observations: a body that is not JSON is an argument error, sent nowhere", code == SCIENCE_ERROR_ARGUMENT, science_last_error());
}

static char *read_all(FILE *in) {
    size_t size = 0, cap = 1 << 16;
    char *buf = malloc(cap);
    size_t got;
    while (buf && (got = fread(buf + size, 1, cap - size - 1, in)) > 0) {
        size += got;
        if (cap - size < 2) buf = realloc(buf, cap *= 2);
    }
    if (buf) buf[size] = '\0';
    return buf;
}

int main(int argc, char **argv) {
    offline();
    if (argc < 2 || strcmp(argv[1], "--offline") == 0) return 0;
    char *ref = read_all(stdin);
    science_client *api = science_connect(argv[1], argc > 3 ? argv[3] : NULL, 60);
    if (!api) {
        affirm("a client of the service", 0, science_last_error());
        return 1;
    }
    char *info = NULL;
    int code = science_info(api, &info);
    char *version = NULL;
    science_json_text(ref, "version", &version);
    affirm("the service's index", code == SCIENCE_OK && text_is(info, "api", "v1") && version && text_is(info, "version", version), science_last_error());
    science_free(info), science_free(version);
    fields(api, ref);
    compute(api, ref);
    observations(api, argc > 2 ? argv[2] : "api-check-c");

    science_client *nowhere = science_connect("http://127.0.0.1:1", NULL, 2);
    code = science_info(nowhere, &info);
    affirm("errors: a service that is not there is SCIENCE_ERROR_UNREACHABLE, with status 0",
           code == SCIENCE_ERROR_UNREACHABLE && science_last_status() == 0, science_last_error());
    science_close(nowhere);
    science_close(api);
    free(ref);
    return 0;
}
