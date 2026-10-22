/* The client of the API in C (api/README.md §3.5, api/openapi.json): the C
 * interface of the C++ client (api/cpp/science.hpp), for programs in C and for
 * any language that calls C, such as a sensor's program on a device running
 * Linux. It is built as libscience.a or libscience.so by api/Makefile; a C
 * program includes this header and links the library with the C++ runtime:
 *
 *     cc -std=c99 -I api/c mine.c api/build/libscience.a -lstdc++ -lm -o mine
 *
 *     science_client *api = science_connect("http://127.0.0.1:8800", NULL, 60);
 *     double lon = 139.69, lat = 35.69, kelvin;
 *     if (science_sample(api, "gfs-forecast", "tmp2m", "now+1d", &lon, &lat, 1, &kelvin) != SCIENCE_OK)
 *         fprintf(stderr, "%s\n", science_last_error());
 *     science_close(api);
 *
 * The conventions:
 *
 * - Every call returns SCIENCE_OK (0) or an error code; science_last_error()
 *   says why, and science_last_status() gives the HTTP status of an answer that
 *   was an error (0 when none came). Both are kept per thread, as errno is.
 * - What the library allocates for the caller (JSON text, arrays of values, a
 *   grid's values) is the caller's, released with science_free (or
 *   science_grid_free); buffers the caller passes in (the values of a sample)
 *   are the caller's throughout.
 * - Text is UTF-8. JSON goes in and out as text: the arguments of a call, a
 *   body of readings, an answer.
 * - A client holds only its address, its token and its time limit: each call
 *   opens its own connection, so a client may be created before the network is
 *   up and outlives a network that comes and goes. One client may serve several
 *   threads at once.
 * - Positions are degrees on WGS 84; times are UTC, as ISO 8601 text
 *   (2026-10-23T06:00:00Z), "now" or "now+6h", or NULL for the service's
 *   default, now; values are in the units the product states, SI unless it
 *   says otherwise; NaN where there is no value.
 * - Plain http:// only: for HTTPS, put the service behind a TLS proxy on a
 *   trusted hop, or use the Python client.
 */

#ifndef SCIENCE_H
#define SCIENCE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SCIENCE_VERSION "0.1.0"

enum science_code {
    SCIENCE_OK = 0,
    SCIENCE_ERROR_ARGUMENT = 1,    /* a null pointer, a bad address, arguments that are not JSON */
    SCIENCE_ERROR_UNREACHABLE = 2, /* no answer came: the service could not be reached, was silent, or broke the connection */
    SCIENCE_ERROR_HTTP = 3,        /* the service answered an error: science_last_status() says which (400, 401, 404, 503, ...) */
    SCIENCE_ERROR_ANSWER = 4,      /* an answer this client cannot read */
    SCIENCE_ERROR_REFUSED = 5,     /* a reading the store refused, with its reason */
    SCIENCE_ERROR_FAILED = 6       /* anything else: memory, a file that cannot be written, a run that failed */
};

typedef struct science_client science_client;

/* A grid, or a box of it: nlat rows of nlon values, row-major, rows north to
 * south; the point (row, col) lies at lon0 + col dlon, lat0 - row dlat. */
typedef struct science_grid {
    double lon0, dlon, lat0, dlat;
    size_t nlon, nlat;
    float *values;     /* nlon * nlat, NaN where there is no value; released by science_grid_free */
    char time[32];     /* the instant, ISO 8601 UTC */
    char units[32];
} science_grid;

/* Called for each event of a stream with its id, its source and its readings
 * as JSON text ([{station, lon, lat, time, values: {...}}]); and at each
 * heartbeat (every 15 s) with all three NULL, so that a quiet stream can be
 * stopped. Returns nonzero to go on, 0 to stop. The texts live until it returns. */
typedef int (*science_event_callback)(const char *id, const char *source, const char *readings, void *user);

/* ---- the client ------------------------------------------------------------------------------ */

/* A client of the service at url (NULL: $SCIENCE_API_URL, else
 * http://127.0.0.1:8800; with or without /api/v1), with the token its writes
 * need (NULL: $SCIENCE_API_TOKEN), and the time limit in seconds on connecting
 * and on each wait for an answer's bytes (0: 60). No request is made here.
 * NULL when the address is not one this client can reach. */
science_client *science_connect(const char *url, const char *token, double timeout_seconds);
void science_close(science_client *client);

const char *science_last_error(void);
int science_last_status(void);
void science_free(void *memory);

/* ---- any request ------------------------------------------------------------------------------- */

/* A GET by its path under /api/v1, with its query: "/events/usgs-earthquakes-week?limit=5". */
int science_get(science_client *client, const char *path, char **json);
/* The version, the endpoints and the write policy. */
int science_info(science_client *client, char **json);

/* ---- fields ---------------------------------------------------------------------------------------- */

/* Values of a variable of a field at count places and an instant, into values (count doubles). */
int science_sample(science_client *client, const char *field, const char *variable, const char *time,
                   const double *lon, const double *lat, size_t count, double *values);

/* The value at a place at every time of the field: *count times, as seconds
 * since 1970-01-01T00:00:00Z (UTC), and values; both released by science_free. */
int science_series(science_client *client, const char *field, const char *variable, double lon, double lat,
                   size_t *count, double **times, double **values);

/* The variable on its grid at a stored time ("0", "1", ...) or an instant, in
 * the box bbox = {west, south, east, north} (NULL: the whole grid), every
 * stride-th point (0 or 1: every one). */
int science_grid_read(science_client *client, const char *field, const char *variable, const char *time,
                      const double *bbox, int stride, science_grid *grid);
void science_grid_free(science_grid *grid);

/* ---- observations ------------------------------------------------------------------------------------ */

/* Defines a sensor source: its title and its variables, as "name:units"
 * ("air_temperature:degree_C", "relative_humidity:%"); the store converts the
 * values to SI on arrival. Defining it again updates it: the variables named
 * are described anew, the others and the stations are kept. A write. */
int science_define(science_client *client, const char *source, const char *title, const char *const *variables, size_t count);

/* Posts one reading of a station: its place (NaN, NaN for a station already
 * described), its time (NULL: now, by this computer's clock, which must be set),
 * and count values of the variables named. SCIENCE_ERROR_REFUSED, with the
 * store's reason, when the reading is not plausible. A write. */
int science_post_reading(science_client *client, const char *source, const char *station, double lon, double lat,
                         const char *time, const char *const *names, const double *values, size_t count);

/* Posts a body {readings: [...], define?: {...}} as JSON text; the answer
 * {accepted, rejected: [{index, reason}]} as JSON text (NULL: not wanted). A write. */
int science_post(science_client *client, const char *source, const char *body, char **answer);

/* Follows new readings as they arrive, of one source (NULL: every source),
 * after the event last_event_id (NULL: from now), until the callback returns 0.
 * After a dropped connection it reconnects by itself, asking for what it missed. */
int science_stream(science_client *client, const char *source, const char *last_event_id,
                   science_event_callback callback, void *user);

/* ---- the domains' functions -------------------------------------------------------------------------- */

/* Calls a pure function of a domain or foundation with its arguments as a JSON
 * array ("[2460000.5]"); its result as JSON text ("{\"ra\": ...}"). */
int science_compute(science_client *client, const char *unit, const char *function, const char *args, char **result);

/* ---- reading JSON ------------------------------------------------------------------------------------ */

/* The value at a path of a JSON text, so that a C program reads an answer
 * without a JSON library of its own. A path is keys and indices joined by '.':
 * "ra", "grid.nlon", "0.values.air_temperature"; "" is the whole text; a key
 * with a '.' in it cannot be named. SCIENCE_ERROR_ARGUMENT when the text is not
 * JSON, or nothing of that kind is at the path. Each call reads the whole text. */

/* A number; null reads as NaN. */
int science_json_number(const char *json, const char *path, double *value);
/* An array of numbers: its first `capacity` into values (null as NaN), its length into *count. */
int science_json_numbers(const char *json, const char *path, double *values, size_t capacity, size_t *count);
/* A string as it is, or any other value as its JSON text; released by science_free. */
int science_json_text(const char *json, const char *path, char **text);
/* The elements of an array, or the members of an object; 0 for any other value. */
int science_json_length(const char *json, const char *path, size_t *length);

#ifdef __cplusplus
}
#endif

#endif /* SCIENCE_H */
