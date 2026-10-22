/* Following the stream of readings from C: each reading posted to the service
 * is printed as it arrives. The callback receives each event's id, its source
 * and its readings as JSON text, read here with science_json_text and
 * science_json_number; it is called with NULLs at each heartbeat, and returns 0
 * to stop. After a dropped connection the library reconnects by itself and
 * receives what was missed.
 *
 *   api/build/c-stream [source, or "" for every source] [count, 0 for ever] [url]
 */

#include "science.h"

#include <stdio.h>
#include <stdlib.h>

struct seen {
    long count, limit;
};

static int on_event(const char *id, const char *source, const char *readings, void *user) {
    struct seen *s = user;
    if (!readings) return 1; /* a heartbeat: go on listening */
    size_t n = 0;
    science_json_length(readings, "", &n); /* the readings are an array */
    for (size_t k = 0; k < n; k++) {
        char path[64], *station = NULL, *time = NULL, *values = NULL;
        snprintf(path, sizeof path, "%zu.station", k);
        science_json_text(readings, path, &station);
        snprintf(path, sizeof path, "%zu.time", k);
        science_json_text(readings, path, &time);
        snprintf(path, sizeof path, "%zu.values", k);
        science_json_text(readings, path, &values); /* {variable: value}, in SI units */
        printf("[%s] %s %s at %s: %s\n", id, source, station ? station : "?", time ? time : "?", values ? values : "{}");
        science_free(station);
        science_free(time);
        science_free(values);
    }
    fflush(stdout);
    return s->limit == 0 || ++s->count < s->limit;
}

int main(int argc, char **argv) {
    const char *source = argc > 1 && argv[1][0] ? argv[1] : NULL;
    struct seen s = {0, argc > 2 ? atol(argv[2]) : 0};
    science_client *api = science_connect(argc > 3 ? argv[3] : NULL, NULL, 60);
    if (!api) {
        fprintf(stderr, "%s\n", science_last_error());
        return 1;
    }
    printf("Listening to %s\n", source ? source : "every source");
    int code = science_stream(api, source, NULL, on_event, &s);
    if (code != SCIENCE_OK) fprintf(stderr, "%s\n", science_last_error());
    science_close(api);
    return code == SCIENCE_OK ? 0 : 1;
}
