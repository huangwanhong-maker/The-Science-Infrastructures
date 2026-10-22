/* The forecast at a place over time, from C: GFS's temperature at 2 m at every
 * time of the product, as a series (the value at the place, bilinear between
 * the four grid points around it), and the first time from now on.
 *
 *   make -C api && api/build/c-forecast [lon lat] [url]     (the service: $SCIENCE_API_URL or http://127.0.0.1:8800)
 */

#define _POSIX_C_SOURCE 200809L

#include "science.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv) {
    double lon = argc > 2 ? atof(argv[1]) : 139.69; /* degrees east, WGS 84 (Tokyo by default) */
    double lat = argc > 2 ? atof(argv[2]) : 35.69;
    science_client *api = science_connect(argc > 3 ? argv[3] : NULL, NULL, 60);
    if (!api) {
        fprintf(stderr, "%s\n", science_last_error());
        return 1;
    }
    size_t count = 0;
    double *times = NULL, *kelvin = NULL; /* times: seconds since 1970-01-01T00:00:00Z */
    if (science_series(api, "gfs-forecast", "tmp2m", lon, lat, &count, &times, &kelvin) != SCIENCE_OK) {
        fprintf(stderr, "%s\n", science_last_error());
        science_close(api);
        return 1;
    }
    time_t now = time(NULL);
    printf("GFS at %.2f, %.2f (times UTC)\n", lon, lat);
    int marked = 0;
    for (size_t k = 0; k < count; k++) {
        time_t t = (time_t)times[k];
        struct tm utc;
        char text[32];
        gmtime_r(&t, &utc);
        strftime(text, sizeof text, "%Y-%m-%d %H:%M", &utc);
        int first_ahead = !marked && t >= now;
        marked |= first_ahead;
        printf("  %s  %6.1f °C%s\n", text, kelvin[k] - 273.15, first_ahead ? "   <- the next time" : "");
    }
    science_free(times);
    science_free(kelvin);
    science_close(api);
    return 0;
}
