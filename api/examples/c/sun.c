/* A domain's functions from C: where the Sun stands, now, computed by the
 * astronomy domain (domains/astronomy/js/sun.js) through the service. The
 * arguments go as a JSON array, the result comes back as JSON text, and the
 * result of one call is an argument of the next:
 *
 *   julianDate(ms)              the Julian date of a UTC instant (ms since 1970);
 *   sunPosition(jd)             the Sun's apparent right ascension and declination
 *                               (degrees, equinox of date) and distance (au), from
 *                               the Astronomical Almanac's low-precision formulae (0.01°);
 *   subsolarPoint(ms)           where the Sun is at the zenith, {lon, lat};
 *   sunAltAz(lon, lat, sub)     its altitude and azimuth (from north through east) at a place;
 *   daylight(altitude)          day, a twilight, or night.
 *
 *   api/build/c-sun [lon lat] [url]
 */

#define _POSIX_C_SOURCE 200809L

#include "science.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* A call whose result is wanted; NULL, with the reason printed, when it fails. */
static char *call(science_client *api, const char *function, const char *args) {
    char *result = NULL;
    if (science_compute(api, "astronomy", function, args, &result) != SCIENCE_OK) {
        fprintf(stderr, "%s: %s\n", function, science_last_error());
        return NULL;
    }
    return result;
}

int main(int argc, char **argv) {
    double lon = argc > 2 ? atof(argv[1]) : 139.69, lat = argc > 2 ? atof(argv[2]) : 35.69;
    science_client *api = science_connect(argc > 3 ? argv[3] : NULL, NULL, 60);
    if (!api) {
        fprintf(stderr, "%s\n", science_last_error());
        return 1;
    }
    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    double ms = now.tv_sec * 1e3 + now.tv_nsec / 1e6; /* UTC, since 1970 */

    char args[512];
    snprintf(args, sizeof args, "[%.17g]", ms);
    char *jd = call(api, "julianDate", args);
    char *sub = call(api, "subsolarPoint", args);
    if (!jd || !sub) return 1;
    snprintf(args, sizeof args, "[%s]", jd);
    char *sun = call(api, "sunPosition", args);
    snprintf(args, sizeof args, "[%.17g, %.17g, %s]", lon, lat, sub); /* the subsolar point, as it came */
    char *sky = call(api, "sunAltAz", args);
    if (!sun || !sky) return 1;
    double ra, dec, au, sub_lon, sub_lat, altitude, azimuth;
    science_json_number(sun, "ra", &ra);
    science_json_number(sun, "dec", &dec);
    science_json_number(sun, "distanceAu", &au);
    science_json_number(sub, "lon", &sub_lon);
    science_json_number(sub, "lat", &sub_lat);
    science_json_number(sky, "altitude", &altitude);
    science_json_number(sky, "azimuth", &azimuth);
    snprintf(args, sizeof args, "[%.17g]", altitude);
    char *light = call(api, "daylight", args), *name = NULL;
    if (light) science_json_text(light, "", &name); /* the JSON string, unquoted */

    printf("JD %s (UTC)\n", jd);
    printf("  the Sun: right ascension %.3f°, declination %+.3f° (equinox of date), %.5f au away\n", ra, dec, au);
    printf("  overhead at %.2f°, %+.2f°\n", sub_lon, sub_lat);
    printf("  at %.2f°, %+.2f°: altitude %+.1f°, azimuth %.1f° (%s)\n", lon, lat, altitude, azimuth, name ? name : "?");
    science_free(jd), science_free(sub), science_free(sun), science_free(sky), science_free(light), science_free(name);
    science_close(api);
    return 0;
}
