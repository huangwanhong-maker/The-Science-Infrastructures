/* A grid read as float32 from C, and its mean: GFS's temperature at 2 m over a
 * box at an instant, its mean over the points, and its mean over the area,
 * which weights each point by the cosine of its latitude, the share of the
 * sphere a cell of a regular longitude-latitude grid covers (cpp/grid_mean.cpp
 * says why).
 *
 *   api/build/c-grid_mean [variable] [west south east north] [time] [url]
 *   (by default tmp2m over Japan, 122 24 150 46, now)
 */

#include "science.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    const char *variable = argc > 1 ? argv[1] : "tmp2m";
    double box[4] = {122, 24, 150, 46}; /* west, south, east, north; degrees */
    if (argc > 5)
        for (int k = 0; k < 4; k++) box[k] = atof(argv[2 + k]);
    const char *time = argc > 6 ? argv[6] : "now";
    science_client *api = science_connect(argc > 7 ? argv[7] : NULL, NULL, 60);
    science_grid g;
    if (!api || science_grid_read(api, "gfs-forecast", variable, time, box, 1, &g) != SCIENCE_OK) {
        fprintf(stderr, "%s\n", science_last_error());
        science_close(api);
        return 1;
    }
    const double rad = acos(-1.0) / 180.0;
    double sum = 0, weighted = 0, weights = 0;
    size_t points = 0;
    for (size_t row = 0; row < g.nlat; row++) {
        double w = cos((g.lat0 - row * g.dlat) * rad); /* rows run north to south */
        for (size_t col = 0; col < g.nlon; col++) {
            float x = g.values[row * g.nlon + col];
            if (isnan(x)) continue;
            sum += x;
            weighted += w * x;
            weights += w;
            points++;
        }
    }
    printf("%s at %s: %zu x %zu points, %zu with a value\n", variable, g.time, g.nlon, g.nlat, points);
    if (points) {
        printf("  mean over the points  %.3f %s\n", sum / points, g.units);
        printf("  mean over the area    %.3f %s   (weights cos φ)\n", weighted / weights, g.units);
    }
    science_grid_free(&g);
    science_close(api);
    return 0;
}
