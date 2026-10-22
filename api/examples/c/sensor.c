/* A sensor's program for a device running Linux (a Raspberry Pi, a
 * BeagleBone, any small computer with a network), in C: the template to start
 * from. It defines its source once, then measures and posts a reading every
 * `--every` seconds, keeping the readings it could not send while the network
 * or the service was away and sending them, with their own times, when it
 * returns.
 *
 *   make -C api && api/build/c-sensor [--source garden] [--station garden-1] [--lon 139.69 --lat 35.69]
 *                                     [--every 60] [--count 0] [--w1 PATH] [--url URL] [--token T]
 *
 * What it measures. read_sensor() reads a DS18B20 temperature probe on the
 * 1-Wire bus when --w1 names its file (on a Raspberry Pi with dtoverlay=w1-gpio
 * in /boot/config.txt: /sys/bus/w1/devices/28-xxxxxxxxxxxx/temperature, which
 * holds millidegrees Celsius), and makes a value up otherwise, so that the
 * program runs anywhere. For a device of your own, change read_sensor(), the
 * variables of the source (each "standard_name:units", in the units the device
 * gives; the store converts them to SI) and the station's place.
 *
 * What it needs. The device's clock must be set (NTP), since a reading's time
 * is the device's: the store refuses times before 1990 or more than a quarter
 * of an hour ahead. A service on another computer is started with a token
 * (--token, or $SCIENCE_API_TOKEN), which the device then gives; it crosses the
 * network in the clear over plain HTTP, so keep both on a network you trust, or
 * reach the service through a TLS proxy or tunnel.
 *
 * Run as a service of the system, it starts with the device and restarts after
 * a failure: a unit /etc/systemd/system/science-sensor.service with
 *     [Service]  ExecStart=/home/pi/c-sensor --url http://192.168.1.10:8800 --token ...   Restart=always
 *     [Install]  WantedBy=multi-user.target
 * then `systemctl enable --now science-sensor`.
 */

#define _POSIX_C_SOURCE 200809L

#include "science.h"

#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define KEEP 1024 /* readings kept while the service is away: 17 hours at one a minute */

struct reading {
    char time[32];
    double celsius;
};

static volatile sig_atomic_t stopping = 0;
static void stop(int signal) { (void)signal; stopping = 1; }

/* The measurement, in degrees Celsius; NaN when the sensor gives none. */
static double read_sensor(const char *w1) {
    if (w1) {
        FILE *f = fopen(w1, "r");
        long millidegrees;
        int ok = f && fscanf(f, "%ld", &millidegrees) == 1;
        if (f) fclose(f);
        return ok ? millidegrees / 1000.0 : NAN;
    }
    /* No sensor: a daily cycle peaking in the afternoon (local solar time at 139.69 E), and some noise. */
    double hours = fmod(time(NULL) / 3600.0 + 139.69 / 15.0, 24.0);
    double made_up = 17.0 + 5.0 * sin((hours - 9.0) / 24.0 * 2.0 * acos(-1.0)) + 0.4 * (rand() / (double)RAND_MAX - 0.5);
    return round(made_up * 100) / 100; /* to the hundredth, as a probe reads */
}

static void utc_now(char *out, size_t size) {
    time_t now = time(NULL);
    struct tm utc;
    gmtime_r(&now, &utc);
    strftime(out, size, "%Y-%m-%dT%H:%M:%SZ", &utc);
}

static void pause_for(double seconds) {
    struct timespec t = {(time_t)seconds, (long)((seconds - floor(seconds)) * 1e9)};
    nanosleep(&t, NULL); /* a signal ends it early */
}

int main(int argc, char **argv) {
    const char *source = "garden", *station = "garden-1", *w1 = NULL, *url = NULL, *token = NULL;
    double lon = 139.69, lat = 35.69, every = 60;
    long count = 0;
    for (int i = 1; i + 1 < argc; i += 2) {
        if (!strcmp(argv[i], "--source")) source = argv[i + 1];
        else if (!strcmp(argv[i], "--station")) station = argv[i + 1];
        else if (!strcmp(argv[i], "--lon")) lon = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "--lat")) lat = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "--every")) every = atof(argv[i + 1]);
        else if (!strcmp(argv[i], "--count")) count = atol(argv[i + 1]);
        else if (!strcmp(argv[i], "--w1")) w1 = argv[i + 1];
        else if (!strcmp(argv[i], "--url")) url = argv[i + 1];
        else if (!strcmp(argv[i], "--token")) token = argv[i + 1];
        else {
            fprintf(stderr, "unknown option %s\n", argv[i]);
            return 2;
        }
    }
    signal(SIGINT, stop);
    signal(SIGTERM, stop);
    science_client *api = science_connect(url, token, 20);
    if (!api) {
        fprintf(stderr, "%s\n", science_last_error());
        return 1;
    }

    /* A clock not yet set by NTP reads 1970: wait for it, rather than send readings the store refuses. */
    while (!stopping && time(NULL) < 1577836800) { /* 2020-01-01 */
        fprintf(stderr, "waiting for the clock to be set\n");
        pause_for(10);
    }

    const char *variables[] = {"air_temperature:degree_C"};
    const char *names[] = {"air_temperature"};
    struct reading *pending = malloc(KEEP * sizeof *pending);
    size_t waiting = 0;
    int defined = 0;
    for (long taken = 0; !stopping && (count == 0 || taken < count); taken++) {
        double celsius = read_sensor(w1);
        if (!isnan(celsius)) {
            if (waiting == KEEP) { /* full: the oldest is let go */
                memmove(pending, pending + 1, (KEEP - 1) * sizeof *pending);
                waiting--;
            }
            utc_now(pending[waiting].time, sizeof pending[waiting].time);
            pending[waiting++].celsius = celsius;
        } else {
            fprintf(stderr, "the sensor gave no value\n");
        }

        if (!defined) { /* once, and again after a failure to reach the service */
            int code = science_define(api, source, "Garden station", variables, 1);
            defined = code == SCIENCE_OK;
            if (!defined) fprintf(stderr, "the source is not defined yet: %s\n", science_last_error());
        }
        size_t sent = 0;
        while (defined && sent < waiting) {
            struct reading *r = &pending[sent];
            int code = science_post_reading(api, source, station, lon, lat, r->time, names, &r->celsius, 1);
            if (code == SCIENCE_OK || code == SCIENCE_ERROR_REFUSED) {
                if (code == SCIENCE_OK) printf("%s  %.2f °C  sent\n", r->time, r->celsius);
                else fprintf(stderr, "%s  %.2f °C  refused: %s\n", r->time, r->celsius, science_last_error());
                sent++;
                continue;
            }
            if (code == SCIENCE_ERROR_HTTP && (science_last_status() == 401 || science_last_status() == 403)) {
                fprintf(stderr, "%s\n", science_last_error()); /* no token, or the wrong one: waiting will not help */
                stopping = 1;
            } else {
                fprintf(stderr, "the service is away (%s); %zu readings kept\n", science_last_error(), waiting - sent);
                if (code == SCIENCE_ERROR_HTTP && science_last_status() == 404) defined = 0; /* the source is gone: define it again */
            }
            break;
        }
        memmove(pending, pending + sent, (waiting - sent) * sizeof *pending);
        waiting -= sent;
        fflush(stdout);
        if (!stopping && (count == 0 || taken + 1 < count)) pause_for(every);
    }
    free(pending);
    science_close(api);
    return 0;
}
