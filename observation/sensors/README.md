# sensors: Our Own Instruments

*Status: draft.* A module of [`observation/`](../README.md), an edge. Use: readings that can be relied on as far as their sensors can; Learning: how a home station is sited and reduced.

## 1. Purpose

This module gives the user's own sensors a home: a weather station in the garden, a thermometer on a balcony, a tide pole read by hand, a logger on a boat. Their readings are kept beside the public observations (METAR, NDBC), in the same form, so that the Earth Atlas draws them on the same globe and a plug-in compares them with a forecast in the same way.

A **source** is one set of sensors under one description: its title, its variables with their units, and its stations, each with a place. A source receives **readings**: one station's values at one instant. Readings arrive in three ways:
- posted over HTTP to the application's API, by a device on the network (§3.1, §3.2);
- read from a file that a logger writes, by `ingest.py` (§3.3);
- made up by `simulate.py`, the template a real device copies (§3.4).

The store keeps every reading, checked and in SI units, and builds from the last hours an observations product that the atlas draws (§4).

## 2. Readings and Units

A reading is one JSON object:

```json
{"station": "garden-1", "lon": 139.70, "lat": 35.68, "elevation": 40,
 "time": "2026-10-22T09:00:00Z",
 "values": {"air_temperature": 289.35, "relative_humidity": 62.0, "air_pressure": 100870}}
```

| Field | Meaning |
|---|---|
| `station` | The station's identifier: letters, digits, `_`, `.`, `-`, at most 64 |
| `lon`, `lat` | Degrees east and north, WGS 84. They may be left out once the station is known; a station that moves (a boat) gives them each time. |
| `elevation` | Optional: metres above mean sea level |
| `time` | ISO 8601, UTC unless an offset is given; or milliseconds since 1970-01-01 UTC |
| `values` | The variables measured, by name; a missing one is left out or `null` |

**Names.** A variable is named by its [CF standard name](https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html) where one exists, as the fields and the public observations are. The store knows these, and keeps each in its SI unit:

| Variable | SI unit | Typical sensor |
|---|---|---|
| `air_temperature`, `dew_point_temperature` | K | BME280, SHT31, DS18B20 in a screen |
| `relative_humidity` | % (with respect to water) | BME280, SHT31 |
| `air_pressure` (at the station), `air_pressure_at_sea_level` | Pa | BME280, BMP390 |
| `wind_speed`, `wind_speed_of_gust` | m s-1 | cup or sonic anemometer |
| `wind_from_direction` | degree, clockwise from true north | wind vane |
| `precipitation_amount` (since the last reading) | kg m-2 (1 mm of water) | tipping-bucket gauge |
| `visibility_in_air` | m | |
| `sea_surface_temperature` | K | a thermistor afloat |
| `sea_surface_wave_significant_height` | m | |
| `water_surface_height_above_reference_datum` | m (state the datum in the long name) | tide pole, pressure or radar gauge |
| `surface_downwelling_shortwave_flux_in_air` | W m-2 | pyranometer |
| `soil_temperature` | K | |
| `mass_concentration_of_pm2p5_ambient_aerosol_particles_in_air`, `..._pm10_...` | kg m-3 | PMS5003, SDS011 |
| `mole_fraction_of_carbon_dioxide_in_air` | 1 (mol mol-1) | SCD30, SCD41 |

Any other variable may be declared with units of its own (`battery_voltage:V`); it is kept as it comes.

**Units.** Values are stored in SI. A source may declare that its values arrive in other units, and they are converted on arrival where the conversion is affine and known: `degree_C` and `degree_F` to K; `hPa`, `mbar`, `kPa`, `inHg`, `mmHg` to Pa; `km h-1`, `kt`, `mph` to m s-1; `mm` and `in` of water to kg m-2; `ug m-3` to kg m-3; `ppm` and `ppb` to 1; a fraction to per cent. The conversions are exact by definition (BIPM SI Brochure, 2019; NIST SP 811, 2008), except the conventional inch and millimetre of mercury (NIST SP 811, B.8). `store.CONVERSIONS` lists them all.

**Checks.** A reading is refused, with the reason, when:
- a variable is not one of the source's;
- a value is not a finite number, or lies beyond the extremes ever observed (the gross-error limits of `data/products/observation_set.py`, from the WMO Archive of Weather and Climate Extremes);
- its place is not on the Earth;
- its time is before 1990 (a clock not set, or seconds sent where milliseconds are expected) or more than a quarter of an hour ahead.

## 3. Feeding In a Device

### 3.1 A Raspberry Pi with a BME280, over HTTP

The BME280 (Bosch Sensortec, datasheet BST-BME280-DS002) measures temperature (±1 K), relative humidity (±3 %) and the pressure at the station (±1 hPa absolute). On a Raspberry Pi it sits on I²C at address 0x76 or 0x77. This loop posts one reading a minute to the application's API, with the standard library's `urllib` (and two small packages for the sensor: `pip install smbus2 RPi.bme280`):

```python
import datetime as dt, json, time, urllib.request
import bme280, smbus2

URL = "http://192.168.1.10:8800/api/v1/observations/garden"
TOKEN = "..."                                  # the API's token for this source
bus, address = smbus2.SMBus(1), 0x76
calibration = bme280.load_calibration_params(bus, address)

while True:
    s = bme280.sample(bus, address, calibration)          # degrees Celsius, hPa, per cent
    reading = {"station": "garden-1", "lon": 139.70, "lat": 35.68, "elevation": 40,
               "time": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
               "values": {"air_temperature": s.temperature + 273.15,
                          "relative_humidity": s.humidity,
                          "air_pressure": s.pressure * 100}}
    request = urllib.request.Request(URL, data=json.dumps({"readings": [reading]}).encode(),
                                     headers={"Content-Type": "application/json",
                                              "Authorization": f"Bearer {TOKEN}"})
    try:
        urllib.request.urlopen(request, timeout=10).close()
    except OSError as error:                   # no network: try again with the next reading
        print("not sent:", error)
    time.sleep(60)
```

The source is defined once beforehand, unless the API defines it from the first post (§6):

```
python observation/sensors/ingest.py garden --define title="Garden station" \
    --variable air_temperature:K --variable relative_humidity:% --variable air_pressure:Pa
```

*Remark.* A barometer measures the pressure where it is. To compare it with the pressure at sea level that weather maps show, reduce it as an altimeter setting (QNH) is reduced, through the ICAO standard atmosphere (ICAO Doc 7488, 1993): p₀ = p (1 − L h / T₀)^(−g M / (R L)), with L = 0.0065 K m⁻¹, T₀ = 288.15 K and the exponent g M / (R L) = 5.25588; h is the barometer's height in metres. At 40 m this adds 4.8 hPa. Post the result as `air_pressure_at_sea_level`, or post the station's pressure and let a reader reduce it.

*Caveat.* A temperature sensor in the sun, or on the board of a warm computer, reads the sensor, not the air. WMO's practice (WMO-No. 8, 2018, volume I, chapter 2) is a ventilated screen, shaded and white, 1.25 to 2 m above short grass, away from walls. A BME280 heats itself by a few tenths of a kelvin if sampled continuously; sample it in forced mode once a minute.

### 3.2 An ESP32 or Arduino, Posting JSON

An ESP32 (Arduino core, with the Adafruit BME280 library) takes its time from NTP, so that each reading carries UTC, and posts the same JSON:

```cpp
#include <WiFi.h>
#include <HTTPClient.h>
#include <Adafruit_BME280.h>
#include <time.h>

const char* SSID = "...";
const char* PASSWORD = "...";
const char* URL = "http://192.168.1.10:8800/api/v1/observations/balcony";
const char* TOKEN = "...";
Adafruit_BME280 bme;

void setup() {
  WiFi.begin(SSID, PASSWORD);
  while (WiFi.status() != WL_CONNECTED) delay(500);
  configTime(0, 0, "pool.ntp.org");                  // UTC
  bme.begin(0x76);
}

void loop() {
  time_t now = time(nullptr);
  if (now > 1700000000) {                            // the clock has been set
    char when[21];
    strftime(when, sizeof when, "%Y-%m-%dT%H:%M:%SZ", gmtime(&now));
    String body = String("{\"readings\":[{\"station\":\"balcony-1\",\"lon\":139.70,\"lat\":35.68,\"time\":\"")
      + when + "\",\"values\":{\"air_temperature\":" + String(bme.readTemperature() + 273.15, 2)
      + ",\"relative_humidity\":" + String(bme.readHumidity(), 1)
      + ",\"air_pressure\":" + String(bme.readPressure(), 0) + "}}]}";   // K, per cent, Pa
    HTTPClient http;
    http.begin(URL);
    http.addHeader("Content-Type", "application/json");
    http.addHeader("Authorization", String("Bearer ") + TOKEN);
    http.POST(body);
    http.end();
  }
  delay(60000);
}
```

A device that cannot keep a clock may send `"time"` as milliseconds since 1970 from a real-time clock module; one with neither should not guess, since the store refuses times before 1990.

### 3.3 A Weather Station or Logger Writing CSV

A logger, or a weather station's software, that writes CSV with a header `station,lon,lat,time,<variable>,...` (optionally `name` and `elevation`) is read by `ingest.py`. With `--follow` it watches the file as it grows, ingests each complete new line, rebuilds the product at most once a minute, and starts again from the top if the file is rotated:

```
station,lon,lat,time,air_temperature,relative_humidity
garden-1,139.70,35.68,2026-10-22T09:00:00Z,16.2,62
```

```
python observation/sensors/ingest.py garden logger.csv --follow \
    --define title="Garden station" --variable air_temperature:degree_C --variable relative_humidity:%
```

JSON lines (`.jsonl`), one reading a line as in §2, are read the same way. A file in another layout (a station's own export) needs a few lines of Python that rename its columns, or a call to `store.append`.

### 3.4 The Simulated Sensor

`simulate.py` is a weather sensor in software. Its temperature has a diurnal cycle in local solar time (with the equation of time), a mean that falls with latitude and height and follows the seasons, and weather that wanders; its humidity falls by day and rises by night; its pressure carries the atmosphere's semidiurnal tide, S2, with maxima near 10:00 and 22:00 local solar time (Haurwitz 1956). It posts each reading to the API, keeping those the API could not take and sending them with the next, or, without `--url`, appends straight to the store:

```
python observation/sensors/simulate.py garden --lon 139.70 --lat 35.68 --elevation 40
python observation/sensors/simulate.py garden --lon 139.70 --lat 35.68 --every 10 --count 30
python observation/sensors/simulate.py garden --lon 139.70 --lat 35.68 \
    --url http://127.0.0.1:8800/api/v1/observations/garden --token $TOKEN
```

Its POST body is `{"title", "variables", "readings": [...]}`: the description travels with the readings, so that a server may define a new source from its first post.

## 4. How Readings Appear

Each source becomes an observations product, `data/cache/products/observations-sensor-<source>.json` (schema `observations/1`, [data/README.md §4](../../data/README.md)), registered in the index of products as `sensor-<source>`. It holds the readings of the 48 hours before the newest, one row per station and time, with the station's name and place and the values in SI units. `ingest.py` builds it after ingesting, `simulate.py` after each reading, and the API after each post.

- The Earth Atlas draws it as an observations layer, beside METAR and NDBC: each station a dot coloured by one variable, at its reading nearest the time shown.
- A plug-in reads it with `context.observations('sensor-<source>')`, for example to compare the readings with a forecast or to correct one.

## 5. Privacy

*Caveat.* A product carries its stations' positions, as given. Anyone who can open the atlas, or the folder of products, sees where the sensors are, and a sensor in a garden says where someone lives. Give a home station's position rounded on the device (to 0.01°, about 1 km), since the store keeps what it is given. Keep the API's tokens secret, and the API on the home network unless it is meant to be reached from outside.

## 6. The Store's Interface

`store.py` is called by the scripts here and by the application's API; every function is safe to call from several threads.

| Function | Meaning |
|---|---|
| `sources()` | Every source: name, title, variables (name, units, long name, standard name), stations, count, time of the last reading |
| `define(name, title, variables, stations=None)` | Creates or updates a source; variables as `{name, units}` or `"name:units"`, where units are those the readings arrive in |
| `append(name, readings, rejected=None)` | Checks and appends readings, records new or moved stations, tells the subscribers; returns the number accepted, and lists the refused with their reasons in `rejected` |
| `read(name, since, until, bbox, limit)` | Readings in order of time, within a time window, inside a box (west, south, east, north; across the antimeridian when west > east), the last `limit` |
| `subscribe(callback)` | `callback(name, readings)` after each append in this process; returns the function that ends it |
| `build_product(name, hours=48)` | Writes and registers `observations-sensor-<name>.json` |

Source names are lowercase, digits and single hyphens (`garden-station`), so that none can leave the store. The store is `data/cache/observations/`, one directory per source with `source.json` and `readings.jsonl`; the environment variable `SCIENCE_OBSERVATIONS` moves it.

## 7. Files

| File | Purpose |
|---|---|
| `store.py` | The store: descriptions, checks, conversions, readings, subscribers, the product |
| `ingest.py` | Readings from a CSV or JSON-lines file, once or following it as it grows; or a source defined |
| `simulate.py` | A simulated weather sensor, posting to the API or appending to the store |

Checks: `data/tests/check_observations.py` (run by `node data/tests/run-tests.mjs`) defines, appends, reads, subscribes and builds in a temporary store.

## 8. References

- Alduchov, O. A. and Eskridge, R. E. (1996), Improved Magnus form approximation of saturation vapor pressure, *Journal of Applied Meteorology* 35, 601–609.
- Bosch Sensortec, *BME280 Combined humidity and pressure sensor*, datasheet BST-BME280-DS002.
- Chapman, S. and Lindzen, R. S. (1970), *Atmospheric Tides*, Reidel.
- Duffie, J. A. and Beckman, W. A. (2013), *Solar Engineering of Thermal Processes*, 4th edition, Wiley (the equation of time, eq. 1.5.3, after Spencer 1971).
- Haurwitz, B. (1956), The geographical distribution of the solar semidiurnal pressure oscillation, *Meteorological Papers* 2(5), New York University.
- ICAO (1993), *Manual of the ICAO Standard Atmosphere*, Doc 7488, 3rd edition.
- WMO (2018), *Guide to Instruments and Methods of Observation*, WMO-No. 8, volume I.
