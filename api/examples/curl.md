# The API from curl

Every request of the API from a shell, against the standalone service (`python api/service/server.py`, port 8800). On Windows, run them in Git Bash, or replace the single quotes around JSON bodies by double quotes with the inner ones escaped. [../README.md](../README.md) explains the concepts; [../openapi.json](../openapi.json) describes every parameter.

## 1. The Service and Its Catalogue

```
curl http://127.0.0.1:8800/api/v1/                     # the version, the endpoints, the write policy
curl http://127.0.0.1:8800/api/v1/catalog              # the registry's datasets and the index of products
curl http://127.0.0.1:8800/api/v1/openapi.json         # the description
```

## 2. Fields

```
B=http://127.0.0.1:8800/api/v1
curl "$B/fields"                                                                   # the field products
curl "$B/fields/gfs-forecast"                                                      # grid, times, variables
curl "$B/fields/gfs-forecast/sample?var=tmp2m&lon=139.69,-0.13&lat=35.69,51.51&time=now+1d"
curl "$B/fields/gfs-forecast/series?var=tmp2m&lon=139.69&lat=35.69"
curl "$B/fields/gfs-forecast/grid?var=tmp2m&time=0&bbox=122,24,150,46&stride=2"     # JSON, rows north to south
curl -D - -o grid.f32 "$B/fields/gfs-forecast/grid?var=tmp2m&time=0&format=f32"     # float32; the grid in X-Grid
curl -OJ "$B/fields/gfs-forecast/export?var=tmp2m&time=now&format=tif&bbox=122,24,150,46"   # a GeoTIFF, named
curl -OJ "$B/fields/gfs-forecast/export?var=tmp2m&time=all&format=nc"                # every time, as netCDF
curl -X POST "$B/fields/gfs-forecast/sample" -H "Content-Type: application/json" \
     -d '{"var": "tmp2m", "lon": [139.69, 135.5], "lat": [35.69, 34.69], "time": "now"}'
```

## 3. Observations and Their Stream

```
curl "$B/observations"                                                     # the sets and the sensor sources
curl "$B/observations/metar?bbox=139,35,140.5,36.5&limit=20"               # aerodromes around Tokyo
curl "$B/observations/metar?bbox=139,35,140.5,36.5&format=geojson"
curl -N "$B/stream"                                                        # every new reading, as it arrives
curl -N "$B/observations/garden/stream"                                    # one source's
curl -X POST "$B/observations/garden" -H "Content-Type: application/json" -d '{
  "define": {"title": "Garden station", "variables": ["air_temperature:degree_C"]},
  "readings": [{"station": "g1", "lon": 139.69, "lat": 35.69, "time": "2026-10-23T06:00:00Z", "values": {"air_temperature": 17.5}}]}'
curl -X POST "$B/observations/garden" -H "Authorization: Bearer $SCIENCE_API_TOKEN" \
     -H "Content-Type: application/json" -d '[{"station": "g1", "time": "2026-10-23T06:10:00Z", "values": {"air_temperature": 17.8}}]'
```

The second post is for a service started with a token; a station already described needs no place.

## 4. Events, Tracks, Regions and Verification

```
curl "$B/events"
curl "$B/events/usgs-earthquakes-week?bbox=120,20,150,50&limit=10"
curl "$B/events/gcmt?since=2026-01-01&format=geojson"
curl "$B/tracks/gdacs-month?format=geojson"
curl "$B/regions/countries/at?lon=139.69,2.35&lat=35.69,48.86"
curl "$B/verification/earth-atlas"
```

## 5. The Domains' Functions

```
curl "$B/compute"                                                          # every unit, function and constant
curl "$B/compute/astronomy"
curl -X POST "$B/compute/astronomy/sunPosition" -H "Content-Type: application/json" -d '{"args": [2461000.5]}'
curl -X POST "$B/compute/meteorology/windFromDirection" -H "Content-Type: application/json" -d '[-3, 4]'
curl -X POST "$B/compute/geodesy/distanceKm" -H "Content-Type: application/json" -d '[35.69, 139.69, 51.51, -0.13]'
```

## 6. Plug-ins

```
curl "$B/plugins"
curl "$B/plugins/describe?app=earth-atlas&file=plugins/vorticity.py"
curl -X POST "$B/plugins/run" -H "Content-Type: application/json" \
     -d '{"app": "earth-atlas", "file": "plugins/vorticity.py", "chosen": {"wind": "gfs-levels-forecast|500 hPa"}}'
curl "$B/jobs/<id>"                                                        # until "state" is "done"
curl "$B/fields/plugin-vorticity-py/sample?var=zeta&lon=139.69&lat=35.69"
```
