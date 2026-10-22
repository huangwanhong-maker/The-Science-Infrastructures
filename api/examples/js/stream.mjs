// Readings posted and streamed from JavaScript: this script listens to a sensor
// source's stream, posts three readings to it as a device would, prints them as
// they come back, and stops. In a browser, the same Client streams through
// EventSource.
//
//   node api/examples/js/stream.mjs [source]       (the service: $SCIENCE_API_URL; writes may need $SCIENCE_API_TOKEN)

import { Client } from '../../js/science.mjs';

const source = process.argv[2] || 'example-garden';
const api = new Client();
const control = new AbortController();
let seen = 0;
const listening = (async () => {
  for await (const event of api.stream(source, { signal: control.signal })) {
    for (const r of event.readings) console.log(`received  ${r.time}  ${r.station}  ${JSON.stringify(r.values)}`);
    if ((seen += event.readings.length) >= 3) control.abort();
  }
})();

await new Promise((r) => setTimeout(r, 500));
const define = { title: 'Example garden station', variables: ['air_temperature:degree_C', 'relative_humidity:%'] };
for (let i = 0; i < 3; i++) {
  const reading = { station: 'garden-js', lon: 139.69, lat: 35.69, time: new Date(), values: { air_temperature: 18 + i / 10, relative_humidity: 60 } };
  const answer = await api.postObservations(source, [reading], { define: i === 0 ? define : undefined });
  console.log(`posted    accepted ${answer.accepted}`);
  await new Promise((r) => setTimeout(r, 300));
}
await listening;
