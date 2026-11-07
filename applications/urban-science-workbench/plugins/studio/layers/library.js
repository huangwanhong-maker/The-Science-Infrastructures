// A record of the data library on the map, through time: all a plug-in need do to draw what happened where and when is to ask ctx.st to show it.

import { SOURCES } from '../library.js';

export async function layer(ctx) {
  const source = SOURCES[ctx.spec.source];
  const shown = ctx.st.show(await source.read(), { ...source.style, credit: source.credit });
  ctx.st.timeline.play();
  return {
    dispose() { shown.remove(); },
  };
}
