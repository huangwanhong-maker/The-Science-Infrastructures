// The colours of the populations: eight hues in a fixed order, Dusty Orchid first, validated
// on the page's plum surface (#1d1520) with the six checks of the dataviz method (OKLCH
// lightness 0.48-0.67 and chroma >= 0.10; adjacent pairs CVD Delta E >= 8 under protanopia and
// deuteranopia, worst 8.3; normal-vision Delta E >= 15, worst 22.3; contrast >= 3:1). A
// population keeps its colour by its place in the specification; a ninth takes the first again
// only in the editor's swatch, never in a chart (charts fold the rest into "other").
export const POPULATION_COLOURS = ['#bf68ae', '#6ba04b', '#5f83d4', '#bb881a', '#8f68cb', '#cb7229', '#009d82', '#d56054'];

export const colourOf = (index) => POPULATION_COLOURS[index % POPULATION_COLOURS.length];
