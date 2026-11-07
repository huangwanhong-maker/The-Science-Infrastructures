// Questions in plain words, put to the tables of a place. Each is a tool's form filled in, so that it can be taken further in the
// form. `draw` names what goes on the map first, by a part of its title. Their claims are the tool's, on the records named; what
// each does and does not show is said in its `read`.

export const QUESTIONS = [
  // ---- Melbourne ------------------------------------------------------------------------------------------------------------
  {
    where: 'melbourne', ask: 'Where did the walkers come back after the lockdowns?', table: 'melbourne-counters-month', tool: 'ols', filter: { compare: 1 }, each: ['counter'],
    roles: { y: 'walkers_a_day', x: ['after'] }, params: { logs: 'y', covariance: 'HC1' }, draw: 'after · estimate', title: 'Walkers in 2025 against 2019, a log ratio',
    read: 'Each counter on its own: its months of 2025 against its months of 2019. The estimate is a difference of logarithms, so that −0.2 is about 18% fewer walkers.',
  },
  {
    where: 'melbourne', ask: 'Where do the centre\'s uses gather?', table: 'melbourne-blocks', tool: 'specialisation', filter: { census_year: 2024 },
    roles: { columns: ['office', 'retail', 'homes', 'hotels', 'education', 'health', 'culture', 'workshops', 'parking'] }, params: { focus: 'office' }, draw: 'Location quotient',
    read: 'A use\'s share of a block\'s floor space over its share of the whole municipality, in the census of 2024; choose another use in the form. Cells the census withheld count as none.',
  },
  {
    where: 'melbourne', ask: 'Which blocks hold more businesses than their floor space explains?', table: 'melbourne-blocks', tool: 'ols', filter: { census_year: 2024 }, byRow: true,
    roles: { y: 'establishments', x: ['floor_total'] }, params: { logs: 'both', covariance: 'HC1' }, draw: 'Residuals',
    read: 'Logarithms of each on the other, across the blocks in 2024. A block standing high holds many small premises for its floor space, and a tower\'s block stands low.',
  },
  {
    where: 'melbourne', ask: 'Where do food and drink cluster?', table: 'melbourne-blocks', tool: 'moran', filter: { census_year: 2024 },
    roles: { value: 'food_and_drink', lon: 'place_lon', lat: 'place_lat' }, params: { k: 6, permutations: 499 }, draw: 'Clusters',
    read: 'Cafés, restaurants, pubs and takeaways counted in the census of 2024, each block against its six nearest. Hot and cold spots at 5%, uncorrected for testing every block.',
  },
  {
    where: 'melbourne', ask: 'Has the spread of businesses over the blocks changed?', table: 'melbourne-blocks', tool: 'inequality', each: ['census_year'],
    roles: { value: 'establishments' }, params: { epsilon: 1 }, draw: 'Gini',
    read: 'The Gini coefficient of the blocks\' counts of businesses, census by census. Blocks with no business count.',
  },
  {
    where: 'melbourne', ask: 'Does rain keep walkers away?', table: 'melbourne-counters-day', tool: 'panel', byRow: true,
    roles: { unit: 'counter', time: 'day', y: 'walkers', x: ['rain_mm', 'temperature_c', 'weekend', 'holiday'] }, params: { effects: 'one-way', random: false, logs: 'y' }, draw: 'unit\'s effect',
    read: 'Each counter against itself since October 2024. A coefficient is a change in the logarithm of the day\'s walkers. The weather is one reading for the whole city, so that its errors are shared across counters on a day.',
  },
  {
    where: 'melbourne', ask: 'Which streets does rain empty most?', table: 'melbourne-counters-day', tool: 'ols', each: ['counter'],
    roles: { y: 'walkers', x: ['rain_mm', 'temperature_c', 'weekend', 'holiday'] }, params: { logs: 'y', covariance: 'HC1' }, draw: 'rain_mm · estimate', title: 'Walkers for a millimetre of rain, a change in their logarithm',
    read: 'A regression for each counter on the same day\'s rain and temperature: the change in the logarithm of its walkers for a millimetre of rain.',
  },
  // ---- Freetown ------------------------------------------------------------------------------------------------------------------
  {
    where: 'freetown', ask: 'How concentrated is the loss from floods among the wards?', table: 'freetown-wards', tool: 'concentration',
    roles: { value: 'annual_loss', label: 'name' }, params: { k: 4 },
    read: 'The assessment\'s annual average loss, ward by ward: how much the largest few wards carry.',
  },
  {
    where: 'freetown', ask: 'Where does the loss from floods cluster?', table: 'freetown-wards', tool: 'moran',
    roles: { value: 'loss_per_km2', lon: 'place_lon', lat: 'place_lat' }, params: { k: 5, permutations: 499 }, draw: 'Clusters',
    read: 'The annual loss for each km² of a ward, each ward against its five nearest.',
  },
  // ---- anywhere ------------------------------------------------------------------------------------------------------------------
  {
    where: 'anywhere', ask: 'What here is brighter than its people explain?', table: 'cells-here', tool: 'ols', byRow: true,
    roles: { y: 'lights', x: ['people'] }, params: { logs: 'both', covariance: 'HC1' }, draw: 'Residuals',
    read: 'The lights of 2024 on the people of 2024, in logarithms, across the cells of the 4° square about the point looked at. Neighbours are tested beside the fit, since they rarely stand apart.',
  },
  {
    where: 'anywhere', ask: 'How unevenly are the lights spread here?', table: 'cells-here', tool: 'inequality',
    roles: { value: 'lights' }, params: { epsilon: 1 },
    read: 'The Gini coefficient and its family, of the radiance on each cell of the 4° square about the point looked at.',
  },
  {
    where: 'anywhere', ask: 'Where do people cluster here?', table: 'cells-here', tool: 'moran',
    roles: { value: 'people', lon: 'place_lon', lat: 'place_lat' }, params: { k: 6, permutations: 199 }, draw: 'Clusters',
    read: 'Each cell against its six nearest, over the 4° square about the point looked at.',
  },
  {
    where: 'anywhere', ask: 'Where has the light grown fastest?', table: 'cells-years', tool: 'ols', each: ['cell'],
    roles: { y: 'lights', x: ['year'] }, params: { logs: 'none', covariance: 'classical' }, draw: 'year · estimate', title: 'Radiance gained a year, 2015–2024',
    read: 'A straight line through ten years of radiance, 2015 to 2024, for each cell of the 4° square about the point looked at: its slope a year.',
  },
  // ---- the world ------------------------------------------------------------------------------------------------------------------
  {
    where: 'world', ask: 'Do poorer economies grow faster?', table: 'economies', tool: 'convergence',
    roles: { unit: 'iso3', time: 'year', value: 'gdp_per_capita' }, params: { from: 2000, to: 2019 }, draw: 'Growth',
    read: 'Growth a year from 2000 to 2019 against the income of 2000, in current US dollars.',
  },
  {
    where: 'world', ask: 'What does income leave of life expectancy?', table: 'economies', tool: 'ols', filter: { year: 2019 }, byRow: true,
    roles: { y: 'life_expectancy', x: ['gdp_per_capita'] }, params: { logs: 'x', covariance: 'HC1' }, draw: 'Residuals',
    read: 'Life expectancy in 2019 on the logarithm of income. A place standing high lives longer than its income predicts.',
  },
  {
    where: 'world', ask: 'How concentrated are each economy\'s partners in trade?', tool: 'trade-concentration', params: {}, draw: 'Partner',
    read: 'The Herfindahl–Hirschman index of an economy\'s exports across its partners, in the latest year of the OECD\'s balanced trade.',
  },
  {
    where: 'world', ask: 'What does distance leave of trade?', tool: 'trade-gravity', params: { method: 'ppml' }, draw: 'Exporter',
    read: 'Gravity with exporter and importer effects, by Poisson pseudo-maximum likelihood: what each economy sends beyond what its distances explain.',
  },
];

export function questionsHere(area) {
  const groups = [];
  const take = (title, where) => { const items = QUESTIONS.filter((q) => q.where === where); if (items.length) groups.push({ title, items }); };
  if (area === 'melbourne') take('Melbourne', 'melbourne');
  if (area === 'freetown') take('Freetown', 'freetown');
  take('About the point looked at', 'anywhere');
  take('The world', 'world');
  return groups;
}
