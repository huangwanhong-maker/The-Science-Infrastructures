// Guided tours through the brain: a few steps each, with the structures shown, cut and
// coloured as the text needs. A step names what to show (`show`, with the gross
// structures always there), the cortex's style and colouring, a cut, a selection, and
// a view: the camera's target, its direction from the target, and its distance (mm).
//
// The statements follow the knowledge base of the neuroscience domain and its sources
// (Kandel et al. 2021; Standring 2020; Vanderah and Gould 2021; Catani and Thiebaut de
// Schotten 2008; Papez 1937; Scoville and Milner 1957).

const CENTRE = [0, -16, 12];
const LEFT = [-1, 0.12, 0.18], RIGHT = [1, 0.12, 0.18], TOP = [0.0001, -0.02, 1], FRONT_LEFT = [-0.78, 0.42, 0.46], BELOW = [-0.2, 0.15, -1];

export const TOURS = [
  {
    id: 'surface',
    title: 'The face of the brain',
    summary: 'The lobes, the two great furrows that bound them, and the gyri a neurologist names first.',
    steps: [
      { title: 'Four lobes and a hidden fifth', text: 'Each hemisphere is divided by convention into lobes named for the skull bones above them: frontal, parietal, temporal and occipital. A fifth, the insula, lies hidden in the depth of the lateral sulcus, and the limbic lobe rings the medial wall.', colouring: 'lobes', view: { target: CENTRE, direction: FRONT_LEFT, distance: 560 } },
      { title: 'The central sulcus', text: 'The central sulcus runs from the top of the hemisphere down and forward almost to the lateral sulcus. In front of it, the precentral gyrus is the primary motor cortex; behind it, the postcentral gyrus is the primary somatosensory cortex. Each holds a map of the body, the leg at the top and the face below.', select: 'precentral-gyrus', side: 'left', view: { target: [-30, -18, 40], direction: LEFT, distance: 470 } },
      { title: 'The lateral sulcus', text: 'The lateral (Sylvian) sulcus separates the temporal lobe below from the frontal and parietal lobes above. Its banks are the opercula, the lids over the insula; on its lower bank lie the planum temporale and Heschl’s gyrus, the primary auditory cortex.', select: 'superior-temporal-gyrus', side: 'left', view: { target: [-40, -10, 6], direction: LEFT, distance: 440 } },
      { title: 'The insula, uncovered', text: 'Cut away the cortex at the level of the lateral sulcus and the insula appears: a triangle of short and long gyri ringed by the circular sulcus. It reads the state of the body, from the heartbeat to taste and pain, and is part of the brain’s system for what matters now.', cortex: 'solid', select: 'insula', side: 'left', cut: { axis: 'sagittal', position: -38, keep: 1, cap: true }, view: { target: [-38, -4, 4], direction: LEFT, distance: 420 } },
      { title: 'The medial wall', text: 'A cut between the hemispheres shows the medial surface: the cingulate gyrus arching over the corpus callosum, the precuneus behind it, and at the back the calcarine sulcus, whose banks hold the primary visual cortex, with the cuneus above it and the lingual gyrus below.', select: 'cingulate-gyrus', cut: { axis: 'sagittal', position: 2, keep: 1, cap: true }, view: { target: [2, -20, 14], direction: LEFT, distance: 520 } },
    ],
  },
  {
    id: 'deep',
    title: 'Into the deep nuclei',
    summary: 'Through a glass cortex to the basal ganglia, the thalamus, the hippocampus and the ventricles.',
    steps: [
      { title: 'The cortex as glass', text: 'With the cortex made glass, the grey nuclei buried in the white matter come into view. They are not cortex but clusters of cell bodies, each with its own connections, and they are where many diseases of movement and memory begin.', cortex: 'glass', view: { target: [0, -10, 8], direction: FRONT_LEFT, distance: 470 } },
      { title: 'The striatum', text: 'The caudate nucleus curves along the wall of the lateral ventricle; the putamen lies lateral to it, the two joined by bridges of grey matter through the internal capsule, whence the name striatum. They take in nearly the whole cortex and loop back to it through the pallidum and the thalamus.', cortex: 'glass', select: 'caudate-nucleus', side: 'left', show: ['putamen', 'caudate-nucleus'], view: { target: [-14, 2, 8], direction: FRONT_LEFT, distance: 330 } },
      { title: 'A coronal section', text: 'A coronal cut near the anterior commissure shows the order of things from the midline outward: ventricle, caudate, internal capsule, globus pallidus, putamen, external capsule, claustrum, extreme capsule, insula. The face of the cut is the template’s image, washed in the colour of each structure.', cortex: 'solid', cut: { axis: 'coronal', position: 2, keep: -1, cap: true }, select: 'putamen', side: 'left', view: { target: [0, 2, 8], direction: [-0.12, 1, 0.1], distance: 480 } },
      { title: 'The thalamus', text: 'Between the hemispheres, the thalamus is the gateway to the cortex: nearly every pathway to the cortex but smell relays here first. Its nuclei sort what passes, the lateral geniculate for vision and the medial for hearing, the ventral posterior for touch, the ventral lateral for movement.', cortex: 'glass', select: 'thalamus', show: ['thalamus'], view: { target: [-2, -18, 6], direction: FRONT_LEFT, distance: 300 } },
      { title: 'The hippocampus and the amygdala', text: 'In the medial temporal lobe, the hippocampus curls along the floor of the temporal horn of the ventricle, with the almond-shaped amygdala just in front. The first is needed to lay down new memories of events, the second to weigh what is threatening or rewarding.', cortex: 'glass', select: 'hippocampus', side: 'left', view: { target: [-24, -18, -12], direction: [-0.75, 0.35, 0.56], distance: 320 } },
      { title: 'The ventricles', text: 'The ventricles are the fluid-filled cavities within: two lateral ventricles, one in each hemisphere, joined to the third between the thalami, which drains by the cerebral aqueduct to the fourth in front of the cerebellum. The choroid plexus lining them makes about half a litre of cerebrospinal fluid a day.', cortex: 'glass', select: { system: 'csf' }, view: { target: [0, -14, 6], direction: FRONT_LEFT, distance: 420 } },
    ],
  },
  {
    id: 'language',
    title: 'The paths of language',
    summary: 'Broca’s and Wernicke’s regions, and the bundles of fibres that join them.',
    steps: [
      { title: 'Broca’s region', text: 'In most people language is lateralised to the left hemisphere. In the inferior frontal gyrus, the pars opercularis and pars triangularis form Broca’s region; damage to it and the cortex around it leaves speech halting and effortful, while the comprehension of simple sentences is relatively spared.', select: 'pars-opercularis', side: 'left', view: { target: [-46, 18, 8], direction: LEFT, distance: 430 } },
      { title: 'Wernicke’s region', text: 'Behind the auditory cortex, the posterior superior temporal gyrus and the planum temporale hold Wernicke’s region. Damage here leaves speech fluent but empty of sense, and comprehension poor. The left planum temporale is usually larger than the right.', select: 'superior-temporal-gyrus-posterior', side: 'left', view: { target: [-52, -34, 8], direction: LEFT, distance: 430 } },
      { title: 'The arcuate fasciculus', text: 'A great arc of fibres curves around the end of the lateral sulcus, joining the temporal and frontal language regions: the arcuate fasciculus. Tractography divides it into a direct long segment and an indirect pathway of anterior and posterior segments through the inferior parietal lobule (Catani et al. 2005).', cortex: 'glass', bundles: true, show: ['arcuate-fasciculus'], select: 'arcuate-fasciculus', side: 'left', view: { target: [-38, -18, 14], direction: LEFT, distance: 440 } },
      { title: 'The ventral stream', text: 'A second, ventral route runs through the temporal lobe: the inferior fronto-occipital fasciculus and the uncinate fasciculus join the temporal cortex to the frontal lobe. Stimulating them during surgery disturbs the meaning of words rather than their sound (Duffau et al. 2005).', cortex: 'glass', bundles: true, show: ['inferior-fronto-occipital-fasciculus', 'uncinate-fasciculus'], select: 'inferior-fronto-occipital-fasciculus', side: 'left', view: { target: [-30, 0, -2], direction: LEFT, distance: 470 } },
      { title: 'The network', text: 'Together these make the language network: the system card names its parts in this atlas. Each part has a card of its own; choose one to read it.', cortex: 'glass', bundles: true, select: { system: 'language' }, view: { target: [-30, -10, 8], direction: LEFT, distance: 480 } },
    ],
  },
  {
    id: 'memory',
    title: 'The circuit of memory',
    summary: 'From the hippocampus by the fornix and the thalamus to the cingulate gyrus: Papez’s circuit.',
    steps: [
      { title: 'Patient H.M.', text: 'In 1953 both medial temporal lobes of a young man with intractable epilepsy were removed. He lost the ability to form new memories of events and facts while keeping his intelligence, most of his older memories and the ability to learn new skills (Scoville and Milner 1957). The hippocampus became the key to memory.', cortex: 'glass', select: 'hippocampus', view: { target: [0, -20, -6], direction: FRONT_LEFT, distance: 400 } },
      { title: 'Inside the hippocampus', text: 'The hippocampus is a rolled sheet of cortex: the dentate gyrus, the fields CA3 to CA1 of the cornu Ammonis, and the subiculum. Signals from the entorhinal cortex pass through them in order, the trisynaptic circuit, and return to the cortex.', cortex: 'glass', show: ['ca1', 'ca23', 'dentate-gyrus', 'subiculum', 'presubiculum', 'uncus'], hide: ['amygdala'], select: 'ca1', side: 'left', view: { target: [-26, -22, -12], direction: [-0.6, 0.2, 0.78], distance: 260 } },
      { title: 'The fornix', text: 'The output of the hippocampus gathers into the fornix, which arches up under the corpus callosum and down again in front of the thalamus to the mammillary bodies of the hypothalamus.', cortex: 'glass', select: 'fornix', show: ['fornix'], view: { target: [0, -14, 6], direction: [-0.85, 0.2, 0.5], distance: 340 } },
      { title: 'Papez’s circuit', text: 'From the mammillary bodies the path runs to the anterior thalamus, on to the cingulate gyrus, by the cingulum to the parahippocampal gyrus, and back to the hippocampus. Papez proposed it in 1937 as a circuit of emotion; it is now known chiefly for its part in memory.', cortex: 'glass', bundles: true, select: { system: 'memory' }, view: { target: [0, -18, 10], direction: LEFT, distance: 480 } },
    ],
  },
  {
    id: 'vision',
    title: 'The way of sight',
    summary: 'From the optic chiasm through the lateral geniculate nucleus to the visual cortex.',
    steps: [
      { title: 'The optic chiasm', text: 'The optic nerves meet beneath the hypothalamus in the chiasm, where the fibres from the nasal half of each retina cross. Each side of the brain thereby receives the opposite half of the visual field from both eyes; a pituitary tumour pressing up on the chiasm takes away the outer half of each eye’s field.', cortex: 'glass', select: 'optic-chiasm', view: { target: [0, 2, -12], direction: BELOW, distance: 360 } },
      { title: 'The lateral geniculate nucleus', text: 'The optic tracts end in the lateral geniculate nucleus of the thalamus, six layers that keep the two eyes apart. From it the optic radiation fans back through the temporal and parietal white matter to the occipital lobe.', cortex: 'glass', show: ['lateral-geniculate-nucleus'], select: 'lateral-geniculate-nucleus', view: { target: [-20, -24, -4], direction: [-0.7, 0.3, 0.65], distance: 300 } },
      { title: 'The primary visual cortex', text: 'The primary visual cortex lines the calcarine sulcus on the medial occipital lobe. The upper visual field is mapped below the sulcus, in the lingual gyrus, and the lower field above it, in the cuneus; the centre of gaze takes the large share at the occipital pole.', select: 'intracalcarine-cortex', cut: { axis: 'sagittal', position: 2, keep: 1, cap: true }, view: { target: [4, -78, 6], direction: LEFT, distance: 380 } },
      { title: 'Two streams', text: 'Beyond the primary cortex, a ventral stream through the occipital and temporal lobes recognises what is seen, faces in the fusiform gyrus among them, and a dorsal stream into the parietal lobe guides action toward where things are (Goodale and Milner 1992).', select: { system: 'visual' }, view: { target: [0, -60, 0], direction: [-0.6, -0.55, -0.55], distance: 520 } },
    ],
  },
];
