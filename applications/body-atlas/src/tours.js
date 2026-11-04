// Guided tours: each a sequence of views of the body, with what to notice. A step names the systems it
// shows, the part it chooses (by FMA concept), the parts it isolates, where it looks from, and, in the
// sections, the plane and the vertebra at which it cuts. Levels are read from the model's own vertebrae.

const ORGANS = { heart: 'solid', respiratory: 'solid', digestive: 'solid', urinary: 'solid', genital: 'solid', endocrine: 'solid', lymphoid: 'solid', eye: 'solid', ear: 'solid', brain: 'solid' };
const ONLY = (styles) => ({ only: true, ...styles });

export const TOURS = [
  {
    key: 'layers',
    title: 'The body in layers',
    about: 'From the skin inward: the skeleton, the muscles, the organs, the vessels and nerves.',
    steps: [
      {
        title: 'The anatomical position',
        text: 'BodyParts3D\'s man stands about 1.73 m tall, upright and facing forward, his arms at his sides. Every direction in the atlas is his own: his left is on your right as he faces you, and the compass in the corner names his directions as you turn him.',
        systems: ONLY({ skin: 'solid' }), view: 'front', frame: 'body',
      },
      {
        title: 'The skeleton',
        text: 'Bone, cartilage and the ligaments of the joints: 264 of the atlas\'s parts. The vertebrae name the levels of the trunk; the sections and the cards read them from this skeleton.',
        systems: ONLY({ skin: 'glass', skeleton: 'solid' }), view: 'front', frame: 'body',
      },
      {
        title: 'The muscles',
        text: 'Over four hundred muscles and their heads and parts, each a mesh of its own. Seen from behind, the trapezius and latissimus dorsi spread over the back.',
        systems: ONLY({ skin: 'glass', skeleton: 'solid', muscles: 'solid' }), view: 'back', frame: 'body',
      },
      {
        title: 'The organs',
        text: 'With the muscles set aside: the brain, the heart and the bronchial trees in the chest, the liver drawn as its segments, the stomach and intestines, the kidneys behind them.',
        systems: ONLY({ skin: 'glass', skeleton: 'glass', ...ORGANS }), view: 'front', frame: 'trunk',
      },
      {
        title: 'Vessels and nerves',
        text: 'Some six hundred arteries and four hundred veins, and the nerves. Click any of them to read it; long ones have a calibre along their course on their card.',
        systems: ONLY({ skin: 'glass', skeleton: 'glass', heart: 'solid', arteries: 'solid', veins: 'solid', nerves: 'solid', brain: 'solid' }), view: 'front', frame: 'body',
      },
    ],
  },
  {
    key: 'heart',
    title: 'The heart and the great vessels',
    about: 'The pump in the middle mediastinum, the vessels that leave and enter it, and its own arteries.',
    steps: [
      {
        title: 'The heart in the chest',
        text: 'Two thirds of the heart lie to the left of the midline, its apex forward and down to the left. Its card measures where it lies against the vertebrae.',
        systems: ONLY({ skeleton: 'glass', heart: 'solid', arteries: 'solid', veins: 'solid', respiratory: 'glass' }), select: 'FMA7088', view: 'front', frame: 'select',
      },
      {
        title: 'The great vessels',
        text: 'The aorta rises from the left ventricle and arches back over the left main bronchus, beginning and ending its arch at the level of the sternal angle. The pulmonary trunk leaves the right ventricle for the lungs; the venae cavae bring the body\'s blood back to the right atrium.',
        systems: ONLY({ heart: 'solid', arteries: 'solid', veins: 'solid', skeleton: 'glass' }),
        isolate: ['FMA7088', 'FMA3768', 'FMA3736', 'FMA87217', 'FMA8612', 'FMA4720', 'FMA10951'], select: 'FMA3768', view: 'front-left', frame: 'isolate',
      },
      {
        title: 'Its own arteries',
        text: 'The coronary arteries leave the aorta just above its valve. Notice that branches named "right" lie left of the midline here: in a name, right and left may be a chamber\'s or an artery\'s, not the body\'s side. The atlas checks the body\'s sides only far from the midline, where they mean the body\'s.',
        systems: ONLY({ heart: 'glass', arteries: 'solid' }), isolate: ['FMA50039', 'FMA50040', 'FMA7088'], select: 'FMA50039', view: 'back-left', frame: 'isolate',
      },
      {
        title: 'The blood a ventricle holds',
        text: 'BodyParts3D models the cavities of the heart as solids, so that the volume on the card is the blood held at the one instant the model shows. The mitral valve\'s two leaflets close its inflow.',
        systems: ONLY({ heart: 'solid' }), isolate: ['FMA9466', 'FMA7235', 'FMA11359'], select: 'FMA9466', view: 'front-left', frame: 'isolate',
      },
      {
        title: 'A section through the chambers',
        text: 'Cut transversely at the seventh thoracic vertebra and seen from below, as a CT scan is: the subject\'s left is on the right of the plate, his front at the top.',
        mode: 'sections', plate: { plane: 'transverse', vertebra: 'T7' }, systems: ONLY({ skin: 'solid', skeleton: 'solid', muscles: 'solid', heart: 'solid', arteries: 'solid', veins: 'solid', respiratory: 'solid', digestive: 'solid', nerves: 'solid' }),
      },
    ],
  },
  {
    key: 'veins',
    title: 'The way back from the foot',
    about: 'The venous return of the right lower limb, from the great saphenous vein to the right atrium.',
    steps: [
      {
        title: 'The great saphenous vein',
        text: 'The longest vein of the body begins at the dorsal venous arch of the foot, passes in front of the medial malleolus and climbs the inner side of the leg and thigh to the groin. Press Calibre on its card to see its width along its length.',
        systems: ONLY({ veins: 'solid', skeleton: 'glass', skin: 'glass' }), select: 'FMA21379', view: 'front-right', frame: 'select',
      },
      {
        title: 'The deep veins and the calf\'s pump',
        text: 'Deep in the calf the paired tibial veins join behind the knee as the popliteal vein. When the calf contracts it squeezes them, and their valves let the blood go only toward the heart: the muscle pump that returns the blood against gravity.',
        systems: ONLY({ veins: 'solid', skeleton: 'glass' }), select: 'FMA44328', view: 'back-right', frame: 'select',
      },
      {
        title: 'The femoral vein',
        text: 'The popliteal vein becomes the femoral vein, a deep vein despite its old name "superficial femoral". In the groin it receives the great saphenous vein at the saphenofemoral junction, where failing valves most often begin varicose veins.',
        systems: ONLY({ veins: 'solid', skeleton: 'glass' }), select: 'FMA21188', view: 'front-right', frame: 'select',
      },
      {
        title: 'Into the pelvis',
        text: 'Behind the inguinal ligament the femoral vein becomes the external iliac vein, which joins the internal iliac vein as the common iliac vein. On the left side the common iliac vein passes behind the right common iliac artery.',
        systems: ONLY({ veins: 'solid', arteries: 'solid', skeleton: 'glass' }), select: 'FMA21387', view: 'front', frame: 'select',
      },
      {
        title: 'The inferior vena cava',
        text: 'The two common iliac veins join at the fifth lumbar vertebra as the inferior vena cava, which climbs to the right of the aorta, behind the liver, and through the diaphragm at the eighth thoracic vertebra.',
        systems: ONLY({ veins: 'solid', arteries: 'glass', skeleton: 'glass' }), select: 'FMA10951', view: 'front-right', frame: 'select',
      },
      {
        title: 'The right atrium',
        text: 'The journey ends in the right atrium, which the superior vena cava enters from above. The vessels\' surfaces here are where the atlas\'s next work begins: the veins reconstructed from images, and the flow within them.',
        systems: ONLY({ veins: 'solid', heart: 'solid' }), isolate: ['FMA11359', 'FMA9457', 'FMA10951', 'FMA4720'], select: 'FMA11359', view: 'front-right', frame: 'isolate',
      },
    ],
  },
  {
    key: 'sections',
    title: 'A walk down the trunk in sections',
    about: 'Transverse sections at the levels anatomy names, seen as a CT scan sees them.',
    steps: [
      {
        title: 'The fourth thoracic vertebra',
        text: 'Near the sternal angle: the arch of the aorta, the trachea about to divide, the oesophagus behind it, the superior vena cava on the right.',
        mode: 'sections', plate: { plane: 'transverse', vertebra: 'T4' },
      },
      {
        title: 'The heart\'s chambers',
        text: 'Lower, the section crosses the heart: its chambers, the descending aorta beside the vertebral body, the bronchial trees of the lungs on either side.',
        mode: 'sections', plate: { plane: 'transverse', vertebra: 'T7' },
      },
      {
        title: 'The tenth thoracic vertebra',
        text: 'The oesophagus passes the diaphragm at this level, to the left of the midline; the aorta lies behind it on the left and the inferior vena cava on the right, with the liver\'s segments around it.',
        mode: 'sections', plate: { plane: 'transverse', vertebra: 'T10' },
      },
      {
        title: 'The first lumbar vertebra',
        text: 'The transpyloric plane: the kidneys at their hila on either side of the vertebra, the pancreas across the midline, the great vessels in front of the spine.',
        mode: 'sections', plate: { plane: 'transverse', vertebra: 'L1' },
      },
      {
        title: 'The fourth lumbar vertebra',
        text: 'Where the aorta divides into the common iliac arteries; the inferior vena cava is formed just below, at the fifth.',
        mode: 'sections', plate: { plane: 'transverse', vertebra: 'L4' },
      },
    ],
  },
];
