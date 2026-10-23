// Complexes that need no engine to be drawn: surfaces as words, triangulations as facets. Each
// can be edited in "Your complex" and sent to the engine for its invariants.

const ring = (n, f) => Array.from({ length: n }, (_, i) => f(i)).join('; ');

export const EXAMPLES = [
  { id: 'torus', title: 'Torus', kind: 'words', text: 'abAB', about: 'A square with opposite sides glued alike: one vertex, two loops, one face.' },
  { id: 'genus-2', title: 'Genus 2', kind: 'words', text: 'a1b1A1B1a2b2A2B2', about: 'An octagon glued into the surface with two handles.' },
  { id: 'projective-plane', title: 'Projective plane', kind: 'words', text: 'aa', about: 'A disc with antipodal boundary points glued: its first homology is ℤ/2.' },
  { id: 'klein', title: 'Klein bottle', kind: 'words', text: 'abAb', about: 'A square glued with one pair of sides reversed.' },
  { id: 'mobius', title: 'Möbius band', kind: 'words', text: 'abac', about: 'A square with one pair of sides glued with a twist; b and c form its one boundary circle.' },
  { id: 'dunce-hat', title: 'Dunce hat', kind: 'words', text: 'aaA', about: 'A triangle glued as a a a⁻¹: contractible, with no free face to collapse from.' },
  { id: 'sphere', title: 'Sphere', kind: 'words', text: 'aA', about: 'A disc whose boundary is folded shut along one edge.' },
  { id: 'two-faces', title: 'Two tori', kind: 'words', text: 'abAB cdCD', about: 'Two polygons with nothing in common: two components.' },
  { id: 'tetrahedron', title: 'Tetrahedron', kind: 'facets', text: '0 1 2; 0 1 3; 0 2 3; 1 2 3', about: 'The boundary of the tetrahedron: the smallest triangulated sphere.' },
  { id: 'octahedron', title: 'Octahedron', kind: 'facets', text: '0 2 4; 0 4 3; 0 3 5; 0 5 2; 1 2 4; 1 4 3; 1 3 5; 1 5 2', about: 'The boundary of the octahedron, a sphere of eight triangles.' },
  { id: 'rp2-6', title: 'RP² with 6 vertices', kind: 'facets', text: '0 1 2; 0 2 3; 0 3 4; 0 4 5; 0 5 1; 1 2 4; 2 3 5; 3 4 1; 4 5 2; 5 1 3', about: 'The minimal triangulation of the projective plane (the hemi-icosahedron).' },
  { id: 'torus-7', title: 'Torus with 7 vertices', kind: 'facets', text: `${ring(7, (i) => `${i} ${(i + 1) % 7} ${(i + 3) % 7}`)}; ${ring(7, (i) => `${i} ${(i + 2) % 7} ${(i + 3) % 7}`)}`, about: 'Möbius and Császár\'s minimal triangulation of the torus: every two vertices joined by an edge.' },
  { id: 'simplex-4', title: '∂Δ⁴, a 3-sphere', kind: 'facets', text: '1 2 3 4; 0 2 3 4; 0 1 3 4; 0 1 2 4; 0 1 2 3', about: 'The boundary of the 4-simplex: five tetrahedra forming a 3-sphere.' },
];
