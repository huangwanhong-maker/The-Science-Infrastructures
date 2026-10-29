// The marks of fields drawn with three.js (ADR 0005; the infrastructure's ADR 0053): curves coloured
// along their length, cells of a volume as small shaded boxes, a level surface, and a drape's image as a
// canvas. The arithmetic (colours, images, segments) is renderers/raster.js's; this module only builds
// three.js's objects from it.
//
// three.js is passed in (fieldMarks(THREE)), as the globe's renderer takes it, so that the module loads
// in Node for its checks and another backend can build its own objects from the same arithmetic.

export function fieldMarks(THREE) {
  // A box whose faces are shaded by their direction (the top lightest), so that unlit materials show
  // the cells' form: a vertex's colour multiplies the cell's.
  function shadedBox() {
    const g = new THREE.BoxGeometry(1, 1, 1);
    const n = g.attributes.normal, shade = new Float32Array(n.count * 3);
    for (let v = 0; v < n.count; v++) {
      const x = n.getX(v), y = n.getY(v), z = n.getZ(v);
      const s = y > 0.5 ? 1.0 : y < -0.5 ? 0.62 : Math.abs(x) > 0.5 ? 0.8 : z > 0.5 ? 0.9 : 0.72;
      shade[3 * v] = shade[3 * v + 1] = shade[3 * v + 2] = s;
    }
    g.setAttribute('color', new THREE.BufferAttribute(shade, 3));
    return g;
  }

  return {
    // Line segments (positions: two points a segment) with a colour at each end (r, g, b, a), or one
    // colour for all; in the units and axes of the positions, which the caller places.
    lines(positions, colours = null, { colour = 0xffffff, opacity = 1, depthTest = true } = {}) {
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.BufferAttribute(positions, 3));
      if (colours) g.setAttribute('color', new THREE.BufferAttribute(colours, 4));
      const material = new THREE.LineBasicMaterial({ color: colours ? 0xffffff : colour, vertexColors: !!colours, transparent: true, opacity, depthTest, depthWrite: false });
      return new THREE.LineSegments(g, material);
    },

    // Cells as one instanced mesh: each a box of edge `size` at its centre, coloured ([r, g, b] in sRGB).
    cells(centres, colours, { size = 1, opacity = 1 } = {}) {
      const material = new THREE.MeshBasicMaterial({ vertexColors: true, transparent: opacity < 1, opacity, depthWrite: opacity >= 1 });
      const mesh = new THREE.InstancedMesh(shadedBox(), material, Math.max(1, centres.length));
      mesh.count = centres.length;
      const m = new THREE.Matrix4(), c = new THREE.Color();
      centres.forEach(([x, y, z], q) => {
        m.makeScale(size, size, size).setPosition(x, y, z);
        mesh.setMatrixAt(q, m);
        const rgb = colours[q];
        mesh.setColorAt(q, c.setRGB(rgb[0], rgb[1], rgb[2], THREE.SRGBColorSpace));
      });
      mesh.instanceMatrix.needsUpdate = true;
      if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;
      mesh.frustumCulled = false;
      return mesh;
    },

    // A drape on a sphere about the frame's origin: the image of a sphere sampling (renderers/raster.js,
    // north row first, west column first) on a sphere of the given radius, in the field's frame with z
    // toward its north pole. three.js's sphere has y up and its seam at −x, where u = 0; turned by a
    // quarter about x, its u = 0 lies at longitude −180° and its top at +z, as the image's.
    // It is drawn on its outer side, so that from outside its far half does not show through the near one.
    sphereDrape(canvas, radius, { opacity = 1, segments = 128, inside = false } = {}) {
      const texture = new THREE.CanvasTexture(canvas);
      texture.colorSpace = THREE.SRGBColorSpace;
      const material = new THREE.MeshBasicMaterial({ map: texture, transparent: true, opacity, depthWrite: false, side: inside ? THREE.BackSide : THREE.FrontSide });
      const mesh = new THREE.Mesh(new THREE.SphereGeometry(radius, segments, segments / 2), material);
      mesh.rotation.x = Math.PI / 2;
      const group = new THREE.Group();
      group.add(mesh);
      return group;
    },

    // A level surface ({positions, normals, indices}, visualization/js/isosurface.js) in one colour.
    surface(level, rgb = [0.85, 0.55, 0.65], { opacity = 0.88 } = {}) {
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.BufferAttribute(level.positions, 3));
      g.setAttribute('normal', new THREE.BufferAttribute(level.normals, 3));
      g.setIndex(new THREE.BufferAttribute(level.indices, 1));
      const colour = new THREE.Color().setRGB(rgb[0], rgb[1], rgb[2], THREE.SRGBColorSpace);
      return new THREE.Mesh(g, new THREE.MeshLambertMaterial({ color: colour, side: THREE.DoubleSide, transparent: opacity < 1, opacity }));
    },
  };
}

// A drape's image (renderers/raster.js) on a canvas, in the browser: for a texture, or a map's tiles.
export function canvasOf(image) {
  const canvas = document.createElement('canvas');
  canvas.width = image.width;
  canvas.height = image.height;
  canvas.getContext('2d').putImageData(new ImageData(image.data, image.width, image.height), 0, 0);
  return canvas;
}
