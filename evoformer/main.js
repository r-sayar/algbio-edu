import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

// ---------------------------------------------------------------------------
// Three.js scene boot
// ---------------------------------------------------------------------------
const canvas = document.getElementById('c');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(window.devicePixelRatio);
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setClearColor(0x0e1116);

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(
  35, window.innerWidth / window.innerHeight, 0.1, 2000
);

const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true;
controls.dampingFactor = 0.08;

scene.add(new THREE.AmbientLight(0xffffff, 0.55));
const dl = new THREE.DirectionalLight(0xffffff, 0.85);
dl.position.set(8, 14, 6);
scene.add(dl);
const dl2 = new THREE.DirectionalLight(0x88aaff, 0.25);
dl2.position.set(-10, -4, -8);
scene.add(dl2);

window.addEventListener('resize', () => {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
});

// Faint ground grid for reference
const grid = new THREE.GridHelper(80, 40, 0x303644, 0x1c2029);
grid.position.y = -40;
scene.add(grid);


// ---------------------------------------------------------------------------
// Color: blue (negative) -> grey (zero) -> red (positive), clamp at ±2.
// ---------------------------------------------------------------------------
function valueToColor(v, target) {
  const t = Math.max(-1, Math.min(1, v / 2.0));
  if (t < 0) {
    const a = -t;
    target.setRGB(0.50 - 0.30 * a, 0.50 - 0.30 * a, 0.50 + 0.50 * a);
  } else {
    target.setRGB(0.50 + 0.50 * t, 0.50 - 0.20 * t, 0.50 - 0.20 * t);
  }
  return target;
}

// ---------------------------------------------------------------------------
// Build a glassy voxel grid for a 3-D tensor.
// Returns a Group with two children:
//   - .userData.body  : InstancedMesh of transparent cubes (per-cell colored)
//   - .userData.edges : LineSegments of cube wireframes (so the grid is legible)
// ---------------------------------------------------------------------------
const CUBE_SIZE = 0.72;
const CUBE_GAP  = 0.34;
const STEP      = CUBE_SIZE + CUBE_GAP;

function buildTensorMesh(dimSize, data, originXYZ) {
  const [s0, s1, s2] = dimSize;
  const total = s0 * s1 * s2;

  // ---- body: transparent cubes, per-instance color ----
  const bodyGeo = new THREE.BoxGeometry(CUBE_SIZE, CUBE_SIZE, CUBE_SIZE);
  const bodyMat = new THREE.MeshStandardMaterial({
    metalness: 0.0, roughness: 0.30,
    transparent: true, opacity: 0.45,
    depthWrite: false,
  });
  const inst = new THREE.InstancedMesh(bodyGeo, bodyMat, total);
  inst.renderOrder = 1;

  // ---- edges: bake every cube's 12 edges into one LineSegments ----
  const edgeTemplate = new THREE.EdgesGeometry(bodyGeo);
  const edgeVertCount = edgeTemplate.attributes.position.count;  // 24 per cube
  const positions = new Float32Array(total * edgeVertCount * 3);
  const src = edgeTemplate.attributes.position.array;

  const dummy = new THREE.Object3D();
  const color = new THREE.Color();
  let k = 0;
  let p = 0;
  for (let a = 0; a < s0; a++) {
    for (let b = 0; b < s1; b++) {
      for (let cI = 0; cI < s2; cI++) {
        const v = data[k];
        const px = originXYZ[0] + b * STEP;
        const py = originXYZ[1] + cI * STEP;
        const pz = originXYZ[2] + a * STEP;

        dummy.position.set(px, py, pz);
        dummy.updateMatrix();
        inst.setMatrixAt(k, dummy.matrix);
        inst.setColorAt(k, valueToColor(v, color));

        // translate the edge template into this cell's position
        for (let q = 0; q < edgeVertCount; q++) {
          positions[p++] = src[q * 3 + 0] + px;
          positions[p++] = src[q * 3 + 1] + py;
          positions[p++] = src[q * 3 + 2] + pz;
        }
        k++;
      }
    }
  }
  inst.instanceMatrix.needsUpdate = true;
  inst.instanceColor.needsUpdate = true;
  edgeTemplate.dispose();

  const lineGeo = new THREE.BufferGeometry();
  lineGeo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  const lineMat = new THREE.LineBasicMaterial({
    color: 0xeaf0ff, transparent: true, opacity: 0.35,
    depthWrite: false,
  });
  const lines = new THREE.LineSegments(lineGeo, lineMat);
  lines.renderOrder = 2;

  const group = new THREE.Group();
  group.add(inst);
  group.add(lines);
  group.userData = { body: inst, edges: lines };
  return group;
}

// Sprite-based text label.
function makeLabel(text, color = '#ffffff', sizePx = 30, scale = 8) {
  const c = document.createElement('canvas');
  c.width = 1024; c.height = 128;
  const ctx = c.getContext('2d');
  ctx.font = `600 ${sizePx}px -apple-system, system-ui, sans-serif`;
  ctx.fillStyle = color;
  ctx.textBaseline = 'middle';
  ctx.fillText(text, 8, 64);
  const tex = new THREE.CanvasTexture(c);
  tex.minFilter = THREE.LinearFilter;
  tex.anisotropy = 4;
  const mat = new THREE.SpriteMaterial({ map: tex, transparent: true });
  const sp = new THREE.Sprite(mat);
  sp.scale.set(scale, scale * (c.height / c.width), 1);
  return sp;
}

// Arrow from a -> b (vec3), colored.
function makeArrow(from, to, hex = 0x9aa6b6, headLen = 1.0, headWid = 0.6) {
  const dir = new THREE.Vector3().subVectors(to, from);
  const len = dir.length();
  const a = new THREE.ArrowHelper(dir.clone().normalize(), from, len,
                                  hex, headLen, headWid);
  a.line.material.linewidth = 2;
  return a;
}


// ---------------------------------------------------------------------------
// Load tensors.json and build the scene
// ---------------------------------------------------------------------------
async function load() {
  const resp = await fetch('./tensors.json');
  const payload = await resp.json();
  const { N_seq, N_res, c, snapshots } = payload;

  const sceneRoot = new THREE.Group();
  scene.add(sceneRoot);

  // --- spatial layout ---
  // Each snapshot is a "level".  Levels stacked along Y (top -> input,
  // bottom -> output), so visual flow matches the architecture diagram.
  // Within a level, m sits at -X, z sits at +X.
  const m_w = N_res * STEP;             // width of m grid in scene units
  const m_d = N_seq * STEP;
  const z_w = N_res * STEP;
  const z_d = N_res * STEP;
  const cubesH = c * STEP;

  const levelH = cubesH + 8;            // vertical gap between snapshots
  const m_x = -10;
  const z_x = +5;

  const levelInfo = [];                 // for ui & arrows
  const levelGroups = [];

  snapshots.forEach((snap, i) => {
    const yTop = (snapshots.length - 1) * levelH * 0.5 - i * levelH;

    const m_origin = [m_x, yTop, -m_d * 0.5];
    const z_origin = [z_x, yTop, -z_d * 0.5];

    const lvl = new THREE.Group();
    const mGroup = buildTensorMesh(snap.m.shape, snap.m.data, m_origin);
    const zGroup = buildTensorMesh(snap.z.shape, snap.z.data, z_origin);
    lvl.add(mGroup, zGroup);
    sceneRoot.add(lvl);
    levelGroups.push({ lvl, mGroup, zGroup });

    // Per-level labels
    const sublayerLabel = makeLabel(snap.name, '#ffd479', 30, 11);
    sublayerLabel.position.set((m_x + z_x) / 2 + (m_w + z_w) / 4 - 2,
                                yTop + cubesH + 1.3, 0);
    lvl.add(sublayerLabel);

    if (i === 0) {
      const lm = makeLabel(`m  (${snap.m.shape.join(' × ')})`, '#9ec5ff', 26, 5);
      lm.position.set(m_origin[0] + m_w / 2 - 1, yTop + cubesH + 0.4, 0);
      lvl.add(lm);
      const lz = makeLabel(`z  (${snap.z.shape.join(' × ')})`, '#ff9da4', 26, 5);
      lz.position.set(z_origin[0] + z_w / 2 - 1, yTop + cubesH + 0.4, 0);
      lvl.add(lz);
    }

    levelInfo.push({
      yCenter: yTop + cubesH * 0.5,
      yTop: yTop + cubesH,
      yBot: yTop,
      m_origin, z_origin, mGroup, zGroup, name: snap.name,
    });
  });

  // --- vertical flow arrows between consecutive levels ---
  for (let i = 0; i < levelInfo.length - 1; i++) {
    const a = levelInfo[i], b = levelInfo[i + 1];
    // arrow on m stream
    sceneRoot.add(makeArrow(
      new THREE.Vector3(m_x + m_w / 2 - 1, a.yBot - 0.6, 0),
      new THREE.Vector3(m_x + m_w / 2 - 1, b.yTop + 0.6, 0),
      0x6c89ff,
    ));
    // arrow on z stream
    sceneRoot.add(makeArrow(
      new THREE.Vector3(z_x + z_w / 2 - 1, a.yBot - 0.6, 0),
      new THREE.Vector3(z_x + z_w / 2 - 1, b.yTop + 0.6, 0),
      0xff7a85,
    ));
  }

  // --- cross-stream arrows that explain the data flow ---
  // pair -> MSA bias (snap 0 -> sublayer 1 = level index 1)
  if (levelInfo.length >= 2) {
    const lvl1 = levelInfo[1];
    sceneRoot.add(makeArrow(
      new THREE.Vector3(z_x - 0.5, lvl1.yCenter, z_d * 0.6),
      new THREE.Vector3(m_x + m_w + 0.5, lvl1.yCenter, m_d * 0.6),
      0x6dd56d, 1.4, 0.9,
    ));
    const tag = makeLabel('pair → MSA  (bias)', '#6dd56d', 24, 6);
    tag.position.set((m_x + z_x) / 2, lvl1.yCenter + 0.3, z_d * 0.7 + 1);
    sceneRoot.add(tag);
  }
  // MSA -> pair OPM (snap 1 -> snap 2)
  if (levelInfo.length >= 3) {
    const lvl2 = levelInfo[2];
    sceneRoot.add(makeArrow(
      new THREE.Vector3(m_x + m_w + 0.5, lvl2.yCenter, -z_d * 0.6),
      new THREE.Vector3(z_x - 0.5, lvl2.yCenter, -z_d * 0.6),
      0x6dd56d, 1.4, 0.9,
    ));
    const tag = makeLabel('MSA → pair  (OPM)', '#6dd56d', 24, 6);
    tag.position.set((m_x + z_x) / 2, lvl2.yCenter + 0.3, -z_d * 0.7 - 1);
    sceneRoot.add(tag);
  }

  // --- title above the whole tower ---
  const title = makeLabel('Evoformer block — full data flow', '#fff', 36, 16);
  title.position.set((m_x + z_x) / 2 + (m_w + z_w) / 4 - 2,
                     levelInfo[0].yTop + 5, 0);
  scene.add(title);

  // --- camera framing: fit all levels in view ---
  const yHigh = levelInfo[0].yTop + 4;
  const yLow  = levelInfo[levelInfo.length - 1].yBot - 4;
  const cy = (yHigh + yLow) / 2;
  const span = yHigh - yLow;
  camera.position.set(28, cy + span * 0.05, 28);
  controls.target.set(0, cy, 0);
  controls.update();

  // --- step UI: highlight by dimming non-current levels ---
  let idx = 0;
  const stepName = document.getElementById('stepname');
  const stepIdx  = document.getElementById('stepidx');

  function setStep(i) {
    idx = ((i % snapshots.length) + snapshots.length) % snapshots.length;
    levelGroups.forEach(({ mGroup, zGroup }, j) => {
      const focused = j === idx;
      const bodyOp = focused ? 0.55 : 0.13;
      const edgeOp = focused ? 0.55 : 0.10;
      mGroup.userData.body.material.opacity  = bodyOp;
      zGroup.userData.body.material.opacity  = bodyOp;
      mGroup.userData.edges.material.opacity = edgeOp;
      zGroup.userData.edges.material.opacity = edgeOp;
    });
    stepName.textContent = snapshots[idx].name;
    stepIdx.textContent = `(${idx + 1} / ${snapshots.length})`;

    // Pan the orbit target to the highlighted level (smooth)
    targetY = levelInfo[idx].yCenter;
  }

  let targetY = levelInfo[0].yCenter;
  setStep(0);

  document.getElementById('prev').onclick = () => setStep(idx - 1);
  document.getElementById('next').onclick = () => setStep(idx + 1);
  document.addEventListener('keydown', (e) => {
    if (e.key === 'ArrowLeft' || e.key === 'ArrowUp')   setStep(idx - 1);
    if (e.key === 'ArrowRight' || e.key === 'ArrowDown') setStep(idx + 1);
  });

  let playing = false, playT = 0;
  const playBtn = document.getElementById('play');
  playBtn.onclick = () => {
    playing = !playing;
    playBtn.textContent = playing ? 'pause' : 'play';
    playT = 0;
  };
  let prev = performance.now();

  function tick(now) {
    const dt = (now - prev) / 1000; prev = now;
    if (playing) {
      playT += dt;
      if (playT > 1.6) { playT = 0; setStep(idx + 1); }
    }
    // smooth-pan target.y toward focused level
    controls.target.y += (targetY - controls.target.y) * Math.min(1, dt * 2.5);
    controls.update();
    renderer.render(scene, camera);
    requestAnimationFrame(tick);
  }
  requestAnimationFrame(tick);
}

load().catch((e) => {
  console.error(e);
  document.getElementById('stepname').textContent = 'load error: ' + e.message;
});
