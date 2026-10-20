#!/usr/bin/env node
// make_unified_slides.js -- editable PowerPoint version of the unified-annotation
// scheme (how Repeat_Annotation_Unified.gff3 is assembled), for presentations.
//
//   npm install pptxgenjs            # once, anywhere on NODE_PATH
//   node scripts/make_unified_slides.js [figs/unified_annotation_scheme.pptx]
//
// Five slides, native shapes only (every glyph is a named shape in the Selection
// Pane, so the deck can be adjusted by hand):
//   1  tier priority -- six evidence tiers, each base goes to the highest one;
//   2  tandem LTR-RT (LTR_RT_TR container, member copies nested);
//   3  tandem repeat derived from an LTR-RT (TE_origin, satellite wins);
//   4  LTR-RT inserted into a satellite (element wins, satellite split);
//   5  satellite-library hit vs TE hit (rm_tc_tandem_gate).
// The GFF3 attributes of each case are in the speaker notes.
//
// Same visual vocabulary and palette as scripts/make_workflow_figure.R
// (figs/carp_workflow.png). Like that figure it is hand-authored and names tools
// and attributes, never thresholds, so it cannot drift from
// make_unified_annotation.R; the result rows are computed by interval
// subtraction from the drawn input rather than drawn to look right.
//
// Optional: set PPTX_SKILL to a pptx skill directory providing
// scripts/apply_theme.js to also write the palette into the deck's theme.
// Without it the deck is identical except for the theme's accent colours.
const pptxgen = require("pptxgenjs");

const OUT = process.argv[2] || "figs/unified_annotation_scheme.pptx";

const THEME = {
  name: "CARP",
  headFontFace: "Calibri",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "14304A", lt1: "FFFFFF", dk2: "5F6E7A", lt2: "F4F6F8",
    accent1: "0072B2", accent2: "D55E00", accent3: "009E73",
    accent4: "8E7F6E", accent5: "CC79A7", accent6: "14304A",
    hlink: "0072B2", folHlink: "5F6E7A",
  },
};
// Glyph colours are hex (fills need transparency + derived dark shades);
// they equal the theme accents so a restyle stays consistent.
const K = {
  ltr: "0072B2", ltrDk: "003F63", tir: "D55E00", tirDk: "763400",
  tan: "009E73", tanDk: "00573F", ssr: "8E7F6E",
  track: "D8D8D8", rule: "DDE2E7", ink: "14304A", mute: "5F6E7A", card: "F4F6F8",
};

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.333 x 7.5
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "How the unified repeat annotation is built";
pres.author = "CARP";
const S = pres.ShapeType;

pres.defineSlideMaster({
  title: "CASE",
  background: { color: "FFFFFF" },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.5, y: 0.32, w: 12.3, h: 0.7,
        fontSize: 32, bold: true, color: THEME.colors.dk1, valign: "middle", align: "left", margin: 0 }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.5, y: 1.0, w: 12.3, h: 0.45,
        fontSize: 16, color: THEME.colors.dk2, valign: "top", margin: 0 }, text: "" } },
  ],
  slideNumber: { x: 12.4, y: 7.05, w: 0.5, h: 0.3, fontSize: 10, color: THEME.colors.dk2, align: "right" },
});
pres.addSection({ title: "Unified annotation" });

// ---------------------------------------------------------------- primitives
let N = 0;
const nm = (s) => `${s} ${++N}`;
// genome-coordinate axis: u in 0..100 -> inches
const axis = (x0, x1) => (u) => x0 + (u / 100) * (x1 - x0);
const W = (X, a, b) => X(b) - X(a);

function track(sl, X, y) {
  sl.addShape(S.roundRect, { objectName: nm("genome"), x: X(0), y: y - 0.035, w: W(X, 0, 100), h: 0.07,
    fill: { color: K.track }, line: { type: "none" }, rectRadius: 0.035 });
}
function box(sl, x, y, w, h, fill, line, o = {}) {
  sl.addShape(S.roundRect, Object.assign({ x, y: y - h / 2, w, h, rectRadius: Math.min(0.05, h / 2),
    fill: { color: fill, transparency: o.t || 0 },
    line: line ? { color: line, width: o.lw || 1 } : { type: "none" } }, o.name ? { objectName: o.name } : {}));
}
const FADE = 80;
// a similarity hit (flat bar)
function hit(sl, X, a, b, y, col, faded, h = 0.2) {
  if (faded) box(sl, X(a), y, W(X, a, b), h, col, null, { t: FADE, name: nm("hit, claimed above") });
  else box(sl, X(a), y, W(X, a, b), h, col, col, { t: 20, name: nm("hit") });
}
// complete LTR-RT; ltrs = genome positions where LTR blocks start
function ltrElement(sl, X, a, b, y, opts = {}) {
  const h = opts.h || 0.3, lw = opts.lw || 2.2, faded = opts.faded;
  if (faded) { box(sl, X(a), y, W(X, a, b), h, K.ltr, null, { t: FADE, name: nm("LTR-RT, claimed above") }); return; }
  box(sl, X(a), y, W(X, a, b), h, K.ltr, K.ltr, { t: 18, lw: 1.25, name: nm("LTR-RT body") });
  for (const l of opts.ltrs || [a, b - lw])
    box(sl, X(l), y, W(X, l, l + lw), h, K.ltrDk, null, { name: nm("LTR") });
}
function tirElement(sl, X, a, b, y, h = 0.3) {
  box(sl, X(a), y, W(X, a, b), h, K.tir, K.tir, { t: 18, lw: 1.25, name: nm("TIR element") });
  // Terminal inverted repeats as outward-pointing triangles. PowerPoint rotates a
  // shape about its centre, so place the UNROTATED box (tall x deep, swapped) such
  // that after the 90-degree turn it covers [edge, edge + deep] and is centred on y.
  const tall = h * 0.84, deep = 0.14, pad = 0.03;
  const tri = (cx, rot) => sl.addShape(S.triangle, { objectName: nm("TIR"),
    x: cx - tall / 2, y: y - deep / 2, w: tall, h: deep, rotate: rot,
    fill: { color: K.tirDk }, line: { type: "none" } });
  tri(X(a) + pad + deep / 2, 270);   // apex points left, out of the element
  tri(X(b) - pad - deep / 2, 90);    // apex points right
}
// protein domain (pentagon arrow)
function domain(sl, X, a, b, y, col, o = {}) {
  const h = o.h || 0.24;
  sl.addShape(S.homePlate, { objectName: nm("domain"), x: X(a), y: y - h / 2, w: W(X, a, b), h,
    fill: { color: o.fill || col, transparency: o.faded ? FADE : (o.ft ?? 40) },
    line: o.faded ? { type: "none" } : { color: o.line || col, width: 1 } });
  if (o.lab) label(sl, X(a) - 0.2, y + h / 2 + 0.02, W(X, a, b) + 0.4, 0.22, o.lab,
                   { fontSize: 10, color: K.mute, align: "center" });
}
// tandem array: one box per monomer. keep = [[a,b],...] solid intervals (null = all)
function satellite(sl, X, a, b, y, o = {}) {
  const h = o.h || 0.3, n = o.n || Math.max(3, Math.round((b - a) / 1.25)), w = (b - a) / n;
  const gap = Math.min(0.025, W(X, 0, w) * 0.12);
  for (let i = 0; i < n; i++) {
    const m0 = a + i * w, mc = m0 + w / 2;
    const solid = !o.keep || o.keep.some(([p, q]) => mc >= p && mc <= q);
    const x = X(m0) + gap, ww = W(X, 0, w) - 2 * gap;
    if (!solid) box(sl, x, y, ww, h, K.tan, null, { t: FADE, name: nm("monomer, claimed above") });
    else if (o.open) box(sl, x, y, ww, h, "FFFFFF", K.tan, { name: nm("monomer (Satellite/Unknown)") });
    else box(sl, x, y, ww, h, K.tan, K.tan, { t: 15, name: nm("monomer") });
    if (solid && o.ltr) box(sl, x, y, W(X, 0, o.ltr), h, K.tanDk, null, { name: nm("monomer LTR part") });
  }
}
// Level-2 child: thin bar tied to its parent by a dashed line
function child(sl, X, a, b, y, col, yParentBottom, h = 0.13) {
  const xm = X((a + b) / 2);
  sl.addShape(S.line, { objectName: nm("parent link"), x: xm, y: yParentBottom, w: 0, h: y - h / 2 - yParentBottom,
    line: { color: col, width: 1, dashType: "dash" } });
  box(sl, X(a), y, W(X, a, b), h, col, col, { t: 20, name: nm("Level-2 child") });
}
function label(sl, x, y, w, h, text, o = {}) {
  sl.addText(text, Object.assign({ x, y, w, h, isTextBox: true, margin: 0, fontSize: 14,
    color: K.ink, valign: "middle" }, o));
}
function rowLabel(sl, xr, y, main, sub) {
  if (!sub) { label(sl, xr - 2.6, y - 0.16, 2.6, 0.32, main, { align: "right", bold: true }); return; }
  sl.addText([
    { text: main, options: { bold: true, fontSize: 14, color: K.ink, breakLine: true } },
    { text: sub, options: { fontSize: 11, color: K.mute } },
  ], { x: xr - 3.0, y: y - 0.3, w: 3.0, h: 0.6, isTextBox: true, margin: 0, align: "right", valign: "middle" });
}
function hline(sl, x0, x1, y) {
  sl.addShape(S.line, { objectName: nm("divider"), x: x0, y, w: x1 - x0, h: 0, line: { color: K.rule, width: 1 } });
}
function downArrow(sl, x, y0, y1) {
  sl.addShape(S.line, { objectName: nm("arrow"), x, y: y0, w: 0, h: y1 - y0,
    line: { color: K.ink, width: 2, endArrowType: "triangle" } });
}
function legendKey(sl, x, y, col, text, style) {
  if (style === "open") box(sl, x, y, 0.28, 0.2, "FFFFFF", col, {});
  else if (style === "faded") box(sl, x, y, 0.28, 0.2, col, null, { t: FADE });
  else box(sl, x, y, 0.28, 0.2, col, col, { t: 18 });
  label(sl, x + 0.36, y - 0.15, 2.4, 0.3, text, { fontSize: 12, color: K.mute });
}
// right-hand "In the GFF3" card
// The GFF3 attributes of each case go to the speaker notes, not onto the slide.
function gffCard(sl, x, y, w, h, lines) {
  sl._gff = "In the unified GFF3:\n" + lines.map((l) =>
    l.head ? `\n${l.head}` : l.code ? `  ${l.code.trim()}` : `  ${l.text}`).join("\n");
}
// "annotated as" label: what the unified GFF3 calls this feature
function tag(sl, x, y, w, what, extra, align = "center") {
  const runs = [{ text: "annotated as  ", options: { color: K.mute } }, { text: what, options: { bold: true, color: K.ink } }];
  if (extra) runs.push({ text: "  " + extra, options: { color: K.mute } });
  sl.addText(runs, { x, y, w, h: 0.3, isTextBox: true, margin: 0, fontSize: 13, align, valign: "middle",
    objectName: nm("annotation label") });
}
function caseSlide(title, sub) {
  const sl = pres.addSlide({ masterName: "CASE", sectionTitle: "Unified annotation" });
  sl.addText(title, { placeholder: "title" });
  sl.addText(sub, { placeholder: "body" });
  const addNotes = sl.addNotes.bind(sl);
  sl.addNotes = (t) => addNotes(sl._gff ? `${t}\n\n${sl._gff}` : t);
  return sl;
}
// interval subtraction, so the result rows are what the rule gives
function subtract(iv, blockers) {
  let out = [iv];
  for (const [p, q] of blockers) {
    const nxt = [];
    for (const [a, b] of out) {
      if (q <= a || p >= b) { nxt.push([a, b]); continue; }
      if (p > a) nxt.push([a, p]);
      if (q < b) nxt.push([q, b]);
    }
    out = nxt;
  }
  return out.filter(([a, b]) => b - a > 0.2);
}

// ============================================================ slide 1
{
  const sl = caseSlide("Each base goes to the highest tier that claims it",
    "Six evidence tiers, best first. A lower tier keeps only the sequence no higher tier has claimed.");
  const X = axis(3.75, 12.83);
  const R = {
    ltr: [3, 35], ssr_in: [30.8, 32.4], degr: [38.5, 52], orphan: [43, 47.5],
    ssr_free: [53.6, 55.2], tir: [57.2, 65.4], tpase: [59.6, 63],
    sat_tc: [70.5, 85], sat_rm: [67.8, 87.8], th_over: [71.5, 84], th_free: [90.6, 95.8],
  };
  const DOM = [["GAG", 6.2, 10.2], ["PROT", 11, 14.4], ["RT", 15.2, 20.4], ["RH", 21.2, 24.8], ["INT", 25.6, 30]];
  const TY = [1.8, 2.34, 2.88, 3.42, 3.96, 4.50];
  const labs = [["1  Complete elements", "DANTE_LTR · DANTE_TIR · DANTE_LINE"],
                ["2  Protein domains", "DANTE"],
                ["3  Tandem arrays", "TideCluster"],
                ["4  Satellite-library hits", "RepeatMasker, TideCluster library"],
                ["5  Repeat-library hits", "RepeatMasker, TE library + simple repeats"],
                ["6  Tandem leftovers", "TideHunter"]];
  TY.forEach((y, i) => { track(sl, X, y); rowLabel(sl, 3.55, y, labs[i][0], labs[i][1]); });

  const t1 = [R.ltr, R.tir], t3 = [R.sat_tc];
  const t2 = subtract(R.orphan, t1);
  const t4 = subtract(R.sat_rm, [...t1, ...t2, ...t3]);
  const claimed = [...t1, ...t2, ...t3, ...t4];

  ltrElement(sl, X, ...R.ltr, TY[0]);
  tirElement(sl, X, ...R.tir, TY[0]);
  DOM.forEach(([n, a, b]) => domain(sl, X, a, b, TY[1], K.ltr, { faded: true }));
  domain(sl, X, ...R.tpase, TY[1], K.tir, { faded: true });
  domain(sl, X, ...R.orphan, TY[1], K.ltr, { lab: "RT" });
  satellite(sl, X, ...R.sat_tc, TY[2]);
  satellite(sl, X, ...R.sat_rm, TY[3], { h: 0.24, keep: t4 });
  for (const iv of [R.ltr, R.degr, R.tir]) {
    const col = iv === R.tir ? K.tir : K.ltr;
    hit(sl, X, ...iv, TY[4], col, true);
    for (const p of subtract(iv, claimed)) hit(sl, X, ...p, TY[4], col, false);
  }
  hit(sl, X, ...R.ssr_in, TY[4], K.ssr, false);
  hit(sl, X, ...R.ssr_free, TY[4], K.ssr, false);
  satellite(sl, X, ...R.th_over, TY[5], { h: 0.22, n: 9, keep: [] });
  satellite(sl, X, ...R.th_free, TY[5], { h: 0.22, n: 4, open: true });

  // result
  downArrow(sl, X(56), 4.78, 5.18);
  hline(sl, 0.5, X(100), 5.3);
  const L1 = 5.92, L2 = 6.5;
  label(sl, 0.5, 5.36, 3.05, 0.34, "Unified annotation", { bold: true, fontSize: 18, align: "right" });
  track(sl, X, L1);
  rowLabel(sl, 3.55, L1, "Level 1", "non-overlapping");
  rowLabel(sl, 3.55, L2, "Level 2", "nested children");
  ltrElement(sl, X, ...R.ltr, L1);
  tirElement(sl, X, ...R.tir, L1);
  for (const p of subtract(R.degr, claimed)) hit(sl, X, ...p, L1, K.ltr, false);
  domain(sl, X, ...R.orphan, L1, K.ltr);
  hit(sl, X, ...R.ssr_free, L1, K.ssr, false);
  satellite(sl, X, ...R.sat_rm, L1);
  satellite(sl, X, ...R.th_free, L1, { h: 0.22, n: 4, open: true });
  child(sl, X, ...R.ssr_in, L2, K.ssr, L1 + 0.15);
  label(sl, X(R.ssr_in[1]) + 0.12, L2 - 0.15, 4.5, 0.3, "a simple repeat inside an element is nested, not cut out",
        { fontSize: 12, color: K.mute });

  let lx = 0.5; const ly = 7.08;
  for (const [c, t, st, w] of [[K.ltr, "LTR retrotransposon", null, 2.1], [K.tir, "TIR DNA transposon", null, 2.05],
       [K.tan, "tandem repeat", null, 1.6], [K.tan, "Satellite/Unknown", "open", 1.9], [K.ssr, "simple repeat", null, 1.6],
       [K.ink, "claimed by a higher tier", "faded", 2.6]]) { legendKey(sl, lx, ly, c, t, st); lx += w; }

  sl.addNotes(
    "Each row is one evidence layer, ordered by how much we trust it. Structure-based calls (complete elements, " +
    "protein domains) come first, then tandem arrays, then similarity searches, then TideHunter leftovers.\n" +
    "Faded parts are what a higher tier already claimed; only the solid parts reach the unified annotation.\n" +
    "Example: RepeatMasker hits the complete LTR-RT too (its consensus came from it), but that span is already " +
    "claimed. The degraded copy is seen only by RepeatMasker plus one surviving RT domain: the domain keeps its span, " +
    "RepeatMasker the flanks. Simple repeats are the exception: inside an element they become nested Level-2 children.");
}

// ============================================================ slide 2
{
  const sl = caseSlide("Tandem LTR-RT: one container, copies nested",
    "Same-lineage copies head-to-tail share their middle LTR. DANTE_LTR reports every copy, so that LTR would be counted twice.");
  const X = axis(2.35, 12.83);
  const lw = 3.6, e1 = [6, 52], e2 = [52 - lw, 98 - lw];
  let y = 2.35;
  track(sl, X, y + 0.25);
  rowLabel(sl, 2.15, y + 0.25, "DANTE_LTR", "two complete copies");
  ltrElement(sl, X, ...e1, y, { lw });
  ltrElement(sl, X, ...e2, y + 0.5, { lw });
  label(sl, X(e2[0]) - 0.9, y + 0.86, 1.8 + W(X, 0, lw), 0.3, "shared LTR", { fontSize: 12, color: K.ltr, align: "center" });
  downArrow(sl, X(50), 3.55, 4.0);
  hline(sl, 0.5, 12.83, 4.12);
  const L1 = 4.75, L2 = 5.55;
  track(sl, X, L1);
  rowLabel(sl, 2.15, L1, "Level 1", "container");
  rowLabel(sl, 2.15, L2 + 0.1, "Level 2", "member copies");
  ltrElement(sl, X, e1[0], e2[1], L1, { lw, ltrs: [e1[0], e2[0], e2[1] - lw] });
  tag(sl, X(e1[0]), L1 - 0.5, W(X, e1[0], e2[1]), "LTR_RT_TR", "(tandem LTR-RT, copy_number = 2)");
  const ch = 0.13;
  [[e1, L2], [e2, L2 + 0.22]].forEach(([e, yy]) => {
    child(sl, X, ...e, yy, K.ltr, L1 + 0.15, ch);
    for (const l of [e[0], e[1] - lw]) box(sl, X(l), yy, W(X, l, l + lw), ch, K.ltrDk, null, { name: nm("member LTR") });
  });
  tag(sl, X(e1[0]), L2 + 0.42, W(X, e1[0], e2[1]), "LTR-RT", "(member copy, nested under the container)");
  label(sl, 0.5, 6.5, 12.3, 0.6,
    "Detected when adjacent copies have the same lineage and strand and their facing LTRs coincide. A nested insertion does not qualify.",
    { fontSize: 14, color: K.mute, valign: "top" });
  gffCard(sl, 9.2, 1.75, 3.63, 3.75, [
    { head: "Level 1 container" }, { code: "source  DANTE_LTR" }, { code: "structure=LTR_RT_TR" }, { code: "copy_number=2" },
    { head: "Level 2, one per copy" }, { code: "in_structure=LTR_RT_TR" }, { code: "member_of=<container ID>" }, { code: "Parent=<container ID>" },
    { head: "Density tracks" }, { text: "LTR_RT_TR also gets its own roll-up BigWig." },
  ]);
  sl.addNotes(
    "DANTE_LTR emits each complete copy. When same-lineage copies sit head-to-tail sharing an LTR " +
    "(LTR-INT-LTR-INT-LTR), the copies overlap on that LTR. resolve_ltr_tandems.py collapses the chain into one " +
    "container spanning the array; the copies are kept as nested Level-2 children, so nothing is lost and nothing is " +
    "double-counted. A nested insertion overlaps deeply but its LTRs do not line up, so it is not merged.");
}

// ============================================================ slide 3
{
  const sl = caseSlide("Tandem repeat derived from an LTR-RT: the satellite wins",
    "The monomer is a whole LTR-RT, so its protein domains recur once per period. The TE family is kept as an attribute.");
  const X = axis(2.35, 12.83);
  const nmono = 5, x0 = 3, x1 = 97, P = (x1 - x0) / nmono, lw = 3.2;
  let y = 1.95;
  track(sl, X, y); rowLabel(sl, 2.15, y, "TideCluster", "tandem array");
  satellite(sl, X, x0, x1, y, { n: nmono, ltr: lw });
  label(sl, X(60), y - 0.42, W(X, 60, 100), 0.24, "monomer = LTR + internal region",
        { fontSize: 12, color: K.tan, align: "right" });
  y = 2.65; track(sl, X, y); rowLabel(sl, 2.15, y, "DANTE", "domains");
  const present = [[0, 1, 2], [0, 1, 2], [0, 2], [0, 1, 2], [1]];
  present.forEach((js, i) => js.forEach((j) => {
    const a = x0 + i * P + lw + 1.2 + j * 4.6; domain(sl, X, a, a + 3.9, y, K.ltr, { h: 0.2 });
  }));
  y = 3.35; track(sl, X, y); rowLabel(sl, 2.15, y, "DANTE_LTR", "complete copies");
  ltrElement(sl, X, x0, x0 + P + lw, y - 0.17, { lw, h: 0.26 });
  ltrElement(sl, X, x0 + P, x0 + 2 * P + lw, y + 0.17, { lw, h: 0.26 });
  ltrElement(sl, X, x0 + 3 * P, x0 + 4 * P + lw, y, { lw, h: 0.26 });
  label(sl, X(x0), y + 0.36, 4.6, 0.28, "2 copies sharing an LTR = LTR_RT_TR container",
        { fontSize: 12, color: K.ltr, align: "left" });
  downArrow(sl, X(50), 4.1, 4.5);
  hline(sl, 0.5, 12.83, 4.62);
  const L1 = 5.25;
  track(sl, X, L1); rowLabel(sl, 2.15, L1, "Level 1", "unified");
  satellite(sl, X, x0, x1, L1, { n: nmono, ltr: lw });
  tag(sl, X(x0), L1 - 0.52, W(X, x0, x1), "TRC_<n>  tandem repeat", "of LTR-RT origin (TE_origin = Class_I/LTR/Ty3_gypsy/…)");
  label(sl, 0.5, 5.75, 12.3, 1.2,
    [{ text: "Not shown in the unified file: ", options: { bold: true, color: K.ink } },
     { text: "the LTR-RT calls and the container. They stay in DANTE_LTR.gff3. Degraded copies have lost some domains, but domains still occur in most periods.", options: { color: K.mute } }],
    { fontSize: 14, valign: "top" });
  gffCard(sl, 9.2, 1.75, 3.63, 4.6, [
    { head: "Level 1 satellite" }, { code: "source  TideCluster_default" }, { code: "Name=TRC_<n>" },
    { code: "TE_origin=" }, { code: "  Class_I/LTR/Ty3_gypsy/…" }, { code: "TE_origin_structure=" }, { code: "  tandem_LTR_RT" },
    { head: "Required" }, { text: "≥2 complete same-family TEs under the array, and their domains recurring with the tandem period." },
    { head: "Also written to" }, { code: "…te_derived_trc.csv" },
  ]);
  sl.addNotes(
    "This is the only case where a satellite outranks a tier-1 structural call. Two tests must pass: the array covers " +
    "several complete TEs of one family, and their DANTE domains recur with the array's tandem period (the domain-rhythm " +
    "gate). Domains are used because they survive decay of the element, so a degraded tandem still passes.\n" +
    "TE_origin records the shared lineage of the covered TEs. TE_origin_structure=tandem_LTR_RT is added when the array " +
    "also overlaps an LTR_RT_TR container (copies of complete LTR-RTs sharing LTRs). A satellite built from only part of an " +
    "element, e.g. an LTR fragment, never triggers this rule and stays a plain satellite.");
}

// ============================================================ slide 4
{
  const sl = caseSlide("LTR-RT inserted into a satellite: the element wins",
    "Domains occur once, not once per period, so this is an insertion, not a TE-derived tandem. The satellite is split around it.");
  const X = axis(2.35, 12.83);
  const ins = [36, 68], s = [2, 98];
  let y = 2.2;
  track(sl, X, y); rowLabel(sl, 2.15, y, "TideCluster", "one array");
  satellite(sl, X, ...s, y, { n: 40 });
  y = 3.1; track(sl, X, y); rowLabel(sl, 2.15, y, "DANTE_LTR", "complete element");
  ltrElement(sl, X, ...ins, y, { lw: 3 });
  [0, 1, 2, 3].forEach((j) => { const a = ins[0] + 5 + j * 6.2; domain(sl, X, a, a + 5, y, K.ltr, { h: 0.16, fill: "FFFFFF", line: "FFFFFF", ft: 30 }); });
  downArrow(sl, X(50), 3.6, 4.05);
  hline(sl, 0.5, 12.83, 4.17);
  const L1 = 4.8;
  track(sl, X, L1); rowLabel(sl, 2.15, L1, "Level 1", "unified");
  satellite(sl, X, s[0], ins[0], L1, { n: 14 });
  ltrElement(sl, X, ...ins, L1, { lw: 3 });
  satellite(sl, X, ins[1], s[1], L1, { n: 12 });
  label(sl, X(s[0]), L1 + 0.3, W(X, s[0], ins[0]), 0.28, "TRC_<n>", { fontSize: 12, color: K.tan, align: "center" });
  label(sl, X(ins[1]), L1 + 0.3, W(X, ins[1], s[1]), 0.28, "TRC_<n>", { fontSize: 12, color: K.tan, align: "center" });
  label(sl, 0.5, 5.65, 12.3, 1.2,
    "The same rule applies to any structural element (TIR, LINE) or protein domain inside an array.",
    { fontSize: 14, color: K.mute, valign: "top" });
  gffCard(sl, 9.2, 1.75, 3.63, 3.6, [
    { head: "Level 1, three features" }, { code: "TideCluster  Name=TRC_<n>" }, { code: "DANTE_LTR    (the element)" }, { code: "TideCluster  Name=TRC_<n>" },
    { head: "Why not TE_origin" }, { text: "Only one copy, and its domains do not repeat along the array, so the domain-rhythm test fails." },
    { head: "Rule" }, { text: "Plain tier priority: tier 1 beats tier 3." },
  ]);
  sl.addNotes(
    "This is the contrast to the previous slide. The array contains a complete element, but its domains appear only " +
    "once, not with the tandem period, so it is an insertion. Normal tier priority applies: the element keeps its span and " +
    "the satellite keeps the rest, as two pieces that both carry the same TRC name.");
}

// ============================================================ slide 5
{
  const sl = caseSlide("Satellite-library hit vs TE hit: tandem evidence decides",
    "A tier-4 hit outranks a TE hit only where TideHunter confirms a real tandem. Otherwise the TE hit wins.");
  const X = axis(2.35, 12.83);
  const L = [3, 47], R = [53, 97];
  label(sl, X(L[0]), 1.6, W(X, ...L), 0.32, "tandem confirmed", { bold: true, align: "center", fontSize: 15 });
  label(sl, X(R[0]), 1.6, W(X, ...R), 0.32, "no tandem evidence", { bold: true, align: "center", fontSize: 15 });
  sl.addShape(S.line, { objectName: nm("divider"), x: X(50), y: 1.6, w: 0, h: 3.65, line: { color: K.rule, width: 1, dashType: "dash" } });
  let y = 2.3;
  track(sl, X, y); rowLabel(sl, 2.15, y, "RepeatMasker", "TideCluster library");
  satellite(sl, X, L[0] + 3, L[1] - 3, y, { n: 14, h: 0.26 });
  satellite(sl, X, 64, 86, y, { n: 7, h: 0.26 });
  y = 3.0; track(sl, X, y); rowLabel(sl, 2.15, y, "RepeatMasker", "TE library");
  hit(sl, X, L[0] + 12, L[1] - 14, y, K.ltr, false);
  hit(sl, X, R[0] + 2, R[1] - 2, y, K.ltr, false);
  y = 3.7; track(sl, X, y); rowLabel(sl, 2.15, y, "TideHunter", "raw tandems");
  satellite(sl, X, L[0] + 3, L[1] - 3, y, { n: 14, h: 0.22, open: true });
  label(sl, X(R[0]), y - 0.16, W(X, ...R), 0.32, "none", { align: "center", color: K.mute, fontSize: 12, italic: true });
  downArrow(sl, X(25), 4.15, 4.55); downArrow(sl, X(75), 4.15, 4.55);
  hline(sl, 0.5, 12.83, 4.67);
  const L1 = 5.25;
  track(sl, X, L1); rowLabel(sl, 2.15, L1, "Level 1", "unified");
  satellite(sl, X, L[0] + 3, L[1] - 3, L1, { n: 14, h: 0.26 });
  hit(sl, X, R[0] + 2, R[1] - 2, L1, K.ltr, false);
  label(sl, 0.5, 5.75, 12.3, 1.2,
    "Right: a short, AT-rich satellite consensus can match inside a degraded TE copy. Without real tandem structure there, the TE call is kept.",
    { fontSize: 14, color: K.mute, valign: "top" });
  gffCard(sl, 9.2, 1.75, 3.63, 3.85, [
    { head: "Left: satellite wins" }, { code: "source  TideCluster_RM" }, { code: "source_tier=4" },
    { head: "Right: TE hit wins" }, { code: "source  RepeatMasker" }, { code: "source_tier=5" },
    { text: "Any part of the satellite hit outside the TE still fills free sequence." },
    { head: "Config" }, { code: "rm_tc_tandem_gate: True" }, { text: "(default)" },
  ]);
  sl.addNotes(
    "Tier 4 is RepeatMasker run with the TideCluster consensus library. Short satellite consensi can produce hits inside " +
    "TEs that are not tandem at all. So a tier-4 hit keeps its priority over a tier-5 TE hit only where raw TideHunter output " +
    "shows tandem structure. Without that support it is demoted below the TE hit: the TE keeps the contested span and the " +
    "satellite hit keeps only sequence nobody else claimed.");
}

(async () => {
  await pres.writeFile({ fileName: OUT });
  if (process.env.PPTX_SKILL) {
    const { applyTheme } = require(process.env.PPTX_SKILL + "/scripts/apply_theme.js");
    await applyTheme(OUT, THEME);
  }
  console.log("wrote", OUT);
})();
