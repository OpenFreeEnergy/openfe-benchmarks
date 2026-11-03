"use strict";

const PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#ad494a", "#637939", "#8c6d31", "#7b4173", "#3182bd"];
const SYMBOLS = ["circle", "square", "diamond", "triangle-up", "cross", "x", "star", "hexagon", "pentagon", "triangle-down"];
const STATS = ["RMSE", "MUE", "R2", "rho", "KTAU"];
const STAT_LABEL = { RMSE: "RMSE", MUE: "MUE", R2: "R²", rho: "Spearman ρ", KTAU: "Kendall τ" };
const ALL_KEY = "__all__";

const $ = (id) => document.getElementById(id);
const state = {
  index: null, stats: null, subs: new Map(), color: new Map(),
  tags: new Set(), tagAny: false, selected: new Set(),
  systemCache: new Map(), ligCache: new Map(), mapCache: new Map(),
  current: [], // plotted points per trace, used for click handling
  all: [],
};

const fmt = (x, n = 2) => (x === null || x === undefined || Number.isNaN(x) ? "–" : Number(x).toFixed(n));
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

async function getJSON(path) {
  const r = await fetch(path);
  if (!r.ok) throw new Error(`${path}: ${r.status}`);
  return r.json();
}

async function init() {
  try {
    [state.index, state.stats] = await Promise.all([getJSON("data/index.json"), getJSON("data/stats.json").catch(() => null)]);
  } catch (e) {
    $("status").textContent = "Could not load data bundle: " + e.message;
    return;
  }
  state.index.submissions.sort((a, b) => b.date.localeCompare(a.date) || a.id.localeCompare(b.id));
  state.index.submissions.forEach((s, i) => {
    state.subs.set(s.id, s);
    state.color.set(s.id, PALETTE[i % PALETTE.length]);
  });
  const params = new URLSearchParams(location.hash.slice(1));
  (params.get("subs") || "").split(",").filter((x) => state.subs.has(x)).forEach((x) => state.selected.add(x));
  (params.get("tags") || "").split(",").filter(Boolean).forEach((x) => state.tags.add(x));
  if (!state.selected.size) state.index.submissions.slice(0, 2).forEach((s) => state.selected.add(s.id));
  if (params.get("kind")) $("kind").value = params.get("kind");

  renderTags();
  renderSubs();
  for (const id of ["system", "kind", "common", "outlier", "errbars", "colorby", "f-overlap", "f-mbar"]) {
    $(id).addEventListener("change", () => { if (id === "kind") syncKindOptions(); refresh(); });
  }
  $("tag-search").addEventListener("input", renderTags);
  $("tag-any").addEventListener("change", () => { state.tagAny = $("tag-any").checked; renderSubs(); });
  $("tag-clear").addEventListener("click", () => { state.tags.clear(); renderTags(); renderSubs(); });
  $("sub-all").addEventListener("click", () => { visibleSubs().forEach((s) => state.selected.add(s.id)); renderSubs(); });
  $("sub-none").addEventListener("click", () => { state.selected.clear(); renderSubs(); });
  if (params.get("system")) state.wantedSystem = params.get("system");
  syncKindOptions();
  renderSystems();
  refresh();
}

function syncKindOptions() { $("ddg-opts").style.display = $("kind").value === "ddg" ? "flex" : "none"; }

// ---------------------------------------------------------------- sidebar

function allTags() {
  const c = new Map();
  state.index.submissions.forEach((s) => s.tags.forEach((t) => c.set(t, (c.get(t) || 0) + 1)));
  return [...c.entries()].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
}

function renderTags() {
  const q = $("tag-search").value.toLowerCase();
  const box = $("tag-list");
  box.innerHTML = "";
  allTags().filter(([t]) => state.tags.has(t) || t.toLowerCase().includes(q)).forEach(([t, n]) => {
    const el = document.createElement("span");
    el.className = "chip" + (state.tags.has(t) ? " on" : "");
    el.textContent = `${t} (${n})`;
    el.onclick = () => { state.tags.has(t) ? state.tags.delete(t) : state.tags.add(t); renderTags(); renderSubs(); };
    box.appendChild(el);
  });
}

function visibleSubs() {
  const tags = [...state.tags];
  return state.index.submissions.filter((s) => !tags.length ||
    (state.tagAny ? tags.some((t) => s.tags.includes(t)) : tags.every((t) => s.tags.includes(t))));
}

function renderSubs() {
  const box = $("sub-list");
  box.innerHTML = "";
  visibleSubs().forEach((s) => {
    const el = document.createElement("label");
    el.className = "sub-item";
    el.innerHTML = `<input type="checkbox" ${state.selected.has(s.id) ? "checked" : ""}>
      <span class="swatch" style="background:${state.color.get(s.id)}"></span>
      <span><span class="id">${esc(s.id)}</span><br>
      <span class="meta">${esc(s.calculation_type)} · ${esc(s.date)} · ${esc(s.small_molecule_forcefield || "")} ${esc(s.partial_charges || "")}</span></span>`;
    el.title = s.summary;
    el.querySelector("input").onchange = (e) => {
      e.target.checked ? state.selected.add(s.id) : state.selected.delete(s.id);
      renderSystems();
      refresh();
    };
    box.appendChild(el);
  });
  renderSystems();
  refresh();
}

function renderSystems() {
  const sel = $("system");
  const prev = state.wantedSystem || sel.value;
  state.wantedSystem = null;
  const chosen = [...state.selected];
  const systems = state.index.systems.filter((s) => chosen.some((id) => state.subs.get(id).systems.includes(s.key)));
  sel.innerHTML = "";
  if (!systems.length) return;
  sel.add(new Option(`All systems pooled (${systems.length})`, ALL_KEY));
  const groups = new Map();
  systems.forEach((s) => { if (!groups.has(s.group)) groups.set(s.group, []); groups.get(s.group).push(s); });
  for (const [g, items] of groups) {
    const og = document.createElement("optgroup");
    og.label = g;
    items.forEach((s) => og.appendChild(new Option(`${s.name} (${chosen.filter((id) => state.subs.get(id).systems.includes(s.key)).length}/${chosen.length})`, s.key)));
    sel.appendChild(og);
  }
  sel.value = [...sel.options].some((o) => o.value === prev) ? prev : sel.options[Math.min(1, sel.options.length - 1)].value;
}

// ---------------------------------------------------------------- data

async function loadSystem(key) {
  if (!state.systemCache.has(key)) state.systemCache.set(key, getJSON(`data/systems/${key}.json`));
  return state.systemCache.get(key);
}
const loadLigands = (key) => {
  if (!state.ligCache.has(key)) state.ligCache.set(key, getJSON(`data/structures/${key}.ligands.json`).catch(() => ({})));
  return state.ligCache.get(key);
};
const loadMappings = (key) => {
  if (!state.mapCache.has(key)) state.mapCache.set(key, getJSON(`data/structures/${key}.mappings.json`).catch(() => ({})));
  return state.mapCache.get(key);
};

const mean = (xs) => xs.reduce((a, b) => a + b, 0) / xs.length;

// Build plotted points for one submission in one system.
function buildPoints(sid, sys, kind) {
  const sub = sys.submissions[sid];
  if (!sub) return [];
  const asfe = state.subs.get(sid).calculation_type === "asfe";
  const pts = [];
  if (kind === "dg") {
    sub.dg.forEach((r) => {
      const e = sys.exp[r.ligand];
      if (e && r.dg !== null) pts.push({ sid, sysKey: sys.key, id: r.ligand, key: r.ligand, label: r.ligand, exp: e[0], expUnc: e[1], pred: r.dg, unc: r.unc || 0 });
    });
    if (pts.length && !asfe) {
      const shift = mean(pts.map((p) => p.exp - p.pred));
      pts.forEach((p) => { p.pred += shift; });
    }
  } else {
    sub.ddg.forEach((r) => {
      const ea = sys.exp[r.a], eb = sys.exp[r.b];
      if (!ea || !eb || r.ddg === null) return;
      let [a, b, exp, pred] = [r.a, r.b, eb[0] - ea[0], r.ddg];
      if (a > b) { [a, b, exp, pred] = [b, a, -exp, -pred]; }
      pts.push({ sid, sysKey: sys.key, a, b, key: `${a}|${b}`, label: `${a} → ${b}`, exp, expUnc: Math.hypot(ea[1], eb[1]), pred, unc: r.unc || 0, qc: r.qc, flipped: r.a !== a });
    });
  }
  pts.forEach((p) => { p.err = p.pred - p.exp; });
  return pts;
}

async function gatherPoints(kind) {
  const key = $("system").value;
  const sids = [...state.selected];
  const keys = key === ALL_KEY
    ? state.index.systems.filter((s) => sids.some((id) => state.subs.get(id).systems.includes(s.key))).map((s) => s.key)
    : [key];
  const systems = await Promise.all(keys.map(loadSystem));
  const bySub = new Map(sids.map((id) => [id, []]));
  systems.forEach((sys, i) => {
    sys.key = keys[i];
    sids.forEach((id) => bySub.get(id).push(...buildPoints(id, sys, kind)));
  });
  return bySub;
}

function applyFilters(bySub, kind) {
  const common = $("common").checked && bySub.size > 1;
  let out = new Map(bySub);
  if (common) {
    const keyOf = (p) => `${p.sysKey}#${p.key}`;
    const sets = [...bySub.values()].map((pts) => new Set(pts.map(keyOf)));
    const inter = [...sets[0]].filter((k) => sets.every((s) => s.has(k)));
    const keep = new Set(inter);
    out = new Map([...bySub].map(([id, pts]) => [id, pts.filter((p) => keep.has(keyOf(p)))]));
  }
  if (kind === "ddg") {
    const ov = parseFloat($("f-overlap").value), mb = parseFloat($("f-mbar").value);
    out = new Map([...out].map(([id, pts]) => [id, pts.filter((p) =>
      (Number.isNaN(ov) || (p.qc.min_overlap ?? 1) >= ov) && (Number.isNaN(mb) || (p.qc.max_mbar_error ?? 0) <= mb))]));
  }
  return out;
}

// ---------------------------------------------------------------- rendering

let refreshToken = 0;
async function refresh() {
  const token = ++refreshToken;
  $("status").textContent = "";
  const kind = $("kind").value;
  if (!state.selected.size || !$("system").value) {
    Plotly.purge("scatter");
    $("stats-table").innerHTML = ""; $("stats-plots").innerHTML = ""; $("outliers").innerHTML = "";
    $("status").textContent = "Select at least one submission.";
    return;
  }
  let bySub;
  try { bySub = await gatherPoints(kind); } catch (e) { $("status").textContent = e.message; return; }
  if (token !== refreshToken) return;
  const shown = applyFilters(bySub, kind);
  const missing = [...bySub].filter(([, p]) => !p.length).map(([id]) => id);
  const flat = [...shown.values()].flat();
  state.all = flat;
  location.hash = new URLSearchParams({ subs: [...state.selected].join(","), tags: [...state.tags].join(","), system: $("system").value, kind }).toString();
  $("note").textContent = (missing.length ? `No ${kind.toUpperCase()} vs experiment data for: ${missing.join(", ")}. ` : "") +
    (kind === "dg" ? "RBFE DG values are shifted so their mean matches experiment (per system); ASFE values are not shifted. " : "") +
    "Edges are shown in a canonical (alphabetical) orientation.";
  drawScatter(shown, kind);
  drawStats(shown, kind);
  drawOutliers(flat, kind);
}

function hoverText(p, kind) {
  let t = `<b>${esc(p.label)}</b><br>${esc(p.sid)}<br>exp ${fmt(p.exp)} · pred ${fmt(p.pred)} ± ${fmt(p.unc)}<br>error ${fmt(p.err)}`;
  if (p.sysKey && $("system").value === ALL_KEY) t += `<br>${esc(p.sysKey)}`;
  if (kind === "ddg" && p.qc) t += `<br>min overlap ${fmt(p.qc.min_overlap, 3)} · max MBAR err ${fmt(p.qc.max_mbar_error, 2)}`;
  return t;
}

function metric(p, name) { return name === "abs_error" ? Math.abs(p.err) : p.qc?.[name]; }

function drawScatter(shown, kind) {
  const colorBy = kind === "ddg" ? $("colorby").value : "";
  const traces = [];
  const all = [...shown.values()].flat();
  const lim = all.length ? [Math.min(...all.map((p) => Math.min(p.exp, p.pred))), Math.max(...all.map((p) => Math.max(p.exp, p.pred)))] : [-1, 1];
  const pad = (lim[1] - lim[0]) * 0.08 + 0.5;
  const lo = lim[0] - pad, hi = lim[1] + pad;
  const band = (d, dash) => ({ x: [lo, hi], y: [lo + d, hi + d], mode: "lines", line: { color: "#999", dash, width: 1 }, hoverinfo: "skip", showlegend: false });
  traces.push(band(0, "solid"), band(1, "dot"), band(-1, "dot"), band(2, "dash"), band(-2, "dash"));
  state.current = [];
  let i = 0, cbShown = false;
  for (const [sid, pts] of shown) {
    if (!pts.length) { i++; continue; }
    state.current.push(pts);
    const marker = { size: 8, color: state.color.get(sid), symbol: "circle", opacity: 0.85, line: { width: 0.5, color: "#fff" } };
    if (colorBy) {
      marker.color = pts.map((p) => metric(p, colorBy) ?? NaN);
      marker.symbol = SYMBOLS[i % SYMBOLS.length];
      marker.colorscale = colorBy === "min_overlap" || colorBy === "min_replica_mixing" ? "RdYlGn" : "YlOrRd";
      marker.showscale = !cbShown;
      marker.colorbar = { title: colorBy.replace(/_/g, " "), thickness: 12 };
      cbShown = true;
    }
    traces.push({
      x: pts.map((p) => p.exp), y: pts.map((p) => p.pred), mode: "markers", type: "scatter", name: `${sid} (n=${pts.length})`,
      marker, text: pts.map((p) => hoverText(p, kind)), hoverinfo: "text",
      error_y: $("errbars").checked ? { type: "data", array: pts.map((p) => p.unc), visible: true, thickness: 0.8, color: state.color.get(sid) + "88" } : undefined,
      meta: state.current.length - 1,
    });
    i++;
  }
  const unit = kind === "dg" ? "ΔG" : "ΔΔG";
  Plotly.react("scatter", traces, {
    margin: { t: 20, l: 60, r: 20, b: 50 }, hovermode: "closest",
    xaxis: { title: `Experimental ${unit} (kcal/mol)`, range: [lo, hi], zeroline: false, scaleanchor: "y" },
    yaxis: { title: `Calculated ${unit} (kcal/mol)`, range: [lo, hi], zeroline: false },
    legend: { orientation: "h", y: -0.2 },
  }, { responsive: true });
  const el = $("scatter");
  el.removeAllListeners?.("plotly_click");
  el.on("plotly_click", (ev) => {
    const pt = ev.points[0];
    if (pt.data.meta === undefined) return;
    showDetail(state.current[pt.data.meta][pt.pointIndex], kind);
  });
}

// ---------------------------------------------------------------- stats

function lookupStats(shown, kind) {
  const sids = [...shown.keys()];
  const sysKey = $("system").value;
  const sorted = [...sids].sort();
  if ($("common").checked && sids.length === 2) {
    const e = state.stats?.pair?.[`${sorted[0]}|${sorted[1]}|${sysKey}|${kind}`];
    if (e) return { rows: e, src: "common subset of the two selected submissions (cinnabar bootstrap)" };
  }
  const rows = {};
  sids.forEach((id) => { const s = state.stats?.single?.[`${id}|${sysKey}|${kind}`]; if (s) rows[id] = s; });
  const extra = $("common").checked && sids.length > 2 ? " – common-subset statistics are precomputed for pairs only" : "";
  return { rows, src: `each submission's full set with experimental data (cinnabar bootstrap, 95% CI)${extra}` };
}

function drawStats(shown, kind) {
  const { rows, src } = lookupStats(shown, kind);
  $("stats-src").textContent = "· " + src;
  const ids = Object.keys(rows);
  if (!ids.length) { $("stats-table").innerHTML = "<p class='hint'>No precomputed statistics for this selection.</p>"; $("stats-plots").innerHTML = ""; return; }
  let h = "<table><tr><th>Submission</th><th>N</th>" + STATS.map((s) => `<th>${STAT_LABEL[s]}</th>`).join("") + "</tr>";
  ids.forEach((id) => {
    const r = rows[id];
    h += `<tr><td><span class="swatch" style="display:inline-block;width:9px;height:9px;background:${state.color.get(id)};margin-right:5px"></span>${esc(id)}</td><td>${r.n}</td>` +
      STATS.map((s) => r[s] ? `<td>${fmt(r[s].mle)} <span class="note">[${fmt(r[s].low)}, ${fmt(r[s].high)}]</span></td>` : "<td>–</td>").join("") + "</tr>";
  });
  $("stats-table").innerHTML = h + "</table>";

  const box = $("stats-plots");
  box.innerHTML = "";
  STATS.forEach((s) => {
    const d = document.createElement("div");
    box.appendChild(d);
    const have = ids.filter((id) => rows[id][s]);
    Plotly.newPlot(d, [{
      type: "scatter", mode: "markers", orientation: "h",
      y: have.map((id) => id), x: have.map((id) => rows[id][s].mle),
      error_x: { type: "data", symmetric: false, array: have.map((id) => rows[id][s].high - rows[id][s].mle), arrayminus: have.map((id) => rows[id][s].mle - rows[id][s].low) },
      marker: { size: 9, color: have.map((id) => state.color.get(id)) },
      hovertemplate: "%{y}<br>%{x:.2f}<extra></extra>",
    }], {
      title: { text: STAT_LABEL[s], font: { size: 13 } }, margin: { t: 30, l: 10, r: 10, b: 30 },
      yaxis: { showticklabels: false, autorange: "reversed" }, xaxis: { zeroline: false },
    }, { responsive: true, displayModeBar: false });
  });
}

// ---------------------------------------------------------------- outliers & detail

async function drawOutliers(points, kind) {
  const thr = parseFloat($("outlier").value) || 0;
  const rows = points.filter((p) => Math.abs(p.err) > thr).sort((a, b) => Math.abs(b.err) - Math.abs(a.err)).slice(0, 25);
  const total = points.filter((p) => Math.abs(p.err) > thr).length;
  if (!rows.length) { $("outliers").innerHTML = "<p class='hint'>No points above the threshold.</p>"; return; }
  const sysKey = $("system").value;
  const ligs = kind === "dg" && sysKey !== ALL_KEY ? await loadLigands(sysKey) : {};
  let h = `<p class="note">${total} points above threshold${total > rows.length ? `, showing top ${rows.length}` : ""}</p><table><tr>` +
    (kind === "dg" && sysKey !== ALL_KEY ? "<th>Structure</th>" : "") +
    `<th>${kind === "dg" ? "Ligand" : "Edge"}</th><th>Submission</th>${sysKey === ALL_KEY ? "<th>System</th>" : ""}<th>Exp</th><th>Calc</th><th>Error</th></tr>`;
  rows.forEach((p, i) => {
    h += `<tr class="clickable" data-i="${i}">` + (kind === "dg" && sysKey !== ALL_KEY ? `<td>${ligs[p.id] || ""}</td>` : "") +
      `<td>${esc(p.label)}</td><td style="text-align:left">${esc(p.sid)}</td>${sysKey === ALL_KEY ? `<td>${esc(p.sysKey)}</td>` : ""}<td>${fmt(p.exp)}</td><td>${fmt(p.pred)}</td><td>${fmt(p.err)}</td></tr>`;
  });
  const box = $("outliers");
  box.innerHTML = h + "</table>";
  box.querySelectorAll("tr.clickable").forEach((tr) => tr.onclick = () => showDetail(rows[+tr.dataset.i], kind));
}

const list = (xs, n = 2) => (xs || []).map((x) => fmt(x, n)).join(", ");

async function showDetail(p, kind) {
  const box = $("detail");
  box.innerHTML = "<p class='hint'>Loading…</p>";
  const ligs = await loadLigands(p.sysKey);
  const struct = (id) => `<div><div class="note">${esc(id)}</div>${ligs[id] || "<span class='hint'>no structure</span>"}</div>`;
  let h = `<h3 style="margin-top:0">${esc(p.label)}</h3><div class="note">${esc(p.sid)} · ${esc(p.sysKey)}</div>
    <table><tr><td>Experiment</td><td>${fmt(p.exp)} ± ${fmt(p.expUnc)}</td></tr>
    <tr><td>Calculated</td><td>${fmt(p.pred)} ± ${fmt(p.unc)}</td></tr><tr><td>Error</td><td>${fmt(p.err)}</td></tr></table>`;
  if (kind === "dg") {
    h += `<div class="pair">${struct(p.id)}</div>`;
  } else {
    const maps = await loadMappings(p.sysKey);
    const fwd = maps[`${p.a}|${p.b}`], rev = maps[`${p.b}|${p.a}`];
    const m = fwd || rev;
    // Show the ligands in the same left/right order as the mapping depiction
    const [first, second] = !fwd && rev ? [p.b, p.a] : [p.a, p.b];
    h += `<div class="pair">${struct(first)}${struct(second)}</div>`;
    h += `<h3>Atom mapping</h3>${m || "<span class='hint'>no mapping available for this edge</span>"}`;
    if (p.flipped) h += `<p class="note">Submission orientation is ${esc(p.b)} → ${esc(p.a)}; values above are sign-flipped to the canonical orientation (per-repeat ΔG values below are in the submitted orientation).</p>`;
    const q = p.qc;
    const rowq = (name, c, s, n = 2) => `<tr><td>${name}</td><td>${list(c, n)}</td><td>${list(s, n)}</td></tr>`;
    h += `<h3>Simulation QC <span class="note">(repeats)</span></h3><table><tr><th></th><th>complex</th><th>solvent</th></tr>
      ${rowq("ΔG (kcal/mol)", q.dgs_complex, q.dgs_solvent)}
      ${rowq("MBAR error", q.complex_mbar_errors, q.solvent_mbar_errors)}
      ${rowq("min MBAR overlap", q.complex_overlap, q.solvent_overlap, 3)}
      ${rowq("min replica mixing", q.complex_replica_mixing, q.solvent_replica_mixing, 3)}
      ${rowq("equil. iterations", q.complex_equil_iters, q.solvent_equil_iters, 0)}
      ${rowq("production iterations", q.complex_prod_iters, q.solvent_prod_iters, 0)}</table>`;
  }
  box.innerHTML = h;
}

init();
