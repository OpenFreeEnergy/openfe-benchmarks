"""Build the static data bundle for the benchmark comparison dashboard.

Writes JSON/SVG assets into ``dashboard/site/data``:

* ``index.json``                      submissions (metadata + tags) and systems
* ``systems/<group>__<name>.json``    per-system DG / DDG results incl. QC metrics
* ``structures/<sys>.ligands.json``   ligand 2D depictions (SVG)
* ``structures/<sys>.mappings.json``  atom-mapping depictions (SVG)
* ``stats.json``                      cinnabar bootstrap statistics

Usage: python dashboard/build_data.py [--nbootstrap 1000] [--no-structures] [--no-stats]
"""

from __future__ import annotations

import argparse
import itertools
import json
import logging
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from cinnabar.stats import bootstrap_statistic
from gufe import LigandNetwork
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

import openfe_benchmarks.results._benchmark_results as _br
from openfe_benchmarks.data import get_data_by_system_name
from openfe_benchmarks.results._benchmark_results import get_benchmark_results

logging.disable(logging.INFO)

SITE_DIR = Path(__file__).parent / "site"
DATA_DIR = SITE_DIR / "data"
STATS = ["RMSE", "MUE", "R2", "rho", "KTAU"]
ALL_KEY = "__all__"
KCAL = "kilocalorie_per_mole"


def _mag(x) -> float | None:
    """Convert a (possibly pint) quantity to a float in kcal/mol."""
    if x is None:
        return None
    if hasattr(x, "to"):
        try:
            return float(x.to(KCAL).magnitude)
        except Exception:
            return float(x.magnitude)
    if hasattr(x, "magnitude"):
        return float(x.magnitude)
    return float(x)


def _mags(xs) -> list[float | None]:
    return [None if v is None else round(_mag(v), 4) for v in (xs or [])]


def _r(x, n=4):
    """Round to ``n`` digits; non-finite values become None so the output stays valid JSON."""
    return None if x is None or not np.isfinite(x) else round(float(x), n)


def system_key(group: str, name: str) -> str:
    return f"{group}__{name}"


# ----------------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------------


def discover_submissions() -> list[str]:
    root = Path(_br.__file__).parent
    return sorted(p.parent.name for p in root.glob("*/submission.yaml"))


def _exp_from_femap(femap) -> dict[str, tuple[float, float]]:
    df = femap.get_absolute_dataframe()
    df = df[df["source"] == "experimental"]
    return {
        r["label"]: (float(r["DG (kcal/mol)"]), float(r["uncertainty (kcal/mol)"]))
        for _, r in df.iterrows()
    }


def load_submission(sid: str) -> dict | None:
    """Return per-system DG/DDG records for one submission."""
    try:
        res = get_benchmark_results(sid)
    except Exception as e:  # noqa: BLE001
        print(f"  skipping {sid}: {e}")
        return None

    raw = res.raw_results or {}
    systems: dict[tuple[str, str], dict] = defaultdict(lambda: {"dg": [], "ddg": [], "exp": {}})

    # experimental values, taken from the FEMaps built by the existing loaders
    for kind in ("dg", "ddg"):
        if kind not in raw:
            continue
        try:
            femaps = res.dg_femaps() if kind == "dg" else res.ddg_femaps()
        except Exception as e:  # noqa: BLE001
            print(f"  {sid}: could not build {kind} femaps: {e}")
            continue
        for key, femap in femaps.items():
            systems[key]["exp"].update(_exp_from_femap(femap))

    for row in raw.get("dg", []):
        key = (row["system_group"], row["system_name"])
        label = row.get("ligand") or row.get("label")
        systems[key]["dg"].append(
            {"ligand": label, "dg": _r(_mag(row["dg"])), "unc": _r(_mag(row.get("dg_uncertainty")))}
        )

    for row in raw.get("ddg", []):
        key = (row["system_group"], row["system_name"])
        c_ov = _mags(row.get("complex_smallest_mbar_overlaps"))
        s_ov = _mags(row.get("solvent_smallest_mbar_overlaps"))
        c_er = _mags(row.get("complex_mbar_errors"))
        s_er = _mags(row.get("solvent_mbar_errors"))
        c_rm = _mags(row.get("complex_smallest_replica_mixing"))
        s_rm = _mags(row.get("solvent_smallest_replica_mixing"))
        valid = lambda xs: [v for v in xs if v is not None]  # noqa: E731
        systems[key]["ddg"].append(
            {
                "a": row["ligand_a"],
                "b": row["ligand_b"],
                "ddg": _r(_mag(row["ddg"])),
                "unc": _r(_mag(row.get("ddg_uncertainty"))),
                "qc": {
                    "dgs_complex": _mags(row.get("dgs_complex")),
                    "dgs_solvent": _mags(row.get("dgs_solvent")),
                    "complex_mbar_errors": c_er,
                    "solvent_mbar_errors": s_er,
                    "complex_overlap": c_ov,
                    "solvent_overlap": s_ov,
                    "complex_replica_mixing": c_rm,
                    "solvent_replica_mixing": s_rm,
                    "complex_equil_iters": _mags(row.get("complex_equilibration_iterations")),
                    "solvent_equil_iters": _mags(row.get("solvent_equilibration_iterations")),
                    "complex_prod_iters": _mags(row.get("complex_production_iterations")),
                    "solvent_prod_iters": _mags(row.get("solvent_production_iterations")),
                    "min_overlap": min(valid(c_ov + s_ov), default=None),
                    "max_mbar_error": max(valid(c_er + s_er), default=None),
                    "min_replica_mixing": min(valid(c_rm + s_rm), default=None),
                },
            }
        )

    meta = {
        "id": sid,
        "title": res.title,
        "summary": (res.summary or "").strip(),
        "date": str(res.date),
        "calculation_type": res.calculation_type,
        "tags": list(res.tags),
        "openfe_version": getattr(res, "openfe_version", None),
        "forcefield": getattr(res, "forcefield", None),
        "small_molecule_forcefield": getattr(res, "small_molecule_forcefield", None),
        "partial_charges": getattr(res, "partial_charges", None),
        "mapper": getattr(res, "mapper", None),
        "authors": [a.get("name") for a in (getattr(res, "authors", None) or [])],
    }
    return {"meta": meta, "systems": dict(systems)}


# ----------------------------------------------------------------------------
# Statistics (cinnabar)
# ----------------------------------------------------------------------------


def _dg_points(sys_data: dict, asfe: bool) -> dict[str, tuple[float, float, float, float]]:
    """ligand -> (exp, exp_unc, pred, pred_unc); pred shifted to the experimental mean for RBFE."""
    pts = {}
    for r in sys_data["dg"]:
        e = sys_data["exp"].get(r["ligand"])
        if e is not None and r["dg"] is not None:
            pts[r["ligand"]] = [e[0], e[1], r["dg"], r["unc"] or 0.0]
    if pts and not asfe:
        shift = np.mean([v[0] - v[2] for v in pts.values()])
        for v in pts.values():
            v[2] += shift
    return {k: tuple(v) for k, v in pts.items()}


def _ddg_points(sys_data: dict) -> dict[tuple[str, str], tuple[float, float, float, float]]:
    """canonical (a, b) -> (exp, exp_unc, pred, pred_unc), orientation-normalised."""
    pts = {}
    for r in sys_data["ddg"]:
        ea, eb = sys_data["exp"].get(r["a"]), sys_data["exp"].get(r["b"])
        if ea is None or eb is None or r["ddg"] is None:
            continue
        exp, eunc = eb[0] - ea[0], float(np.hypot(ea[1], eb[1]))
        pred, unc = r["ddg"], r["unc"] or 0.0
        key = (r["a"], r["b"])
        if key[0] > key[1]:
            key, exp, pred = (key[1], key[0]), -exp, -pred
        pts[key] = (exp, eunc, pred, unc)
    return pts


def _bootstrap(points: list[tuple], nbootstrap: int, seed: int = 0) -> dict | None:
    if len(points) < 3:
        return None
    arr = np.array(points, dtype=float)
    out = {"n": len(points)}
    for stat in STATS:
        np.random.seed(seed)
        try:
            r = bootstrap_statistic(
                arr[:, 0], arr[:, 2], arr[:, 1], arr[:, 3], statistic=stat, nbootstrap=nbootstrap
            )
        except Exception:  # noqa: BLE001
            continue
        out[stat] = {k: _r(r[k]) for k in ("mle", "mean", "low", "high")}
    return out


def compute_stats(subs: dict[str, dict], nbootstrap: int) -> dict:
    """Own-set stats per submission and common-subset stats per pair, per system and pooled."""
    pts: dict[str, dict[str, dict]] = {}  # sid -> kind -> sysk -> points
    for sid, sub in subs.items():
        asfe = sub["meta"]["calculation_type"] == "asfe"
        pts[sid] = {"dg": {}, "ddg": {}}
        for (g, n), sd in sub["systems"].items():
            k = system_key(g, n)
            d = _dg_points(sd, asfe)
            if d:
                pts[sid]["dg"][k] = d
            e = _ddg_points(sd)
            if e:
                pts[sid]["ddg"][k] = e

    def pooled(kind, sid, keys=None):
        data = pts[sid][kind]
        return {(k, p): v for k, d in data.items() if keys is None or k in keys for p, v in d.items()}

    stats: dict = {"single": {}, "pair": {}}
    for sid in subs:
        for kind in ("dg", "ddg"):
            for k, d in pts[sid][kind].items():
                s = _bootstrap(list(d.values()), nbootstrap)
                if s:
                    stats["single"][f"{sid}|{k}|{kind}"] = s
            allp = pooled(kind, sid)
            s = _bootstrap(list(allp.values()), nbootstrap)
            if s:
                stats["single"][f"{sid}|{ALL_KEY}|{kind}"] = s

    for a, b in itertools.combinations(sorted(subs), 2):
        for kind in ("dg", "ddg"):
            shared_sys = set(pts[a][kind]) & set(pts[b][kind])
            common_all = {}
            for k in sorted(shared_sys):
                common = sorted(set(pts[a][kind][k]) & set(pts[b][kind][k]))
                if not common:
                    continue
                sa = _bootstrap([pts[a][kind][k][p] for p in common], nbootstrap)
                sb = _bootstrap([pts[b][kind][k][p] for p in common], nbootstrap)
                if sa and sb:
                    stats["pair"][f"{a}|{b}|{k}|{kind}"] = {a: sa, b: sb}
                common_all[k] = common
            if len(common_all) > 1:
                sa = _bootstrap([pts[a][kind][k][p] for k, c in common_all.items() for p in c], nbootstrap)
                sb = _bootstrap([pts[b][kind][k][p] for k, c in common_all.items() for p in c], nbootstrap)
                if sa and sb:
                    stats["pair"][f"{a}|{b}|{ALL_KEY}|{kind}"] = {a: sa, b: sb}
    return stats


# ----------------------------------------------------------------------------
# Structures
# ----------------------------------------------------------------------------


def _minify_svg(svg: str) -> str:
    """Drop the XML declaration and redundant default style attributes to shrink the SVG."""
    svg = re.sub(r"<\?xml[^>]*\?>\s*", "", svg)
    svg = svg.replace("<!-- END OF HEADER -->\n", "")
    for junk in (";fill-rule:evenodd", ";stroke-linecap:butt", ";stroke-linejoin:miter", ";stroke-opacity:1", ";fill-opacity:1"):
        svg = svg.replace(junk, "")
    return svg.strip()


def _svg(mol, w=260, h=200) -> str:
    d = rdMolDraw2D.MolDraw2DSVG(w, h)
    d.drawOptions().clearBackground = False
    d.drawOptions().useBWAtomPalette()
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return _minify_svg(d.GetDrawingText())


def _ligand_mols(bd) -> dict[str, Chem.Mol]:
    mols = {}
    path = bd.ligands.get("no_charges") if bd is not None else None
    if path is not None:
        for m in Chem.SDMolSupplier(str(path), removeHs=False):
            if m is not None and m.HasProp("_Name"):
                mols[m.GetProp("_Name")] = m
    return mols


def build_structures(system_ligands: dict, system_edges: dict) -> None:
    """system_ligands: sysk -> set(ligand ids); system_edges: sysk -> set((a, b))."""
    from gufe.visualization.mapping_visualization import draw_mapping

    out = DATA_DIR / "structures"
    out.mkdir(parents=True, exist_ok=True)
    for sysk, ligands in sorted(system_ligands.items()):
        group, name = sysk.split("__", 1)
        try:
            bd = get_data_by_system_name(group, name)
        except Exception as e:  # noqa: BLE001
            print(f"  no benchmark data for {sysk}: {e}")
            continue
        mols = _ligand_mols(bd)
        svgs = {}
        for lig in sorted(ligands):
            m = mols.get(lig)
            if m is None:
                continue
            m = Chem.RemoveHs(m)
            rdDepictor.Compute2DCoords(m)
            svgs[lig] = _svg(m)
        (out / f"{sysk}.ligands.json").write_text(json.dumps(svgs, separators=(",", ":")))

        edges = system_edges.get(sysk, set())
        if not edges or not bd.ligand_networks:
            continue
        mappings = {}
        net_path = next(iter(bd.ligand_networks.values()))
        try:
            net = LigandNetwork.from_json(net_path)
        except Exception as e:  # noqa: BLE001
            print(f"  cannot read network for {sysk}: {e}")
            continue
        wanted = {frozenset(e) for e in edges}
        for e in net.edges:
            pair = frozenset((e.componentA.name, e.componentB.name))
            if pair not in wanted:
                continue
            d2d = rdMolDraw2D.MolDraw2DSVG(640, 280, 320, 280)
            d2d.drawOptions().clearBackground = False
            try:
                draw_mapping(e.componentA_to_componentB, e.componentA.to_rdkit(), e.componentB.to_rdkit(), d2d)
                mappings[f"{e.componentA.name}|{e.componentB.name}"] = _minify_svg(d2d.GetDrawingText())
            except Exception as ex:  # noqa: BLE001
                print(f"  mapping draw failed {sysk} {pair}: {ex}")
        (out / f"{sysk}.mappings.json").write_text(json.dumps(mappings, separators=(",", ":")))
        print(f"  {sysk}: {len(svgs)} ligands, {len(mappings)}/{len(edges)} mappings")


# ----------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nbootstrap", type=int, default=1000)
    ap.add_argument("--no-structures", action="store_true")
    ap.add_argument("--no-stats", action="store_true")
    ap.add_argument("--submission", action="append", help="restrict to these submission ids")
    args = ap.parse_args()

    (DATA_DIR / "systems").mkdir(parents=True, exist_ok=True)
    subs: dict[str, dict] = {}
    for sid in args.submission or discover_submissions():
        print(f"loading {sid}")
        s = load_submission(sid)
        if s is not None:
            subs[sid] = s

    per_system: dict[str, dict] = defaultdict(lambda: {"submissions": {}})
    ligs, edges = defaultdict(set), defaultdict(set)
    for sid, sub in subs.items():
        for (g, n), sd in sub["systems"].items():
            k = system_key(g, n)
            ps = per_system[k]
            ps.update(group=g, name=n)
            ps.setdefault("exp", {}).update({lig: [_r(v[0]), _r(v[1])] for lig, v in sd["exp"].items()})
            ps["submissions"][sid] = {"dg": sd["dg"], "ddg": sd["ddg"]}
            ligs[k].update(r["ligand"] for r in sd["dg"])
            ligs[k].update(sd["exp"])
            for r in sd["ddg"]:
                edges[k].add((r["a"], r["b"]))
                ligs[k].update((r["a"], r["b"]))

    for k, ps in per_system.items():
        (DATA_DIR / "systems" / f"{k}.json").write_text(json.dumps(ps, separators=(",", ":")))

    index = {
        "submissions": [s["meta"] | {"systems": sorted(system_key(*key) for key in s["systems"])} for s in subs.values()],
        "systems": [
            {"key": k, "group": v["group"], "name": v["name"], "n_submissions": len(v["submissions"])}
            for k, v in sorted(per_system.items())
        ],
    }
    (DATA_DIR / "index.json").write_text(json.dumps(index, separators=(",", ":")))

    if not args.no_stats:
        print("computing bootstrap statistics with cinnabar")
        stats = compute_stats(subs, args.nbootstrap)
        (DATA_DIR / "stats.json").write_text(json.dumps(stats, separators=(",", ":")))
    if not args.no_structures:
        print("building structures")
        build_structures(ligs, edges)


if __name__ == "__main__":
    main()
