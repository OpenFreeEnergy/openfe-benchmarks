"""Tests for the dashboard build helpers (no results data required)."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import build_data as bd  # noqa: E402


def _sys(exp, dg=(), ddg=()):
    return {"exp": exp, "dg": list(dg), "ddg": list(ddg)}


def test_ddg_points_orientation_normalised():
    sd = _sys(
        {"a": (-8.0, 0.1), "b": (-9.0, 0.2)},
        ddg=[{"a": "b", "b": "a", "ddg": 1.5, "unc": 0.1, "qc": {}}],
    )
    pts = bd._ddg_points(sd)
    exp, _, pred, _ = pts[("a", "b")]
    # submitted b->a: exp ddg = -8 - (-9) = +1 ; flipped to canonical a->b gives -1 / -1.5
    assert exp == -1.0 and pred == -1.5


def test_dg_points_shift_rbfe_not_asfe():
    sd = _sys(
        {"a": (-8.0, 0.0), "b": (-10.0, 0.0)},
        dg=[{"ligand": "a", "dg": 1.0, "unc": 0.1}, {"ligand": "b", "dg": -1.0, "unc": 0.1}],
    )
    rbfe = bd._dg_points(sd, asfe=False)
    assert np.isclose(np.mean([v[2] - v[0] for v in rbfe.values()]), 0.0)
    asfe = bd._dg_points(sd, asfe=True)
    assert asfe["a"][2] == 1.0


def test_bootstrap_requires_three_points_and_returns_stats():
    assert bd._bootstrap([(0, 0, 0, 0)] * 2, 10) is None
    rng = np.random.default_rng(0)
    x = rng.normal(size=20)
    pts = [(v, 0.0, v + 0.1, 0.0) for v in x]
    res = bd._bootstrap(pts, 50)
    assert res["n"] == 20
    assert res["RMSE"]["mle"] < 0.2
    assert res["R2"]["mle"] > 0.9


def test_minify_svg_strips_xml_header():
    svg = "<?xml version='1.0'?>\n<svg>\n<!-- END OF HEADER -->\n<path style='fill:#000;fill-rule:evenodd'/></svg>"
    out = bd._minify_svg(svg)
    assert out.startswith("<svg>") and "evenodd" not in out
