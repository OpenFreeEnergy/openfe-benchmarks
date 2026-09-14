"""Focused tests for prepare_metadata_submission_v2 summary and charge labeling."""

from openfe_benchmarks.scripts import prepare_metadata_submission as pms_v2
from openfe_benchmarks.results import BenchmarkResults
from openfe_benchmarks.results._benchmark_results import Archive


def test_rbfe_summary_includes_separate_ff_and_charge_clauses() -> None:
    """RBFE summary should keep separate protein/solvent and ligand/cofactor wording."""
    metadata = pms_v2._Metadata(mode="rbfe", network_mode="alchemicalarchive")
    metadata.forcefield = [
        (("ff14SB", "lipid17_merged", "opc3_standard", "phosaa10"), ["edge1"])
    ]
    metadata.small_molecule_forcefield = [("openff-2.3.0", ["edge1"])]
    metadata.partial_charges = [("nagl_openff-gnn-am1bcc-1.0.0.pt", ["edge1"])]
    metadata.system_order = [("jacs_set", "tyk2")]
    metadata.n_transformations = 42
    metadata.systems[("jacs_set", "tyk2")] = pms_v2._SystemRecord(
        system_group="jacs_set",
        system_name="tyk2",
        network_key="AlchemicalNetwork-abc",
        ligands={"ligA", "ligB"},
        proteins={"tyk2"},
        cofactors=set(),
        solvents={"solvent"},
    )

    summary = pms_v2._build_content_summary(
        metadata=metadata,
        mode_spec=pms_v2._mode_spec("rbfe"),
        used_alchemiscale=False,
    )

    assert "for proteins and solvents" in summary
    assert "for ligands, solutes, and cofactors" in summary
    assert "openff-2.3.0" in summary
    assert "nagl_openff-gnn-am1bcc-1.0.0.pt" in summary


def test_normalize_charge_method_from_provenance_prefers_nagl_model() -> None:
    """NAGL provenance should include the model filename in the output tag."""
    provenance = {
        "charge_method": "NAGL",
        "nagl_model": "/path/to/openff-gnn-am1bcc-1.0.0.pt",
    }

    assert (
        pms_v2._normalize_charge_method_from_provenance(provenance)
        == "nagl_openff-gnn-am1bcc-1.0.0.pt"
    )


def test_make_zenodo_description_includes_requested_detail_sections() -> None:
    """Zenodo markdown should include the requested detailed metadata sections."""
    benchmark_results = BenchmarkResults(
        submission_id="2026-09-09-example",
        title="OpenFE RBFE Example",
        summary="Example summary.",
        tags=["rbfe", "openfe"],
        calculation_type="rbfe",
        authors=[{"name": "Example Author"}],
        date="2026-09-09",
        results="computational_results.json",
        archive=Archive(doi="10.1234/example", archive_provider="zenodo"),
        license="CC-BY-4.0",
        openfe_version="1.0.0",
        openmm_version="8.0.0",
        openff_toolkit_version="0.10.0",
        partial_charges="am1bcc_at",
        benchmark_data={"source_repository": "https://example.invalid/repo"},
        protocol_settings=[],
    )

    markdown = pms_v2._make_zenodo_description(
        benchmark_results=benchmark_results,
        network_mode="alchemicalarchive",
        used_alchemiscale=True,
    )

    assert "\n        ##" not in markdown
    assert "## Submission Snapshot" in markdown
    assert "## Systems Covered" in markdown
    assert "## Repository Reference" in markdown
    assert "## Software Versions" in markdown
    assert "## Recommended Descriptors" in markdown
    assert "## BenchmarkData Provenance" in markdown
    assert "## Protocol Settings" in markdown


def test_build_protocol_payload_uses_multiline_bullets_for_applies_notes() -> None:
    """Protocol notes should render each applies entry on its own line."""
    protocol_settings = [
        (
            {
                "protocol": "RelativeHybridTopologyProtocol",
                "protocol_library": "openfe",
                "notes": "",
            },
            [
                "edge-a",
                "edge-b",
            ],
        )
    ]

    payload = pms_v2._build_protocol_payload(protocol_settings)
    notes = payload[0]["notes"]

    assert isinstance(notes, str)
    assert "Applies to 2 edges:" in notes
    assert "\n- edge-a" in notes
    assert "\n- edge-b" in notes


def test_partial_charge_extraction_handles_small_molecule_type_name_fallback() -> None:
    """Charge provenance extraction should work when class identity differs but name matches."""

    FakeSmallMoleculeComponent = type(
        "SmallMoleculeComponent",
        (),
        {
            "__init__": lambda self, payload: setattr(self, "_payload", payload),
            "to_dict": lambda self: self._payload,
        },
    )

    class FakeChemicalSystem:
        def __init__(self, components):
            self.components = components

    class FakeTransformation:
        def __init__(self):
            payload = {
                "molprops": {
                    "charge_provenance": (
                        '{"charge_method": "NAGL", '
                        '"nagl_model": "/tmp/openff-gnn-am1bcc-1.0.0.pt"}'
                    )
                }
            }
            ligand = FakeSmallMoleculeComponent(payload)
            self.stateA = FakeChemicalSystem({"ligand": ligand})
            self.stateB = FakeChemicalSystem({"ligand": ligand})

    method = pms_v2._partial_charge_from_transformation(
        FakeTransformation(),
        pms_v2._mode_spec("rbfe"),
    )

    assert method == "nagl_openff-gnn-am1bcc-1.0.0.pt"
