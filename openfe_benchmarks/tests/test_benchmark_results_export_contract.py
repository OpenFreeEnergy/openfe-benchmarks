"""Tests for BenchmarkResults submission export contract checks."""

from __future__ import annotations

import yaml

import openfe_benchmarks.results._benchmark_results as br_module
from openfe_benchmarks.results import BenchmarkResults
from openfe_benchmarks.results._benchmark_results import Archive, LiteralStr

_INTERNAL_NON_EXPORT_FIELDS = [
    "raw_results",
    "results_file",
    "submission_file",
    "_dg_femaps_cache",
    "_ddg_femaps_cache",
    "_dg_femaps_source",
    "_ddg_femaps_source",
]


def _make_minimal_benchmark_results() -> BenchmarkResults:
    """Create a minimal valid BenchmarkResults object for export tests."""
    return BenchmarkResults(
        submission_id="test-submission",
        title="Test Submission",
        summary="Test summary",
        tags=["test"],
        calculation_type="rbfe",
        authors=[{"name": "Test Author"}],
        date="2026-01-01",
        results="computational_results.json",
        archive=Archive(doi="10.1234/example", archive_provider="zenodo"),
        license="CC-BY-4.0",
        openfe_version="1.0.0",
        openmm_version="8.0.0",
        openff_toolkit_version="0.10.0",
        partial_charges="am1bcc",
        benchmark_data={"source_repository": "https://example.invalid/repo"},
        protocol_settings=[],
    )


def test_submission_field_classification_is_exhaustive() -> None:
    """All dataclass fields should be categorized as export or internal."""
    dataclass_fields = set(BenchmarkResults.__dataclass_fields__)
    categorized_fields = (
        set(br_module._SUBMISSION_EXPORT_FIELDS)
        | set(br_module._SUBMISSION_OPTIONAL_FIELDS)
        | set(_INTERNAL_NON_EXPORT_FIELDS)
    )

    assert dataclass_fields == categorized_fields


def test_to_submission_dict_includes_optional_fields_as_none() -> None:
    """Optional export fields are always emitted with explicit None when unset."""
    benchmark_results = _make_minimal_benchmark_results()

    exported = benchmark_results.to_submission_dict()

    for field_name in br_module._SUBMISSION_OPTIONAL_FIELDS:
        assert field_name in exported
        assert exported[field_name] is None


def test_to_submission_dict_excludes_internal_fields() -> None:
    """Internal runtime-only fields should never appear in exported metadata."""
    benchmark_results = _make_minimal_benchmark_results()

    exported = benchmark_results.to_submission_dict()

    for field_name in _INTERNAL_NON_EXPORT_FIELDS:
        assert field_name not in exported


def test_exported_dict_and_yaml_keys_match_declared_export_contract() -> None:
    """Exported dict/YAML keys must exactly match required+optional export lists."""
    benchmark_results = _make_minimal_benchmark_results()
    expected_keys = set(br_module._SUBMISSION_EXPORT_FIELDS) | set(
        br_module._SUBMISSION_OPTIONAL_FIELDS
    )

    exported_dict = benchmark_results.to_submission_dict()
    assert set(exported_dict.keys()) == expected_keys

    exported_yaml = benchmark_results.to_submission_yaml()
    parsed_yaml = yaml.safe_load(exported_yaml)
    assert set(parsed_yaml.keys()) == expected_keys


def test_to_submission_yaml_coerces_non_primitive_metadata_values() -> None:
    """Token-like objects in metadata should be converted to strings for YAML export."""

    class _FakeKey:
        def __str__(self) -> str:
            return "AlchemicalNetwork-32c2e28da705b7c006ca4f1c71bf6cf9"

    benchmark_results = _make_minimal_benchmark_results()
    benchmark_results.benchmark_data = {
        "source_repository": "https://example.invalid/repo",
        "jacs_set": {"tyk2": _FakeKey()},
    }

    exported_yaml = benchmark_results.to_submission_yaml()
    parsed_yaml = yaml.safe_load(exported_yaml)

    assert parsed_yaml["benchmark_data"]["jacs_set"]["tyk2"] == str(_FakeKey())


def test_to_submission_yaml_emits_literal_block_for_literalstr_fields() -> None:
    """LiteralStr values should be emitted as YAML block scalars."""
    benchmark_results = _make_minimal_benchmark_results()
    benchmark_results.protocol_settings = [
        {
            "protocol": "RelativeHybridTopologyProtocol",
            "notes": LiteralStr("Applies to 2 edges:\n- edge-a\n- edge-b"),
        }
    ]

    exported_yaml = benchmark_results.to_submission_yaml()

    assert "notes: |-" in exported_yaml or "notes: |" in exported_yaml
    assert "- edge-a" in exported_yaml
    assert "- edge-b" in exported_yaml
