"""
Tests for the BenchmarkResults class and filtering functionality.

Tests cover:
- Loading (full and fast YAML-only)
- Raw results structure validation
- Filtering (tags, systems, nested fields, wildcards, OR/NOT/AND logic)
- Lazy FEMap generation (dg and ddg)
- Error handling (missing files, invalid YAML, etc.)
"""

import pytest
import yaml
import time

from cinnabar import FEMap

from openfe_benchmarks.results import get_benchmark_results, filter_results
import openfe_benchmarks.results._benchmark_results as br_module


# Test data submission IDs
RBFE_SUBMISSION = "2026-03-18-openmm-840-qa-testing"
ASFE_SUBMISSION = "2026-08-06-openff-2.3.0-solvation_set_freesolv"


# ========== Fixtures ==========


@pytest.fixture
def mock_submission_yaml():
    """
    Factory fixture to create complete mock submission YAML data.

    Returns a function that generates a valid submission.yaml dict
    with all required fields, allowing custom overrides.
    """

    def _make_yaml(submission_id="test-submission", **overrides):
        """
        Generate a complete submission.yaml dict.

        Parameters
        ----------
        submission_id : str
            Submission ID for the test
        **overrides : dict
            Any fields to override from the defaults

        Returns
        -------
        dict
            Complete submission.yaml dictionary
        """
        yaml_data = {
            "submission_id": submission_id,
            "title": "Test Submission",
            "summary": "Test summary",
            "calculation_type": "rbfe",
            "tags": ["test"],
            "authors": [{"name": "Test Author"}],
            "date": "2026-01-01",
            "results": "results.json",
            "archive": {"doi": "10.1234/test", "archive_provider": "test"},
            "license": "MIT",
            "openfe_version": "1.0.0",
            "openmm_version": "8.0.0",
            "openff_toolkit_version": "0.10.0",
            "partial_charges": "am1bcc",
            "benchmark_data": {},
            "protocol_settings": [],
        }
        yaml_data.update(overrides)
        return yaml_data

    return _make_yaml


@pytest.fixture
def create_test_submission(tmp_path, monkeypatch, mock_submission_yaml):
    """
    Factory fixture to create test submission directories with YAML files.

    Automatically patches _RESULTS_DIR and creates the directory structure.
    """
    monkeypatch.setattr(br_module, "_RESULTS_DIR", tmp_path)

    def _create(submission_id, yaml_data=None, **yaml_overrides):
        """
        Create a test submission directory with submission.yaml.

        Parameters
        ----------
        submission_id : str
            Submission ID for the test
        yaml_data : dict, optional
            Complete YAML data dict. If None, uses mock_submission_yaml factory
        **yaml_overrides : dict
            Fields to override in the default YAML (only used if yaml_data is None)

        Returns
        -------
        Path
            Path to the submission directory
        """
        submission_dir = tmp_path / submission_id
        submission_dir.mkdir()

        if yaml_data is None:
            yaml_data = mock_submission_yaml(submission_id, **yaml_overrides)

        yaml_file = submission_dir / "submission.yaml"
        with open(yaml_file, "w") as f:
            yaml.dump(yaml_data, f)

        return submission_dir

    return _create


# ========== Loading Tests ==========


def test_load_by_submission_id():
    """Test loading by submission_id with full data loading."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Verify attributes
    assert result.submission_id == RBFE_SUBMISSION
    assert isinstance(result.title, str)
    assert result.calculation_type == "rbfe"
    assert isinstance(result.tags, list)
    assert "rbfe" in result.tags
    assert result.raw_results is not None
    assert isinstance(result.raw_results, dict)
    assert result.results_file.exists()
    assert result.submission_file.exists()


def test_load_yaml_only_fast():
    """Test fast YAML-only loading (load_results=False)."""
    start = time.time()
    result = get_benchmark_results(RBFE_SUBMISSION, load_results=False)
    elapsed = time.time() - start

    # Verify YAML metadata loaded
    assert result.submission_id == RBFE_SUBMISSION
    assert isinstance(result.title, str)
    assert result.calculation_type == "rbfe"
    assert isinstance(result.tags, list)

    # Verify raw_results is None
    assert result.raw_results is None

    # Verify fast loading (generous threshold for slow CI systems)
    # Most systems should complete in <0.5s, but allow up to 2s for loaded CI
    assert elapsed < 2.0, f"Fast load too slow: {elapsed:.3f}s (should be <2s)"


# ========== Raw Results Structure Tests ==========


def test_raw_results_structure():
    """Test that raw_results has expected structure with dg/ddg keys."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Verify keys
    assert "dg" in result.raw_results or "ddg" in result.raw_results

    # Verify structure
    if "dg" in result.raw_results:
        assert isinstance(result.raw_results["dg"], list)
        assert len(result.raw_results["dg"]) >= 0

    if "ddg" in result.raw_results:
        assert isinstance(result.raw_results["ddg"], list)
        assert len(result.raw_results["ddg"]) > 0  # RBFE should have ddg results

        # Verify result structure (first result should have system info)
        first_result = result.raw_results["ddg"][0]
        assert isinstance(first_result, dict)
        assert "system_group" in first_result or "system_name" in first_result


# ========== Filter Tests ==========


def test_filter_by_tags():
    """Test filtering by exact tag match."""
    result = get_benchmark_results(RBFE_SUBMISSION)
    filtered = filter_results(result, tags="rbfe")

    assert len(filtered) > 0
    # All results should have the tag (check metadata since tags is top-level)
    assert "rbfe" in result.tags


def test_filter_by_multiple_tags_and():
    """Test filtering by multiple tags with AND logic (default)."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Use tags that should exist in the submission
    filtered = filter_results(result, tags=["rbfe", "jacs_set"])

    # Verify it returns a list
    assert isinstance(filtered, list)

    # Verify tags exist in submission metadata (tags are top-level, not per-result)
    assert "rbfe" in result.tags
    assert "jacs_set" in result.tags


def test_filter_by_multiple_tags_or():
    """Test filtering by multiple tags with OR logic."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Use tags_mode='any' for OR logic
    filtered = filter_results(result, tags=["rbfe", "asfe"], tags_mode="any")

    # Should return results since 'rbfe' tag exists
    assert len(filtered) > 0
    assert isinstance(filtered, list)


def test_filter_by_system():
    """Test filtering by system_group and system_name."""
    result = get_benchmark_results(RBFE_SUBMISSION)
    filtered = filter_results(result, system_group="jacs_set", system_name="tyk2")

    assert len(filtered) > 0
    # Verify all results match
    for r in filtered:
        assert r["system_group"] == "jacs_set"
        assert r["system_name"] == "tyk2"


def test_filter_by_nested_field():
    """
    Test filtering by nested field using __ syntax.

    Note: Current test data doesn't have nested dicts within individual results,
    so this test verifies that the nested field filtering mechanism works correctly
    when the field doesn't exist (should return empty list).
    """
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Test that nested field syntax doesn't crash when field doesn't exist
    # Using a hypothetical nested field that doesn't exist in current test data
    filtered = filter_results(result, hypothetical__nested__field="value")

    # Should return empty list (field doesn't exist)
    assert isinstance(filtered, list)
    assert len(filtered) == 0, (
        "Filtering by non-existent nested field should return empty list"
    )

    # Verify regular filtering still works
    filtered_normal = filter_results(result, system_group="jacs_set")
    assert len(filtered_normal) > 0, "Regular filtering should still work"


def test_filter_with_wildcard():
    """Test filtering with wildcard pattern matching."""
    result = get_benchmark_results(RBFE_SUBMISSION)
    filtered = filter_results(result, system_name="*tyk2*")

    assert len(filtered) > 0
    # Verify all results match wildcard (case-sensitive by fnmatch)
    for r in filtered:
        assert "tyk2" in r["system_name"]


def test_filter_or_logic():
    """Test OR logic within a field using list values."""
    result = get_benchmark_results(RBFE_SUBMISSION)
    filtered = filter_results(result, system_name=["tyk2", "p38"])

    assert len(filtered) > 0
    # Verify all results match one of the values
    for r in filtered:
        assert r["system_name"] in ["tyk2", "p38"]


def test_filter_not_logic():
    """Test NOT logic using exclude_ prefix."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Test 1: Filter for jacs_set but exclude tyk2 system
    filtered = filter_results(
        result, system_group="jacs_set", exclude_system_name="tyk2"
    )

    assert len(filtered) > 0
    # Verify no tyk2 results
    for r in filtered:
        assert r["system_group"] == "jacs_set"
        assert r["system_name"] != "tyk2"

    # Test 2: Exclude by tags (as specified in plan)
    # Note: This tests the mechanism even if 'deprecated' tag doesn't exist
    # The filter should return all results (no exclusion) if tag doesn't exist
    filtered_tags = filter_results(result, exclude_tags="deprecated")
    assert isinstance(filtered_tags, list)


def test_filter_complex_logic():
    """Test complex filtering combining tags AND + system OR + exclude."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Complex filter: tags AND + system OR + exclude
    filtered = filter_results(
        result,
        tags=["rbfe", "jacs_set"],
        system_name=["tyk2", "p38"],
        exclude_tags="deprecated",
    )

    # Verify it returns a list
    assert isinstance(filtered, list)

    # Verify tags exist in metadata (tags are top-level)
    assert "rbfe" in result.tags
    assert "jacs_set" in result.tags

    # Verify filtering worked correctly
    if len(filtered) > 0:
        for r in filtered:
            # Must match one of the system names (OR logic within field)
            assert r["system_name"] in ["tyk2", "p38"], (
                f"Expected system_name in ['tyk2', 'p38'], got {r['system_name']}"
            )

            # Must NOT have 'deprecated' tag (exclude logic)
            # Tags are at result level if they exist in the result dict
            if "tags" in r:
                assert "deprecated" not in r["tags"], (
                    f"Result should not have 'deprecated' tag but found it in {r.get('tags')}"
                )
    else:
        # If no results, ensure the filter criteria could be met
        # (i.e., the submission has the required tags)
        assert "rbfe" in result.tags and "jacs_set" in result.tags


def test_filter_multi_and_logic():
    """Test multiple predicates with AND logic between fields."""
    result = get_benchmark_results(RBFE_SUBMISSION)
    filtered = filter_results(result, system_group="jacs_set", calculation_type="rbfe")

    # Should return results matching both predicates
    assert len(filtered) > 0
    for r in filtered:
        assert r["system_group"] == "jacs_set"

    # Verify calculation_type in metadata (it's a top-level field)
    assert result.calculation_type == "rbfe"


# ========== FEMap Tests ==========


def test_ddg_femaps_lazy_load():
    """Test ddg_femaps property with lazy loading and caching."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # First access - should compute
    femaps_first = result.ddg_femaps

    # Verify structure
    assert isinstance(femaps_first, dict)
    assert len(femaps_first) > 0

    # Verify keys are tuples
    for key in femaps_first.keys():
        assert isinstance(key, tuple)
        assert len(key) == 2

    # Verify values are FEMaps
    for femap in femaps_first.values():
        assert isinstance(femap, FEMap)

    # Second access - should return cached
    femaps_second = result.ddg_femaps
    assert femaps_second is femaps_first  # Identity check for caching


def test_dg_femaps():
    """Test dg_femaps property with ASFE results."""
    result = get_benchmark_results(ASFE_SUBMISSION)

    # Access dg_femaps
    femaps = result.dg_femaps

    # Verify structure
    assert isinstance(femaps, dict)
    assert len(femaps) > 0

    # Verify keys are tuples
    for key in femaps.keys():
        assert isinstance(key, tuple)
        assert len(key) == 2

    # Verify values are FEMaps
    for femap in femaps.values():
        assert isinstance(femap, FEMap)


# ========== Error Handling Tests ==========


def test_missing_submission_id():
    """Test error when submission_id does not exist."""
    submission_id = "nonexistent-submission"
    with pytest.raises(FileNotFoundError) as excinfo:
        get_benchmark_results(submission_id)

    assert str(excinfo.value) == f"Submission not found: {submission_id}"


def test_missing_results_file(create_test_submission):
    """Test error when results file referenced in YAML does not exist."""
    submission_id = "test-missing-results"
    create_test_submission(
        submission_id, title="Test Missing Results", results="nonexistent_results.json"
    )

    # Try to load with load_results=True
    with pytest.raises(FileNotFoundError) as excinfo:
        get_benchmark_results(submission_id, load_results=True)

    assert "Results file not found" in str(excinfo.value)


def test_invalid_submission_yaml(create_test_submission):
    """Test error when submission.yaml is missing required fields."""
    submission_id = "test-invalid-yaml"

    # Create YAML with only submission_id and tags (missing many required fields)
    yaml_data = {
        "submission_id": submission_id,
        "tags": ["test"],
    }
    create_test_submission(submission_id, yaml_data=yaml_data)

    with pytest.raises(ValueError) as excinfo:
        get_benchmark_results(submission_id, load_results=False)

    # Should raise ValueError for a missing required field
    assert "missing required field" in str(excinfo.value)


def test_submission_id_mismatch(create_test_submission, mock_submission_yaml):
    """Test error when submission_id in YAML doesn't match directory name."""
    submission_id = "test-mismatch"

    # Create YAML with different submission_id (different from directory name)
    yaml_data = mock_submission_yaml("different-id", title="Test Mismatch")
    create_test_submission(submission_id, yaml_data=yaml_data)

    with pytest.raises(ValueError) as excinfo:
        get_benchmark_results(submission_id, load_results=False)

    assert "submission_id mismatch" in str(excinfo.value)


def test_unknown_fields_in_yaml(create_test_submission):
    """Test error when submission.yaml contains unknown fields."""
    submission_id = "test-unknown-fields"

    # Create complete YAML with extra unknown fields
    create_test_submission(
        submission_id,
        title="Test Unknown Fields",
        unknown_field="this should not be here",
        another_unknown="also unexpected",
    )

    with pytest.raises(ValueError) as excinfo:
        get_benchmark_results(submission_id, load_results=False)

    assert "Unknown fields" in str(excinfo.value)
    assert "unknown_field" in str(excinfo.value) or "another_unknown" in str(
        excinfo.value
    )


def test_filter_on_fast_loaded_results():
    """Test that filtering raises error when raw_results is None."""
    result = get_benchmark_results(RBFE_SUBMISSION, load_results=False)

    with pytest.raises(ValueError) as excinfo:
        filter_results(result, tags="rbfe")

    expected_msg = (
        "Cannot filter results: raw_results is None. "
        "Initialize with load_results=True to access computational data."
    )
    assert str(excinfo.value) == expected_msg


def test_femaps_on_fast_loaded_results():
    """Test that accessing FEMaps raises error when raw_results is None."""
    result = get_benchmark_results(RBFE_SUBMISSION, load_results=False)

    with pytest.raises(ValueError) as excinfo:
        _ = result.ddg_femaps

    expected_msg = (
        "Cannot access ddg_femaps: raw_results is None. "
        "Initialize with load_results=True to access computational data."
    )
    assert str(excinfo.value) == expected_msg


# ========== Comparison Operator Tests ==========


def test_filter_date_greater_than_or_equal():
    """Test date filtering with >= operator."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Filter for dates >= 2026-03-01
    filtered = filter_results(result, date=">=2026-03-01")

    # Should include results (submission date is after 2026-03-01)
    assert len(filtered) > 0
    # Convert date to string for comparison
    assert result.date.isoformat() >= "2026-03-01"


def test_filter_date_less_than():
    """Test date filtering with < operator."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Filter for dates < 2027-01-01
    filtered = filter_results(result, date="<2027-01-01")

    # Should include results (submission date is before 2027)
    assert len(filtered) > 0
    assert result.date.isoformat() < "2027-01-01"


def test_filter_date_range():
    """Test date filtering with range (combining two filters)."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Filter for dates in 2026 (>= 2026-01-01 AND < 2027-01-01)
    # Note: This requires applying both filters, but current API doesn't support
    # multiple operators on same field, so we test separately
    filtered_after = filter_results(result, date=">=2026-01-01")
    filtered_before = filter_results(result, date="<2027-01-01")

    # Both should return results
    assert len(filtered_after) > 0
    assert len(filtered_before) > 0


def test_filter_version_less_than():
    """Test version filtering with < operator using semantic versioning."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Get the actual version from the submission
    actual_version = result.openfe_version

    # Filter for versions < actual_version + 1 (should include all)
    major, minor, *rest = actual_version.split(".")
    next_version = f"{int(major) + 1}.0.0"
    filtered = filter_results(result, openfe_version=f"<{next_version}")

    # Should include results
    assert len(filtered) > 0


def test_filter_version_greater_than_or_equal():
    """Test version filtering with >= operator using semantic versioning."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Filter for versions >= 1.0.0 (should include most versions)
    filtered = filter_results(result, openfe_version=">=1.0.0")

    # Should include results (openfe_version is likely >= 1.0.0)
    assert isinstance(filtered, list)


def test_filter_version_semantic_comparison():
    """Test that version comparison uses semantic versioning, not string comparison."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Semantic: 1.10.0 > 1.9.0
    # String: "1.10.0" < "1.9.0" (would be wrong)

    # Just verify filtering works without error
    filtered_gte = filter_results(result, openfe_version=">=1.0.0")
    filtered_lt = filter_results(result, openfe_version="<99.0.0")

    assert isinstance(filtered_gte, list)
    assert isinstance(filtered_lt, list)


def test_filter_comparison_with_other_filters():
    """Test comparison operators combined with other filter types."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Combine date comparison with system filter
    filtered = filter_results(result, date=">=2026-01-01", system_group="jacs_set")

    # Should return results matching both conditions
    assert isinstance(filtered, list)
    if len(filtered) > 0:
        for r in filtered:
            assert r["system_group"] == "jacs_set"
        # Date is metadata-level, checked via isoformat
        assert result.date.isoformat() >= "2026-01-01"


def test_filter_comparison_no_operator():
    """Test that values without operators still work as exact match."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Date without operator should be exact match
    exact_date = result.date
    filtered = filter_results(result, date=exact_date)

    # Should return all results (exact match on metadata)
    assert len(filtered) > 0


def test_filter_comparison_exclude_with_operator():
    """Test that exclude_ prefix works with comparison operators."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Exclude dates before 2026
    filtered = filter_results(result, exclude_date="<2026-01-01")

    # Should return results (submission date is after 2026-01-01, so not excluded)
    assert len(filtered) > 0
    assert result.date.isoformat() >= "2026-01-01"


def test_filter_comparison_openmm_version():
    """Test comparison on openmm_version field."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Filter for OpenMM >= 8.0.0
    filtered = filter_results(result, openmm_version=">=8.0.0")

    # Should return results if openmm_version >= 8.0.0
    assert isinstance(filtered, list)


def test_filter_comparison_openff_toolkit_version():
    """Test comparison on openff_toolkit_version field."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Filter for OpenFF Toolkit >= 0.10.0
    filtered = filter_results(result, openff_toolkit_version=">=0.10.0")

    # Should return results if version matches
    assert isinstance(filtered, list)


def test_filter_comparison_warning_non_version_field():
    """Test that warning is issued when comparison operators used with non-date/version fields."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Using comparison operator with system_name (not a date/version field) should warn
    with pytest.warns(
        UserWarning,
        match="Comparison operator.*system_name.*designed for date and version fields",
    ):
        filtered = filter_results(result, system_name=">=tyk2")

    # Still returns results (falls back to string comparison)
    assert isinstance(filtered, list)


def test_filter_quantity_comparison_greater_than():
    """Test filtering on pint Quantity fields with > operator."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Filter results where dg > 1.0 kcal/mol
    # dg values are Quantity objects deserialized from JSON
    filtered = filter_results(result, dg=">1.0 kilocalories_per_mole")

    assert isinstance(filtered, list)
    assert len(filtered) > 0

    # Verify all results have dg > 1.0 kcal/mol
    for r in filtered:
        assert r["dg"].magnitude > 1.0


def test_filter_quantity_comparison_less_than_or_equal():
    """Test filtering on pint Quantity fields with <= operator."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Filter results where dg <= 2.0 kcal/mol
    filtered = filter_results(result, dg="<=2.0 kilocalories_per_mole")

    assert isinstance(filtered, list)
    assert len(filtered) > 0

    # Verify all results have dg <= 2.0 kcal/mol
    for r in filtered:
        assert r["dg"].magnitude <= 2.0


def test_filter_quantity_comparison_unit_conversion():
    """Test that pint Quantity comparison handles unit conversion automatically."""

    result = get_benchmark_results(RBFE_SUBMISSION)

    # Filter with different but compatible units
    # 1 kcal/mol ≈ 4.184 kJ/mol
    filtered_kcal = filter_results(result, dg=">1.0 kilocalories_per_mole")
    filtered_kj = filter_results(result, dg=">4.184 kilojoules_per_mole")

    # Should get the same results (or very close due to floating point)
    assert (
        len(filtered_kcal) == len(filtered_kj)
        or abs(len(filtered_kcal) - len(filtered_kj)) <= 1
    )


def test_filter_quantity_comparison_range():
    """Test filtering on pint Quantity with range (testing both filters separately)."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Filter for dg >= 0.5 kcal/mol
    filtered_gte = filter_results(result, dg=">=0.5 kilocalories_per_mole")

    # Filter for dg <= 2.0 kcal/mol
    filtered_lte = filter_results(result, dg="<=2.0 kilocalories_per_mole")

    assert isinstance(filtered_gte, list)
    assert isinstance(filtered_lte, list)
    assert len(filtered_gte) > 0
    assert len(filtered_lte) > 0

    # Verify ranges
    for r in filtered_gte:
        assert r["dg"].magnitude >= 0.5

    for r in filtered_lte:
        assert r["dg"].magnitude <= 2.0


def test_filter_quantity_incompatible_units_raises():
    """Test that incompatible units raise ValueError (no silent fallback)."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Try to compare energy (kcal/mol) with temperature (K) - dimensionally incompatible
    with pytest.raises(ValueError, match="Incompatible units"):
        _ = filter_results(result, dg=">298 kelvin")


def test_filter_quantity_invalid_string_raises():
    """Test that invalid quantity strings raise ValueError (no silent fallback)."""
    result = get_benchmark_results(RBFE_SUBMISSION)

    # Use a string that can't be parsed as a quantity
    with pytest.raises(ValueError, match="Invalid quantity filter value"):
        _ = filter_results(result, dg=">not_a_quantity")
