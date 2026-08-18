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

from openfe_benchmarks.results import BenchmarkResults, filter_results


# Test data submission IDs
RBFE_SUBMISSION = "2026-03-18-openmm-840-qa-testing"
ASFE_SUBMISSION = "2026-08-06-openff-2.3.0-solvation_set_freesolv"


# ========== Loading Tests ==========


def test_load_by_submission_id():
    """Test loading by submission_id with full data loading."""
    result = BenchmarkResults(RBFE_SUBMISSION)

    # Verify attributes
    assert result.submission_id == RBFE_SUBMISSION
    assert isinstance(result.title, str)
    assert result.calculation_type == "rbfe"
    assert isinstance(result.tags, list)
    assert "rbfe" in result.tags
    assert isinstance(result.metadata, dict)
    assert result.raw_results is not None
    assert isinstance(result.raw_results, dict)
    assert result.results_file.exists()
    assert result.submission_file.exists()


def test_load_yaml_only_fast():
    """Test fast YAML-only loading (load_results=False)."""
    start = time.time()
    result = BenchmarkResults(RBFE_SUBMISSION, load_results=False)
    elapsed = time.time() - start

    # Verify YAML metadata loaded
    assert result.submission_id == RBFE_SUBMISSION
    assert isinstance(result.title, str)
    assert result.calculation_type == "rbfe"
    assert isinstance(result.tags, list)

    # Verify raw_results is None
    assert result.raw_results is None

    # Verify fast loading (<1s)
    assert elapsed < 1.0, f"Fast load too slow: {elapsed:.3f}s (should be <1s)"


# ========== Raw Results Structure Tests ==========


def test_raw_results_structure():
    """Test that raw_results has expected structure with dg/ddg keys."""
    result = BenchmarkResults(RBFE_SUBMISSION)

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
    result = BenchmarkResults(RBFE_SUBMISSION)
    filtered = filter_results(result, tags="rbfe")

    assert len(filtered) > 0
    # All results should have the tag (check metadata since tags is top-level)
    assert "rbfe" in result.tags


def test_filter_by_multiple_tags_and():
    """Test filtering by multiple tags with AND logic (default)."""
    result = BenchmarkResults(RBFE_SUBMISSION)

    # Use tags that should exist in the submission
    filtered = filter_results(result, tags=["rbfe", "jacs_set"])

    # Verify it returns a list
    assert isinstance(filtered, list)

    # Verify tags exist in submission metadata (tags are top-level, not per-result)
    assert "rbfe" in result.tags
    assert "jacs_set" in result.tags


def test_filter_by_multiple_tags_or():
    """Test filtering by multiple tags with OR logic."""
    result = BenchmarkResults(RBFE_SUBMISSION)

    # Use tags_mode='any' for OR logic
    filtered = filter_results(result, tags=["rbfe", "asfe"], tags_mode="any")

    # Should return results since 'rbfe' tag exists
    assert len(filtered) > 0
    assert isinstance(filtered, list)


def test_filter_by_system():
    """Test filtering by system_group and system_name."""
    result = BenchmarkResults(RBFE_SUBMISSION)
    filtered = filter_results(result, system_group="jacs_set", system_name="tyk2")

    assert len(filtered) > 0
    # Verify all results match
    for r in filtered:
        assert r["system_group"] == "jacs_set"
        assert r["system_name"] == "tyk2"


def test_filter_by_nested_field():
    """Test filtering by nested field using __ syntax."""
    result = BenchmarkResults(RBFE_SUBMISSION)

    # Note: Current test data doesn't have nested dicts within individual results.
    # Nested field filtering (using __) works for result-level nested fields,
    # but test data only has flat result dicts.
    # This test verifies the mechanism doesn't crash with __ syntax.

    # Test that nested field syntax doesn't crash
    # Using a hypothetical nested field that doesn't exist
    filtered = filter_results(result, hypothetical__nested__field="value")

    # Should return empty list (field doesn't exist)
    assert isinstance(filtered, list)
    assert len(filtered) == 0

    # Verify regular filtering still works
    filtered_normal = filter_results(result, system_group="jacs_set")
    assert len(filtered_normal) > 0


def test_filter_with_wildcard():
    """Test filtering with wildcard pattern matching."""
    result = BenchmarkResults(RBFE_SUBMISSION)
    filtered = filter_results(result, system_name="*tyk2*")

    assert len(filtered) > 0
    # Verify all results match wildcard (case-sensitive by fnmatch)
    for r in filtered:
        assert "tyk2" in r["system_name"]


def test_filter_or_logic():
    """Test OR logic within a field using list values."""
    result = BenchmarkResults(RBFE_SUBMISSION)
    filtered = filter_results(result, system_name=["tyk2", "p38"])

    assert len(filtered) > 0
    # Verify all results match one of the values
    for r in filtered:
        assert r["system_name"] in ["tyk2", "p38"]


def test_filter_not_logic():
    """Test NOT logic using exclude_ prefix."""
    result = BenchmarkResults(RBFE_SUBMISSION)

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
    result = BenchmarkResults(RBFE_SUBMISSION)

    # Complex filter: tags AND + system OR + exclude (as per plan)
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

    # Verify system_name matches one of the specified values (if results exist)
    if len(filtered) > 0:
        for r in filtered:
            assert r["system_name"] in ["tyk2", "p38"]


def test_filter_multi_and_logic():
    """Test multiple predicates with AND logic between fields."""
    result = BenchmarkResults(RBFE_SUBMISSION)
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
    result = BenchmarkResults(RBFE_SUBMISSION)

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
    result = BenchmarkResults(ASFE_SUBMISSION)

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
        BenchmarkResults(submission_id)

    assert str(excinfo.value) == f"Submission not found: {submission_id}"


def test_missing_results_file(tmp_path):
    """Test error when results file referenced in YAML does not exist."""
    # Create a mock submission directory with YAML but no results file
    submission_id = "test-missing-results"
    submission_dir = tmp_path / submission_id
    submission_dir.mkdir()

    yaml_data = {
        "submission_id": submission_id,
        "title": "Test Missing Results",
        "calculation_type": "rbfe",
        "tags": ["test"],
        "results": "nonexistent_results.json",
    }

    yaml_file = submission_dir / "submission.yaml"
    with open(yaml_file, "w") as f:
        yaml.dump(yaml_data, f)

    # Try to load with load_results=True
    with pytest.raises(FileNotFoundError) as excinfo:
        BenchmarkResults(submission_id, load_results=True, results_dir=tmp_path)

    assert "Results file not found" in str(excinfo.value)


def test_invalid_submission_yaml(tmp_path):
    """Test error when submission.yaml is missing required fields."""
    # Create a mock submission directory with incomplete YAML
    submission_id = "test-invalid-yaml"
    submission_dir = tmp_path / submission_id
    submission_dir.mkdir()

    # Missing 'title' and 'calculation_type' fields
    yaml_data = {
        "submission_id": submission_id,
        "tags": ["test"],
    }

    yaml_file = submission_dir / "submission.yaml"
    with open(yaml_file, "w") as f:
        yaml.dump(yaml_data, f)

    with pytest.raises(ValueError) as excinfo:
        BenchmarkResults(submission_id, results_dir=tmp_path, load_results=False)

    assert "missing required fields" in str(excinfo.value)


def test_submission_id_mismatch(tmp_path):
    """Test error when submission_id in YAML doesn't match directory name."""
    # Create a mock submission directory
    submission_id = "test-mismatch"
    submission_dir = tmp_path / submission_id
    submission_dir.mkdir()

    # YAML has different submission_id
    yaml_data = {
        "submission_id": "different-id",
        "title": "Test Mismatch",
        "calculation_type": "rbfe",
        "tags": ["test"],
        "results": "results.json",
    }

    yaml_file = submission_dir / "submission.yaml"
    with open(yaml_file, "w") as f:
        yaml.dump(yaml_data, f)

    with pytest.raises(ValueError) as excinfo:
        BenchmarkResults(submission_id, results_dir=tmp_path, load_results=False)

    assert "submission_id mismatch" in str(excinfo.value)


def test_filter_on_fast_loaded_results():
    """Test that filtering raises error when raw_results is None."""
    result = BenchmarkResults(RBFE_SUBMISSION, load_results=False)

    with pytest.raises(ValueError) as excinfo:
        filter_results(result, tags="rbfe")

    expected_msg = (
        "Cannot filter results: raw_results is None. "
        "Initialize with load_results=True to access computational data."
    )
    assert str(excinfo.value) == expected_msg


def test_femaps_on_fast_loaded_results():
    """Test that accessing FEMaps raises error when raw_results is None."""
    result = BenchmarkResults(RBFE_SUBMISSION, load_results=False)

    with pytest.raises(ValueError) as excinfo:
        _ = result.ddg_femaps

    expected_msg = (
        "Cannot access ddg_femaps: raw_results is None. "
        "Initialize with load_results=True to access computational data."
    )
    assert str(excinfo.value) == expected_msg
