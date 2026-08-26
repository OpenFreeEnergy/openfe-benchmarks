"""
Tests for CI-optimized validation helpers.

Validates fast YAML-only validation, git-aware changed file detection,
random sampling with calculation_type coverage, and hybrid CI workflow timing.

Target: <6 min CI time for 100 submissions with 10% changed.
"""

import time
from unittest.mock import patch
import yaml

from openfe_benchmarks.results._validation import (
    validate_submission_yaml_fast,
    get_all_submission_ids,
    get_changed_submission_ids,
    select_random_sample,
)
from openfe_benchmarks.results import get_benchmark_results
from openfe_benchmarks.results._benchmark_results import _RESULTS_DIR


def _trigger_available_femaps(benchmark_result):
    """Load only FEMap types that exist in raw_results for this submission."""
    loaded_any = False

    if "ddg" in benchmark_result.raw_results:
        _ = benchmark_result.ddg_femaps()
        loaded_any = True

    if "dg" in benchmark_result.raw_results:
        _ = benchmark_result.dg_femaps()
        loaded_any = True

    assert loaded_any, (
        f"Expected at least one of 'dg' or 'ddg' in raw_results for "
        f"{benchmark_result.submission_id}, found keys: "
        f"{list(benchmark_result.raw_results.keys())}"
    )


def test_fast_yaml_validation_rejects_noncanonical_path(tmp_path):
    """Test that validate_submission_yaml_fast rejects non-canonical YAML paths."""
    submission_id = "standalone-submission"
    submission_dir = tmp_path / submission_id
    submission_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = submission_dir / "submission.yaml"

    yaml_data = {
        "submission_id": submission_id,
        "title": "Standalone Test Submission",
        "summary": "Validate explicit path handling",
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
    with open(yaml_path, "w") as f:
        yaml.dump(yaml_data, f)

    result = validate_submission_yaml_fast(yaml_path)

    assert result["valid"] is False
    assert result["submission_id"] is None
    assert result["calculation_type"] is None
    assert any("Non-canonical submission path" in err for err in result["errors"])


def test_fast_yaml_validation_all():
    """
    Test fast YAML-only validation for all submissions.

    Success criteria:
    - All submissions validate without errors
    - Validation completes within reasonable time per submission
    """
    # Get all real submission_ids from results directory
    all_ids = get_all_submission_ids()

    assert len(all_ids) > 0, "Expected at least one submission in results directory"

    print(f"\n{'=' * 60}")
    print(f"Fast YAML Validation: {len(all_ids)} submissions")
    print(f"{'=' * 60}")

    invalid_submissions = []

    for submission_id in all_ids:
        yaml_path = _RESULTS_DIR / submission_id / "submission.yaml"
        result = validate_submission_yaml_fast(yaml_path)

        status = "✓" if result["valid"] else "✗"
        print(f"  {status} {submission_id}")

        if not result["valid"]:
            invalid_submissions.append((submission_id, result["errors"]))
            print(f"    Errors: {result['errors']}")

    # Assert all submissions are valid
    if invalid_submissions:
        error_msg = f"\n{len(invalid_submissions)} invalid submission(s) found:\n"
        for submission_id, errors in invalid_submissions:
            error_msg += f"\n  {submission_id}:\n"
            for error in errors:
                error_msg += f"    - {error}\n"
        assert False, error_msg


def test_fast_yaml_validation_performance():
    """
    Test performance of YAML validation.

    Success criteria:
    - Individual submissions validate in reasonable time (<5s each)
    - Average time per submission is tracked for regression detection
    """
    all_ids = get_all_submission_ids()

    assert len(all_ids) > 0, "Expected at least one submission"

    print(f"\n{'=' * 60}")
    print(f"YAML Validation Performance: {len(all_ids)} submissions")
    print(f"{'=' * 60}")

    total_start = time.time()
    timings = []

    for submission_id in all_ids:
        yaml_path = _RESULTS_DIR / submission_id / "submission.yaml"

        start = time.time()
        _ = validate_submission_yaml_fast(yaml_path)
        elapsed = time.time() - start

        timings.append(elapsed)
        print(f"  {submission_id}: {elapsed:.3f}s")

        # Warn if validation is slow (>1s) but only fail if extremely slow (>5s)
        if elapsed > 1.0:
            print(f"    WARNING: Slower than target 1.0s (actual: {elapsed:.3f}s)")

        assert elapsed < 5.0, (
            f"YAML validation extremely slow for {submission_id}: {elapsed:.3f}s "
            f"(failing threshold: 5s, target: 1s)"
        )

    total_elapsed = time.time() - total_start
    avg_time = sum(timings) / len(timings)
    max_time = max(timings)

    print(f"\n{'=' * 60}")
    print("Performance Summary:")
    print(f"  Total time: {total_elapsed:.3f}s")
    print(f"  Average per submission: {avg_time:.3f}s")
    print(f"  Max per submission: {max_time:.3f}s")
    print(f"{'=' * 60}\n")


def test_fast_yaml_validation_scaling():
    """
    Test that YAML validation scales acceptably for large submission sets.

    Success criteria:
    - Extrapolated time for 100 submissions is reasonable (<500s)
    """
    all_ids = get_all_submission_ids()

    assert len(all_ids) > 0, "Expected at least one submission"

    # Quick timing sample
    timings = []
    for submission_id in all_ids:
        yaml_path = _RESULTS_DIR / submission_id / "submission.yaml"
        start = time.time()
        _ = validate_submission_yaml_fast(yaml_path)
        elapsed = time.time() - start
        timings.append(elapsed)

    avg_time = sum(timings) / len(timings)
    extrapolated_100 = avg_time * 100

    print(f"\n{'=' * 60}")
    print("Scaling Analysis:")
    print(f"  Current submissions: {len(all_ids)}")
    print(f"  Average time: {avg_time:.3f}s")
    print(f"  Extrapolated for 100: {extrapolated_100:.1f}s")
    print("  Target: <100s, Failing threshold: <500s")
    print(f"{'=' * 60}\n")

    # Use 5x tolerance for CI variability (target 100s, fail at 500s)
    if extrapolated_100 > 100:
        print(
            f"  INFO: Extrapolated time ({extrapolated_100:.1f}s) exceeds target of 100s"
        )

    assert extrapolated_100 < 500, (
        f"Fast validation would take {extrapolated_100:.1f}s for 100 submissions "
        f"(failing threshold: 500s, target: 100s)"
    )


def test_changed_files_detection():
    """
    Test git diff-based changed submission detection.

    Success criteria:
    - Correctly extracts submission_ids from git diff output
    - Handles both submission.yaml and computational_results.json changes
    - Returns empty list gracefully when git unavailable
    """
    # Mock git diff output
    mock_git_output = """openfe_benchmarks/results/2026-03-18-openmm-840-qa-testing/submission.yaml
openfe_benchmarks/results/2026-03-18-openmm-840-qa-testing/computational_results.json.bz2
openfe_benchmarks/results/2026-08-06-openff-2.3.0-solvation_set_freesolv/submission.yaml
openfe_benchmarks/data/something_else.py
"""

    # Mock subprocess.run to return our test data
    with patch("openfe_benchmarks.results._validation.subprocess.run") as mock_run:
        mock_run.return_value.stdout = mock_git_output
        mock_run.return_value.returncode = 0

        changed_ids = get_changed_submission_ids(base_branch="origin/main")

        # Verify correct submission_ids extracted
        assert "2026-03-18-openmm-840-qa-testing" in changed_ids
        assert "2026-08-06-openff-2.3.0-solvation_set_freesolv" in changed_ids
        assert len(changed_ids) == 2

        # Verify git command called correctly
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        assert call_args == ["git", "diff", "--name-only", "origin/main...HEAD"]

    # Test graceful failure when git not available
    with patch("openfe_benchmarks.results._validation.subprocess.run") as mock_run:
        mock_run.side_effect = FileNotFoundError("git not found")

        changed_ids = get_changed_submission_ids()
        assert changed_ids == []


def test_random_sampling_reproducibility():
    """
    Test reproducibility of random sampling with seeds.

    Success criteria:
    - Same seed produces identical samples
    - Different seeds produce different samples (when sampling is truly random)
    """
    all_ids = get_all_submission_ids()

    # Filter to only valid submissions
    valid_ids = []
    for submission_id in all_ids:
        yaml_path = _RESULTS_DIR / submission_id / "submission.yaml"
        result = validate_submission_yaml_fast(yaml_path)
        if result["valid"]:
            valid_ids.append(submission_id)

    assert len(valid_ids) >= 2, (
        f"Need at least 2 valid submissions, found {len(valid_ids)}"
    )

    # Test 10% sampling reproducibility
    sample1 = select_random_sample(valid_ids, sample_rate=0.1, seed=42)
    sample2 = select_random_sample(valid_ids, sample_rate=0.1, seed=42)
    sample3 = select_random_sample(valid_ids, sample_rate=0.1, seed=99)

    print(f"\n{'=' * 60}")
    print("Random Sampling Reproducibility Test")
    print(f"{'=' * 60}")
    print(f"Total submissions: {len(all_ids)}")
    print(f"Valid submissions: {len(valid_ids)}")
    print(f"Sample 1 (seed=42): {sample1}")
    print(f"Sample 2 (seed=42): {sample2}")
    print(f"Sample 3 (seed=99): {sample3}")
    print(f"{'=' * 60}\n")

    # Verify reproducibility with same seed
    assert sample1 == sample2, "Same seed should produce identical samples"

    # Count calculation types to determine if sampling is deterministic
    calc_types = set()
    for submission_id in valid_ids:
        benchmark = get_benchmark_results(submission_id, load_results=False)
        calc_types.add(benchmark.calculation_type)

    n_types = len(calc_types)
    target_count = max(int(len(valid_ids) * 0.1), n_types)

    # Only expect different samples if target_count > n_types (true random sampling occurs)
    if target_count > n_types:
        assert sample1 != sample3 or len(sample1) == 1, (
            "Different seeds should usually produce different samples when target > n_types"
        )
    else:
        print(
            f"Note: Deterministic sampling (target={target_count} == n_types={n_types})"
        )


def test_random_sampling_coverage():
    """
    Test that random sampling provides calculation_type coverage.

    Success criteria:
    - Sample size is appropriate (at least 1 per calc_type or 10% of valid_ids)
    - All calculation_types are represented when multiple types exist
    """
    all_ids = get_all_submission_ids()

    # Filter to only valid submissions
    valid_ids = []
    for submission_id in all_ids:
        yaml_path = _RESULTS_DIR / submission_id / "submission.yaml"
        result = validate_submission_yaml_fast(yaml_path)
        if result["valid"]:
            valid_ids.append(submission_id)

    assert len(valid_ids) >= 2, (
        f"Need at least 2 valid submissions, found {len(valid_ids)}"
    )

    # Get sample
    sample = select_random_sample(valid_ids, sample_rate=0.1, seed=42)

    # Count calculation types
    calc_types = set()
    for submission_id in valid_ids:
        benchmark = get_benchmark_results(submission_id, load_results=False)
        calc_types.add(benchmark.calculation_type)

    calc_types_in_sample = set()
    for submission_id in sample:
        benchmark = get_benchmark_results(submission_id, load_results=False)
        calc_types_in_sample.add(benchmark.calculation_type)

    print(f"\n{'=' * 60}")
    print("Random Sampling Coverage Test")
    print(f"{'=' * 60}")
    print(f"Total valid submissions: {len(valid_ids)}")
    print(f"Sample size (10%): {len(sample)}")
    print(f"Calculation types in valid: {sorted(calc_types)}")
    print(f"Calculation types in sample: {sorted(calc_types_in_sample)}")
    print(f"{'=' * 60}\n")

    # Verify sample size
    n_types = len(calc_types)
    target_count = max(int(len(valid_ids) * 0.1), n_types)
    expected_min = max(1, int(target_count * 0.8))
    expected_max = int(target_count * 1.2) + 1

    assert expected_min <= len(sample) <= expected_max, (
        f"Sample size {len(sample)} not in expected range [{expected_min}, {expected_max}] "
        f"(target: {target_count}, n_types: {n_types})"
    )

    # Verify calculation_type coverage (at least one per type if multiple types exist)
    if n_types > 1:
        assert len(calc_types_in_sample) >= min(len(sample), n_types), (
            f"Expected coverage of calculation_types, got {calc_types_in_sample}"
        )


def test_deep_validation_changed_only():
    """
    Test deep validation (YAML+JSON+FEMap) for changed files only.

    Success criteria:
    - Deep validation completes in <5s per submission
    - Loads full BenchmarkResults with FEMaps
    - Timing reported for regression tracking
    """

    all_ids = get_all_submission_ids()

    # Filter to only valid submissions (YAML validation passes)
    valid_ids = []
    for submission_id in all_ids:
        yaml_path = _RESULTS_DIR / submission_id / "submission.yaml"
        result = validate_submission_yaml_fast(yaml_path)
        if result["valid"]:
            valid_ids.append(submission_id)

    # Select a few submissions for deep validation (simulate changed files)
    test_ids = valid_ids[: min(2, len(valid_ids))]

    assert len(test_ids) > 0, (
        "Need at least one valid submission for deep validation test"
    )

    print(f"\n{'=' * 60}")
    print(f"Deep Validation: {len(test_ids)} submissions")
    print(f"{'=' * 60}")

    timings = []

    for submission_id in test_ids:
        start = time.time()

        # Full load: YAML + JSON + FEMap generation
        result = get_benchmark_results(submission_id, load_results=True)

        # Access available FEMaps to trigger generation.
        _trigger_available_femaps(result)

        elapsed = time.time() - start
        timings.append(elapsed)

        print(f"  {submission_id}: {elapsed:.3f}s")

        # Warn if slow but only fail if extremely slow
        if elapsed > 5.0:
            print(f"    WARNING: Slower than target 5.0s (actual: {elapsed:.3f}s)")

        assert elapsed < 25.0, (
            f"Deep validation extremely slow for {submission_id}: {elapsed:.3f}s "
            f"(failing threshold: 25s, target: 5s)"
        )

    avg_time = sum(timings) / len(timings) if timings else 0

    print(f"\n{'=' * 60}")
    print("Summary:")
    print(f"  Average per submission: {avg_time:.3f}s")
    print(f"  Extrapolated for 10 changed: {avg_time * 10:.1f}s")
    print(f"{'=' * 60}\n")


def test_ci_hybrid_workflow():
    """
    Integration test simulating CI build with hybrid validation strategy.

    Success criteria:
    - Fast validation for all submissions
    - Deep validation for changed submissions only
    - Random 10% sample for FEMap validation
    - Total time <30s for 8 submissions (extrapolates to <200s for 100)
    - Timing breakdown reported

    Note: Uses only valid submissions for deep validation phases.
    """

    all_ids = get_all_submission_ids()

    # Use real data (typically 8 submissions as of 2026-08-18)
    assert len(all_ids) >= 1, "Need at least 1 submission for CI workflow test"

    print(f"\n{'=' * 60}")
    print("CI Hybrid Workflow Simulation")
    print(f"{'=' * 60}")
    print(f"Total submissions: {len(all_ids)}")

    workflow_start = time.time()

    # Phase 1: Fast YAML validation for ALL submissions
    print(f"\nPhase 1: Fast YAML validation (all {len(all_ids)} submissions)")
    fast_start = time.time()

    valid_ids = []
    for submission_id in all_ids:
        yaml_path = _RESULTS_DIR / submission_id / "submission.yaml"
        result = validate_submission_yaml_fast(yaml_path)
        if result["valid"]:
            valid_ids.append(submission_id)
        else:
            print(f"  Skipping invalid submission {submission_id} for deep validation")

    fast_elapsed = time.time() - fast_start
    print(f"  Completed in {fast_elapsed:.3f}s")
    print(f"  Valid submissions for deep validation: {len(valid_ids)}/{len(all_ids)}")

    # Phase 2: Deep validation for CHANGED submissions (mock 1 changed)
    # In real CI, this would use get_changed_submission_ids()
    changed_ids = valid_ids[: min(1, len(valid_ids))]

    print(
        f"\nPhase 2: Deep validation (changed files only: {len(changed_ids)} submissions)"
    )
    deep_start = time.time()

    for submission_id in changed_ids:
        result = get_benchmark_results(submission_id, load_results=True)
        _trigger_available_femaps(result)

    deep_elapsed = time.time() - deep_start
    print(f"  Completed in {deep_elapsed:.3f}s")

    # Phase 3: Random 10% sample for FEMap validation (only from valid submissions)
    sample_ids = select_random_sample(valid_ids, sample_rate=0.1, seed=42)
    # Exclude already-validated changed submissions to avoid double-counting
    sample_ids = [sid for sid in sample_ids if sid not in changed_ids]

    print(f"\nPhase 3: Random sample validation (10%: {len(sample_ids)} submissions)")
    sample_start = time.time()

    for submission_id in sample_ids:
        result = get_benchmark_results(submission_id, load_results=True)
        _trigger_available_femaps(result)

    sample_elapsed = time.time() - sample_start
    print(f"  Completed in {sample_elapsed:.3f}s")

    # Calculate total workflow time
    workflow_elapsed = time.time() - workflow_start

    print(f"\n{'=' * 60}")
    print("Timing Breakdown:")
    print(
        f"  Phase 1 (fast all):      {fast_elapsed:6.3f}s ({fast_elapsed / workflow_elapsed * 100:5.1f}%)"
    )
    print(
        f"  Phase 2 (deep changed):  {deep_elapsed:6.3f}s ({deep_elapsed / workflow_elapsed * 100:5.1f}%)"
    )
    print(
        f"  Phase 3 (random sample): {sample_elapsed:6.3f}s ({sample_elapsed / workflow_elapsed * 100:5.1f}%)"
    )
    print(f"  {'─' * 40}")
    print(f"  Total workflow time:     {workflow_elapsed:6.3f}s")
    print(f"{'=' * 60}")

    # Extrapolate to 100 submissions with 10% changed
    # Scale fast validation linearly
    fast_100 = fast_elapsed * (100 / len(all_ids))
    # Assume 10 changed (10% of 100)
    deep_100 = deep_elapsed * (10 / max(len(changed_ids), 1))
    # Assume 10 in sample (10% of 100)
    sample_100 = sample_elapsed * (10 / max(len(sample_ids), 1))
    total_100 = fast_100 + deep_100 + sample_100

    print("\nExtrapolated for 100 submissions (10% changed):")
    print(f"  Phase 1 (fast all):      {fast_100:6.1f}s")
    print(f"  Phase 2 (deep changed):  {deep_100:6.1f}s")
    print(f"  Phase 3 (random sample): {sample_100:6.1f}s")
    print(f"  {'─' * 40}")
    print(f"  Total:                   {total_100:6.1f}s ({total_100 / 60:.1f} min)")
    print("  Target:                  <360s (<6 min)")
    print(f"{'=' * 60}\n")

    # Assert total time for current test data is reasonable
    # For 8 submissions: <30s target (relaxed to allow for deep validation)
    max_time = 30.0 if len(all_ids) <= 8 else 30.0 * (len(all_ids) / 8)
    assert workflow_elapsed < max_time, (
        f"CI workflow too slow: {workflow_elapsed:.1f}s (target: <{max_time:.1f}s)"
    )

    # Assert extrapolated time meets target
    assert total_100 < 360, (
        f"Extrapolated CI time {total_100:.1f}s exceeds 6 min target"
    )
