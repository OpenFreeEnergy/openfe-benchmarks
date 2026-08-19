"""
CI-optimized validation helpers for BenchmarkResults.

This module provides fast validation functions for CI scalability:
- YAML-only validation using get_benchmark_results(load_results=False)
- Git-aware changed file detection
- Random sampling with calculation_type coverage

Target CI performance: <6 min for 100 submissions with 10% changed.
"""

import subprocess
import random
from pathlib import Path

from openfe_benchmarks.results import get_benchmark_results
from openfe_benchmarks.results._benchmark_results import _RESULTS_DIR


def validate_submission_yaml_fast(yaml_path: Path) -> dict:
    """
    Fast YAML-only validation using get_benchmark_results(load_results=False).

    Validates by attempting to construct BenchmarkResults object without
    loading JSON data or generating FEMaps. Delegates validation to the
    actual factory function that will be used.

    Parameters
    ----------
    yaml_path : Path
        Path to submission.yaml file

    Returns
    -------
    dict
        Validation result with keys:
        - valid : bool
            True if BenchmarkResults construction succeeded
        - errors : list[str]
            List of validation errors (empty if valid)
        - submission_id : str or None
            Submission ID if present in YAML
        - calculation_type : str or None
            Calculation type if present in YAML

    Examples
    --------
    >>> result = validate_submission_yaml_fast(Path("results/my-submission/submission.yaml"))
    >>> if not result['valid']:
    ...     print(f"Errors: {result['errors']}")
    """

    result = {
        "valid": True,
        "errors": [],
        "submission_id": None,
        "calculation_type": None,
    }

    # Check file exists
    if not yaml_path.exists():
        result["valid"] = False
        result["errors"].append(f"File not found: {yaml_path}")
        return result

    # Extract submission_id from path
    submission_id = yaml_path.parent.name

    # Attempt to construct BenchmarkResults without loading data
    try:
        benchmark = get_benchmark_results(
            submission_id=submission_id, load_results=False
        )

        # Extract metadata
        result["submission_id"] = benchmark.submission_id
        result["calculation_type"] = benchmark.calculation_type

        # Verify submission_id matches directory name
        if result["submission_id"] and result["submission_id"] != submission_id:
            result["valid"] = False
            result["errors"].append(
                f"submission_id mismatch: YAML has '{result['submission_id']}', "
                f"directory is '{submission_id}'"
            )

    except (FileNotFoundError, ValueError) as e:
        result["valid"] = False
        result["errors"].append(str(e))
    except Exception as e:
        result["valid"] = False
        result["errors"].append(f"Unexpected error: {e}")

    return result


def get_all_submission_ids() -> list[str]:
    """
    Find all submission_ids in results directory.

    Discovers directories containing submission.yaml files.

    Returns
    -------
    list[str]
        List of submission_ids (directory names containing submission.yaml)

    Examples
    --------
    >>> ids = get_all_submission_ids()
    >>> print(f"Found {len(ids)} submissions")
    """

    if not _RESULTS_DIR.exists():
        return []

    submission_ids = []

    # Find all directories with submission.yaml
    for item in _RESULTS_DIR.iterdir():
        if item.is_dir():
            yaml_path = item / "submission.yaml"
            if yaml_path.exists():
                submission_ids.append(item.name)

    return sorted(submission_ids)


def get_changed_submission_ids(base_branch: str = "origin/main") -> list[str]:
    """
    Get submission_ids that changed in current branch.

    Uses git diff to detect changes in results/*/submission.yaml or
    results/*/computational_results.json files.

    Parameters
    ----------
    base_branch : str, default='origin/main'
        Base branch to compare against

    Returns
    -------
    list[str]
        List of submission_ids (extracted from changed file paths)
        Returns empty list if git not available or no git repository

    Examples
    --------
    >>> changed = get_changed_submission_ids()
    >>> print(f"Changed submissions: {changed}")
    """
    try:
        # Run git diff to get changed files
        result = subprocess.run(
            ["git", "diff", "--name-only", f"{base_branch}...HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )

        changed_files = result.stdout.strip().split("\n")

        # Extract submission_ids from changed results paths
        submission_ids = set()
        for file_path in changed_files:
            # Match patterns like: results/SUBMISSION_ID/submission.yaml
            # or results/SUBMISSION_ID/computational_results.json
            parts = Path(file_path).parts

            submission_id = None
            if len(parts) >= 3 and parts[0] == "openfe_benchmarks" and parts[1] == "results":
                submission_id = parts[2]
            elif len(parts) >= 2 and parts[0] == "results":
                submission_id = parts[1]

            if submission_id is None:
                continue

            if parts[-1] == "submission.yaml" or "computational_results" in parts[-1]:
                submission_ids.add(submission_id)

        return sorted(list(submission_ids))

    except (subprocess.CalledProcessError, FileNotFoundError):
        # Git not available or not a git repository
        return []


def select_random_sample(
    all_ids: list[str], sample_rate: float = 0.1, seed: int = 42
) -> list[str]:
    """
    Select random sample of submissions ensuring calculation_type coverage.

    Uses get_benchmark_results(load_results=False) to load metadata, groups by
    calculation_type, ensures at least one submission per type, then randomly
    samples remaining to reach target sample_rate.

    Parameters
    ----------
    all_ids : list[str]
        All submission_ids to sample from
    sample_rate : float, default=0.1
        Fraction to sample (0.1 = 10%, reduced from 20% for CI performance)
    seed : int, default=42
        Random seed for reproducibility

    Returns
    -------
    list[str]
        Selected submission_ids (at least one per calculation_type)

    Examples
    --------
    >>> sample = select_random_sample(all_ids, sample_rate=0.1, seed=42)
    >>> # Verify reproducibility
    >>> sample2 = select_random_sample(all_ids, sample_rate=0.1, seed=42)
    >>> assert sample == sample2
    """
    # Set random seed for reproducibility
    random.seed(seed)

    # Load metadata using get_benchmark_results to get calculation_type for each submission
    calc_type_groups = {}

    for submission_id in all_ids:
        # Use get_benchmark_results to load metadata (consistent with validation)
        # If loading fails, let the exception propagate - silent failures are bad
        benchmark = get_benchmark_results(
            submission_id=submission_id, load_results=False
        )
        calc_type = benchmark.calculation_type

        if calc_type not in calc_type_groups:
            calc_type_groups[calc_type] = []

        calc_type_groups[calc_type].append(submission_id)

    # Ensure at least one per calculation_type
    selected = []
    remaining = []

    for calc_type, ids in calc_type_groups.items():
        if ids:
            # Select one representative for this type
            selected.append(ids[0])
            # Add rest to remaining pool
            remaining.extend(ids[1:])

    # Calculate how many more we need to reach sample_rate
    target_count = max(int(len(all_ids) * sample_rate), len(selected))
    additional_needed = target_count - len(selected)

    if additional_needed > 0 and remaining:
        # Randomly sample from remaining to reach target
        additional = random.sample(remaining, min(additional_needed, len(remaining)))
        selected.extend(additional)

    return sorted(selected)
