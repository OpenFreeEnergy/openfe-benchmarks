"""
Module for loading and filtering computational benchmark results from submission.yaml files.

The BenchmarkResults class is a data container that loads on initialization.
Filtering functions are separate from the class, following functional programming principles.
"""

from pathlib import Path
from dataclasses import dataclass
import yaml
import json
import bz2
import fnmatch
import logging

from gufe.tokenization import JSON_HANDLER

# Get logger for this module - will inherit configuration from parent
logger = logging.getLogger(__name__)

__all__ = [
    "BenchmarkResults",
    "filter_results",
]

_RESULTS_DIR = (
    Path(__file__).resolve().parent.parent.parent / "openfe_benchmarks" / "results"
)


@dataclass
class BenchmarkResults:
    """
    Represents computational results from a submission.yaml file.

    Initialize with a submission_id to load the results. Use filter_results()
    function to filter the raw computational data.

    Attributes
    ----------
    submission_id : str
        Unique identifier from submission.yaml
    title : str
        Descriptive title
    calculation_type : str
        Type of calculation (rbfe, asfe, etc.)
    tags : list[str]
        List of submission tags
    metadata : dict
        All metadata from submission.yaml (authors, date, forcefields, etc.)
    raw_results : dict | None
        Raw computational_results.json loaded as nested dict, or None if load_results=False
    results_file : Path
        Path to the computational_results.json file
    submission_file : Path
        Path to the submission.yaml file

    Examples
    --------
    >>> # Load results with computational data
    >>> results = BenchmarkResults(submission_id='2026-03-18-openmm-840-qa-testing')
    >>>
    >>> # Fast YAML-only load (for CI validation)
    >>> results = BenchmarkResults(submission_id='2026-03-18-openmm-840-qa-testing', load_results=False)
    >>>
    >>> # Filter the results
    >>> rbfe_results = filter_results(results, tags='rbfe')
    >>> tyk2_results = filter_results(results, system_name='tyk2')
    """

    submission_id: str
    title: str
    calculation_type: str
    tags: list[str]
    metadata: dict
    raw_results: dict | None
    results_file: Path
    submission_file: Path

    def __init__(
        self,
        submission_id: str,
        load_results: bool = True,
        results_dir: Path | None = None,
    ):
        """
        Initialize BenchmarkResults by loading from submission_id.

        Parameters
        ----------
        submission_id : str
            Unique submission identifier (e.g., '2026-03-18-openmm-840-qa-testing')
        load_results : bool, default=True
            If True, load computational_results.json. If False, only load YAML metadata
            (for fast CI validation). When False, raw_results will be None.
        results_dir : Path, optional
            Results directory to search. Defaults to openfe_benchmarks/results/

        Raises
        ------
        FileNotFoundError
            If submission_id directory or submission.yaml not found
        ValueError
            If submission.yaml missing required fields or submission_id mismatch
        """
        # Auto-discover results directory
        if results_dir is None:
            results_dir = _RESULTS_DIR
        else:
            results_dir = Path(results_dir)

        # Construct paths
        submission_dir = results_dir / submission_id
        submission_file = submission_dir / "submission.yaml"

        # Validate paths exist
        if not submission_dir.exists():
            raise FileNotFoundError(f"Submission not found: {submission_id}")

        if not submission_file.exists():
            raise FileNotFoundError(f"submission.yaml not found for: {submission_id}")

        # Load YAML
        try:
            with open(submission_file, "r") as f:
                yaml_data = yaml.safe_load(f)
        except Exception as e:
            raise ValueError(f"Error reading submission.yaml for {submission_id}: {e}")

        # Validate required fields
        required_fields = [
            "submission_id",
            "title",
            "calculation_type",
            "tags",
            "results",
        ]
        missing_fields = [field for field in required_fields if field not in yaml_data]
        if missing_fields:
            raise ValueError(
                f"Invalid submission.yaml: missing required fields: {missing_fields}"
            )

        # Validate submission_id matches
        yaml_submission_id = yaml_data["submission_id"]
        if yaml_submission_id != submission_id:
            raise ValueError(
                f"submission_id mismatch: YAML has '{yaml_submission_id}', expected '{submission_id}'"
            )

        # Extract basic fields
        self.submission_id = yaml_data["submission_id"]
        self.title = yaml_data["title"]
        self.calculation_type = yaml_data["calculation_type"]
        self.tags = (
            yaml_data["tags"]
            if isinstance(yaml_data["tags"], list)
            else [yaml_data["tags"]]
        )
        self.metadata = yaml_data

        # Construct results file path
        results_filename = yaml_data["results"]
        self.results_file = submission_dir / results_filename
        self.submission_file = submission_file

        # Conditionally load computational results
        self.raw_results = None
        if load_results:
            # Validate results file exists
            if not self.results_file.exists():
                raise FileNotFoundError(f"Results file not found: {self.results_file}")

            # Load JSON with bz2 support
            open_func = bz2.open if ".bz2" in str(self.results_file) else open
            try:
                with open_func(self.results_file, "rt") as handle:
                    self.raw_results = json.load(handle, cls=JSON_HANDLER.decoder)
            except Exception as e:
                raise ValueError(f"Error loading results file {self.results_file}: {e}")

        logger.debug(
            f"Loaded BenchmarkResults: {self.submission_id} (load_results={load_results})"
        )

    def __repr__(self):
        """Return string representation."""
        return f"BenchmarkResults(submission_id='{self.submission_id}', title='{self.title}')"


# Standalone filtering functions


def filter_results(
    benchmark_results: BenchmarkResults, tags_mode: str = "all", **filters
) -> list[dict]:
    """
    Filter raw computational results by any combination of fields.

    Supports:
    - Top-level metadata: tags='rbfe', calculation_type='rbfe'
    - Nested fields: protocol_settings__temperature='298.15 K' (use __ for nesting)
    - Result fields: system_group='jacs_set', system_name='tyk2'
    - Wildcards: system_name='*tyk2*' (uses fnmatch)
    - OR logic within field: pass list for ANY match (e.g., system_name=['tyk2', 'thrombin'])
    - NOT logic: use exclude_ prefix (e.g., exclude_tags=['deprecated'])
    - AND logic between fields: all filter conditions must match

    Parameters
    ----------
    benchmark_results : BenchmarkResults
        BenchmarkResults instance to filter
    tags_mode : {'all', 'any'}, default='all'
        When filtering by multiple tags:
        - 'all': result must have ALL specified tags (AND logic, default)
        - 'any': result must have ANY specified tag (OR logic)
        Ignored if tags filter is not provided or is a single value.
    **filters : dict
        Field=value pairs to filter by. Values can be:
        - Single value: exact match (or wildcard if contains * or ?)
        - List: matches if ANY value matches (OR logic) - EXCEPT tags which respects tags_mode
        - Use exclude_ prefix for NOT logic (e.g., exclude_tags=['test'])

    Returns
    -------
    list[dict]
        Filtered result dictionaries

    Raises
    ------
    ValueError
        If raw_results is None (load_results=False was used)

    Examples
    --------
    >>> results = BenchmarkResults(submission_id='2026-03-18-openmm-840-qa-testing')
    >>>
    >>> # Exact match
    >>> rbfe_results = filter_results(results, tags='rbfe')
    >>>
    >>> # Tags AND logic (default): must have ALL tags
    >>> validated = filter_results(results, tags=['rbfe', 'openfe', 'validation'])
    >>>
    >>> # Tags OR logic: must have ANY tag
    >>> any_calc = filter_results(results, tags=['rbfe', 'asfe'], tags_mode='any')
    >>>
    >>> # OR within non-tag fields (list = ANY match)
    >>> multi_system = filter_results(results, system_name=['tyk2', 'thrombin'])
    >>>
    >>> # NOT logic (exclude_ prefix)
    >>> no_deprecated = filter_results(results, exclude_tags='deprecated')
    >>> no_test = filter_results(results, exclude_tags=['deprecated', 'test'], tags_mode='any')
    >>>
    >>> # Complex: multiple tags (AND), system OR, exclude
    >>> filtered = filter_results(
    ...     results,
    ...     tags=['rbfe', 'openfe'],           # has BOTH rbfe AND openfe
    ...     system_name=['tyk2', 'thrombin'],  # AND (tyk2 OR thrombin)
    ...     exclude_tags=['deprecated']        # AND NOT deprecated
    ... )
    >>>
    >>> # Nested fields and wildcards
    >>> lambda11 = filter_results(results, protocol_settings__lambda_windows='11')
    >>> tyk2_wildcard = filter_results(results, system_name='*tyk2*')
    """
    if benchmark_results.raw_results is None:
        raise ValueError(
            "Cannot filter results: raw_results is None. "
            "Initialize with load_results=True to access computational data."
        )

    # Collect all results (dg + ddg)
    all_results = []
    if "dg" in benchmark_results.raw_results:
        all_results.extend(benchmark_results.raw_results["dg"])
    if "ddg" in benchmark_results.raw_results:
        all_results.extend(benchmark_results.raw_results["ddg"])

    # Apply filters
    filtered = []
    for result in all_results:
        # Check all filters (AND logic between different filters)
        if all(
            _match_filter(result, benchmark_results.metadata, key, value, tags_mode)
            for key, value in filters.items()
        ):
            filtered.append(result)

    return filtered


def _match_filter(
    result: dict, metadata: dict, filter_key: str, filter_val, tags_mode: str = "all"
) -> bool:
    """
    Check if a result matches a filter predicate.

    Parameters
    ----------
    result : dict
        Result dictionary to check
    metadata : dict
        Metadata from BenchmarkResults (for fallback field access)
    filter_key : str
        Filter key (may have exclude_ prefix, may have __ for nested access)
    filter_val : any
        Filter value (may be single value or list for OR logic)
    tags_mode : str
        Mode for tags filtering ('all' or 'any')

    Returns
    -------
    bool
        True if result matches filter, False otherwise
    """
    # Handle NOT logic (exclude_ prefix)
    negate = False
    if filter_key.startswith("exclude_"):
        negate = True
        filter_key = filter_key[8:]  # Remove 'exclude_' prefix

    # Handle nested field access (e.g., protocol_settings__lambda_windows)
    if "__" in filter_key:
        # Split on __ and traverse nested dict
        keys = filter_key.split("__")
        value = result
        for key in keys:
            if isinstance(value, dict):
                value = value.get(key, None)
            else:
                value = None
                break

        if value is None:
            # Field not found in result
            result_matched = False
        else:
            result_matched = _match_value(value, filter_val, filter_key, tags_mode)
    else:
        # Direct field access
        # First try result dict, then fall back to metadata for top-level fields
        if filter_key in result:
            value = result[filter_key]
        elif filter_key in metadata:
            value = metadata[filter_key]
        else:
            # Field not found
            result_matched = False
            return not result_matched if negate else result_matched

        result_matched = _match_value(value, filter_val, filter_key, tags_mode)

    # Apply negation if needed
    return not result_matched if negate else result_matched


def _match_value(result_value, filter_val, filter_key: str, tags_mode: str) -> bool:
    """
    Check if a result value matches a filter value.

    Parameters
    ----------
    result_value : any
        Value from result to check
    filter_val : any
        Filter value (may be single value or list)
    filter_key : str
        Filter key (for special handling of 'tags')
    tags_mode : str
        Mode for tags filtering ('all' or 'any')

    Returns
    -------
    bool
        True if value matches filter, False otherwise
    """
    # Special handling for tags field with list filter_val
    if filter_key == "tags" and isinstance(filter_val, list):
        # Ensure result_value is a list
        result_tags = result_value if isinstance(result_value, list) else [result_value]

        if tags_mode == "all":
            # AND logic: result must have ALL filter tags
            return all(tag in result_tags for tag in filter_val)
        elif tags_mode == "any":
            # OR logic: result must have ANY filter tag
            return any(tag in result_tags for tag in filter_val)
        else:
            raise ValueError(f"Invalid tags_mode: {tags_mode}. Must be 'all' or 'any'")

    # OR logic for list filter values (non-tags fields)
    if isinstance(filter_val, list):
        # Check if ANY filter value matches
        return any(_match_single_value(result_value, fv) for fv in filter_val)
    else:
        # Single value match
        return _match_single_value(result_value, filter_val)


def _match_single_value(result_value, filter_val) -> bool:
    """
    Check if a single result value matches a single filter value.

    Handles wildcards and list-valued result fields.

    Parameters
    ----------
    result_value : any
        Value from result to check
    filter_val : any
        Single filter value (not a list)

    Returns
    -------
    bool
        True if value matches, False otherwise
    """
    # Handle list-valued result fields (check if ANY element matches)
    if isinstance(result_value, list):
        return any(_match_single_value(item, filter_val) for item in result_value)

    # Wildcard matching for strings
    if isinstance(filter_val, str) and ("*" in filter_val or "?" in filter_val):
        return fnmatch.fnmatch(str(result_value), filter_val)

    # Exact match
    return result_value == filter_val
