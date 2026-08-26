"""
Filtering utilities for BenchmarkResults objects.

This module provides a flexible, extensible filtering system for benchmark submissions
with support for:
- Nested field access (via __ separator)
- Comparison operators (>=, >, <=, <) for versions, dates, and quantities
- Wildcard matching (* and ?)
- Exclusion filters (exclude_ prefix)
- Tags filtering with all/any modes

Main functions include:
    :func:`apply_filter` : Apply a filter to a source object, i.e. BenchmarkResults
    :func:`compare_values` : Compare result value against filter value using a comparison operator.
    :func:`extract_value` : Extract a field value from a source object (i.e., BenchmarkResults) or dict.
"""

from dataclasses import dataclass
from typing import Optional, Any, Union, Protocol
from datetime import date as date_type
import fnmatch
import re
import warnings
import operator as py_operator

from packaging.version import Version, InvalidVersion
import pint

__all__ = [
    "apply_filter",
    "parse_filter",
    "extract_value",
    "match_value",
    "compare_values",
    "ParsedFilter",
]


# Type aliases
FilterValue = Union[str, list, Any]


# ============================================================================
# Comparison Strategy Pattern
# ============================================================================


class ComparisonStrategy(Protocol):
    """Protocol for value comparison strategies."""

    def can_handle(self, value: Any) -> bool:
        """Return True if this strategy can handle the given value type."""
        ...

    def compare(self, left: Any, right: Any, operator: str) -> bool:
        """Perform comparison using the given operator."""
        ...


@dataclass
class DateComparisonStrategy:
    """Handle datetime.date comparisons with ISO string coercion."""

    def can_handle(self, value: Any) -> bool:
        return isinstance(value, date_type)

    def compare(self, left: date_type, right: Any, operator: str) -> bool:
        # Coerce right operand to date if string
        if isinstance(right, str):
            try:
                right = date_type.fromisoformat(right)
            except (ValueError, AttributeError) as e:
                raise ValueError(
                    f"Invalid date filter value '{right}'. "
                    f"Provide a valid ISO format date (YYYY-MM-DD): {e}"
                )

        return _apply_operator(left, right, operator)


@dataclass
class QuantityComparisonStrategy:
    """Handle pint.Quantity comparisons with automatic unit conversion."""

    def can_handle(self, value: Any) -> bool:
        return isinstance(value, pint.Quantity)

    def compare(self, left: pint.Quantity, right: Any, operator: str) -> bool:
        # Coerce right operand to Quantity
        if isinstance(right, str):
            try:
                ureg = left._REGISTRY
                right = ureg.Quantity(right)
            except (ValueError, pint.errors.UndefinedUnitError, AttributeError) as e:
                raise ValueError(
                    f"Invalid quantity filter value '{right}'. "
                    f"Provide a valid quantity string with units: {e}"
                )
        elif not isinstance(right, pint.Quantity):
            raise ValueError(
                f"Invalid filter type '{type(right).__name__}' for Quantity comparison. "
                "Provide a quantity string or pint.Quantity value."
            )

        try:
            return _apply_operator(left, right, operator)
        except pint.errors.DimensionalityError as e:
            raise ValueError(f"Incompatible units for Quantity comparison: {e}") from e


@dataclass
class VersionComparisonStrategy:
    """Handle semantic version comparisons."""

    def can_handle(self, value: Any) -> bool:
        # Always try as fallback (returns False to not auto-match)
        return False

    def compare(self, left: Any, right: Any, operator: str) -> bool:
        try:
            left_ver = Version(str(left))
            right_ver = Version(str(right))
            return _apply_operator(left_ver, right_ver, operator)
        except (InvalidVersion, TypeError) as e:
            raise ValueError(
                f"Invalid version filter value. Cannot parse '{left}' or '{right}' as semantic version. "
                f"Comparison operators with non-version fields are not supported: {e}"
            )


# Registry of comparison strategies (order matters - first match wins)
_COMPARISON_STRATEGIES: list[ComparisonStrategy] = [
    DateComparisonStrategy(),
    QuantityComparisonStrategy(),
    VersionComparisonStrategy(),  # Fallback for version strings
]


def _apply_operator(left: Any, right: Any, operator: str) -> bool:
    """Apply comparison operator to two values."""
    operator_map = {
        "<": py_operator.lt,
        "<=": py_operator.le,
        ">": py_operator.gt,
        ">=": py_operator.ge,
    }
    compare_func = operator_map.get(operator)
    if compare_func is None:
        raise ValueError(f"Unknown comparison operator: {operator}")

    return compare_func(left, right)


def _warn_if_nonsemantic_comparison(
    result_value: Any, operator: str, field_name: str
) -> None:
    """Warn when comparison operators are used on fields that may not compare semantically."""
    semantic_fields = {
        "date",
        "openfe_version",
        "openmm_version",
        "openff_toolkit_version",
        "pontibus_version",
    }
    if (
        field_name
        and field_name not in semantic_fields
        and not field_name.endswith("_version")
        and not isinstance(result_value, pint.Quantity)
    ):
        warnings.warn(
            f"Comparison operator '{operator}' used with field '{field_name}'. "
            "Comparison operators are designed for date and version fields. "
            "Results may not be semantically meaningful for other field types.",
            UserWarning,
            stacklevel=5,
        )


# ============================================================================
# Filter Parsing
# ============================================================================


@dataclass
class ParsedFilter:
    """
    Parsed and normalized filter specification.

    Attributes
    ----------
    field_name : str
        Field name without exclude_ prefix
    filter_value : Any
        Filter value after operator extraction
    comparison_op : Optional[str]
        Comparison operator if present (>=, >, <=, <)
    negate : bool
        Whether this is an exclusion filter (had exclude_ prefix)
    is_nested : bool
        Whether this uses nested field access (contains __)
    """

    field_name: str
    filter_value: Any
    comparison_op: Optional[str]
    negate: bool
    is_nested: bool


def parse_filter(filter_key: str, filter_val: FilterValue) -> ParsedFilter:
    """
    Parse and normalize a filter specification.

    Handles:
    - exclude_ prefix for negation
    - Inline comparison operators (>=, >, <=, <)
    - Nested field detection (__)

    Parameters
    ----------
    filter_key : str
        Filter key (may have exclude_ prefix and/or __ for nesting)
    filter_val : FilterValue
        Filter value (may have inline comparison operator)

    Returns
    -------
    ParsedFilter
        Normalized filter components

    Examples
    --------
    >>> parse_filter("exclude_tags", ["test"])
    ParsedFilter(field_name='tags', filter_value=['test'], comparison_op=None, negate=True, is_nested=False)

    >>> parse_filter("openfe_version", ">=1.0.0")
    ParsedFilter(field_name='openfe_version', filter_value='1.0.0', comparison_op='>=', negate=False, is_nested=False)

    >>> parse_filter("protocol_settings__temperature", "298 K")
    ParsedFilter(field_name='protocol_settings__temperature', filter_value='298 K', comparison_op=None, negate=False, is_nested=True)
    """
    # Extract negation
    negate = filter_key.startswith("exclude_")
    if negate:
        filter_key = filter_key[len("exclude_") :]

    # Parse inline comparison operator
    comparison_op = None
    if isinstance(filter_val, str):
        comparison_op, filter_val = _parse_inline_operator(filter_val)

    # Check for nested field access
    is_nested = "__" in filter_key

    return ParsedFilter(
        field_name=filter_key,
        filter_value=filter_val,
        comparison_op=comparison_op,
        negate=negate,
        is_nested=is_nested,
    )


def _parse_inline_operator(value: str) -> tuple[Optional[str], str]:
    """
    Parse inline comparison operator from a filter value string.

    Parameters
    ----------
    value : str
        Filter value that may contain an inline operator

    Returns
    -------
    tuple[Optional[str], str]
        (operator, remaining_value) where operator is one of >=, >, <=, <, or None

    Examples
    --------
    >>> _parse_inline_operator(">=2026-01-01")
    ('>=', '2026-01-01')
    >>> _parse_inline_operator("2026-01-01")
    (None, '2026-01-01')
    """
    match = re.match(r"^(>=|>|<=|<)(.+)$", value)
    if not match:
        return None, value

    op_str, remaining = match.groups()
    return op_str, remaining.strip()


# ============================================================================
# Value Extraction
# ============================================================================


def extract_value(
    source: Any,
    field_name: str,
    is_nested: bool,
) -> Optional[Any]:
    """
    Extract a field value from a source object or dict.

    Parameters
    ----------
    source : Any
        Source object (dict or object)
    field_name : str
        Field name (may include __ for nesting)
    is_nested : bool
        Whether this is a nested field access

    Returns
    -------
    Optional[Any]
        Extracted value, or None if not found

    Examples
    --------
    >>> result = {"system_name": "tyk2", "ligand_name": "jmc_23"}
    >>> extract_value(result, "system_name", is_nested=False)
    'tyk2'
    """
    if is_nested:
        return _get_nested_value(source, field_name)

    # Direct field access
    if isinstance(source, dict):
        return source.get(field_name)
    else:
        return getattr(source, field_name, None)


def _get_nested_value(source: Any, nested_key: str) -> Optional[Any]:
    """
    Extract value from nested structure using __ separator.

    Handles dicts, objects, and lists of dicts/objects.

    Parameters
    ----------
    source : Any
        Source object or dictionary
    nested_key : str
        Nested key with __ separator (e.g., 'protocol_settings__temperature')

    Returns
    -------
    Optional[Any]
        The nested value if found, None otherwise

    Examples
    --------
    >>> result = {'protocol_settings': {'temperature': '298.15 K'}}
    >>> _get_nested_value(result, 'protocol_settings__temperature')
    '298.15 K'

    >>> print(_get_nested_value(result, 'protocol_settings__missing'))
    None

    >>> # Handles lists of dicts
    >>> result = {'protocol_settings': [{'temp': 298}, {'temp': 300}]}
    >>> _get_nested_value(result, 'protocol_settings__temp')
    [298, 300]
    """
    keys = nested_key.split("__")
    value = source

    for key in keys:
        if isinstance(value, dict):
            value = value.get(key)
        elif isinstance(value, list):
            # Collect values from list of dicts/objects
            next_values = []
            for item in value:
                if isinstance(item, dict):
                    nested = item.get(key)
                else:
                    nested = getattr(item, key, None)

                if nested is not None:
                    next_values.append(nested)

            value = next_values if next_values else None
        else:
            value = getattr(value, key, None)

        if value is None:
            return None

    return value


# ============================================================================
# Value Matching (Equality, Wildcards, Lists)
# ============================================================================


def match_value(
    result_value: Any,
    filter_val: FilterValue,
    filter_key: str,
    tags_mode: str,
) -> bool:
    """
    Check if a result value matches a filter value.

    Handles:
    - Tags with all/any mode
    - List filter values (OR logic)
    - Wildcard matching
    - Exact matching

    Parameters
    ----------
    result_value : Any
        Value from result to check
    filter_val : FilterValue
        Filter value (single value, list, or wildcard pattern)
    filter_key : str
        Filter key name (for special tags handling)
    tags_mode : str
        Mode for tags filtering: 'all' (AND) or 'any' (OR)

    Returns
    -------
    bool
        True if value matches filter

    Raises
    ------
    ValueError
        If tags_mode is not 'all' or 'any'

    Examples
    --------
    >>> match_value(["rbfe", "jacs"], ["rbfe"], "tags", "any")
    True
    >>> match_value("tyk2", "tyk*", "system_name", "all")
    True
    >>> match_value("1.0.0", ["1.0.0", "2.0.0"], "version", "all")
    True
    """
    # Special handling for tags field with list filter_val
    if filter_key == "tags" and isinstance(filter_val, list):
        result_tags = result_value if isinstance(result_value, list) else [result_value]

        if tags_mode == "all":
            return all(tag in result_tags for tag in filter_val)
        elif tags_mode == "any":
            return any(tag in result_tags for tag in filter_val)
        else:
            raise ValueError(f"Invalid tags_mode: {tags_mode}. Must be 'all' or 'any'")

    # OR logic for list filter values (non-tags fields)
    if isinstance(filter_val, list):
        return any(_match_single_value(result_value, fv) for fv in filter_val)

    # Single value match
    return _match_single_value(result_value, filter_val)


def _match_single_value(result_value: Any, filter_val: Any) -> bool:
    """
    Check if a single result value matches a single filter value.

    Handles wildcards and list-valued result fields.

    Parameters
    ----------
    result_value : Any
        Value from result to check
    filter_val : Any
        Single filter value (not a list)

    Returns
    -------
    bool
        True if value matches

    Examples
    --------
    >>> _match_single_value("tyk2", "tyk*")
    True
    >>> _match_single_value(["a", "b", "c"], "b")
    True
    >>> _match_single_value("exact", "exact")
    True
    """
    # Handle list-valued result fields (check if ANY element matches)
    if isinstance(result_value, list):
        return any(_match_single_value(item, filter_val) for item in result_value)

    # Wildcard matching for strings
    if isinstance(filter_val, str) and ("*" in filter_val or "?" in filter_val):
        return fnmatch.fnmatch(str(result_value), filter_val)

    # Exact match
    return result_value == filter_val


# ============================================================================
# Value Comparison (Versions, Dates, Quantities)
# ============================================================================


def compare_values(
    result_value: Any,
    filter_val: Any,
    operator: str,
    field_name: str = "",
) -> bool:
    """
    Compare result value against filter value using a comparison operator.

    Uses strategy pattern to handle different value types:
    - datetime.date objects
    - pint.Quantity objects
    - Semantic versions
    - Fallback to string comparison

    Parameters
    ----------
    result_value : Any
        Value from result to compare
    filter_val : Any
        Filter value to compare against
    operator : str
        Comparison operator: '<', '<=', '>', or '>='
    field_name : str, optional
        Name of the field being compared (for warnings)

    Returns
    -------
    bool
        True if comparison succeeds

    Raises
    ------
    ValueError
        If comparison cannot be performed (e.g., incompatible units)

    Examples
    --------
    >>> from datetime import date
    >>> compare_values(date(2026, 1, 15), "2026-01-01", ">=", "date")
    True

    >>> compare_values("1.0.0", "2.0.0", "<", "openfe_version")
    True

    >>> # With pint quantities
    >>> import pint
    >>> ureg = pint.UnitRegistry()
    >>> temp = ureg.Quantity(298.15, 'kelvin')
    >>> compare_values(temp, "300 K", "<")
    True
    """
    # Handle list-valued nested fields (match if ANY element satisfies comparison)
    if isinstance(result_value, list):
        return any(
            compare_values(item, filter_val, operator, field_name)
            for item in result_value
        )

    _warn_if_nonsemantic_comparison(result_value, operator, field_name)

    # Try each comparison strategy
    for strategy in _COMPARISON_STRATEGIES:
        if strategy.can_handle(result_value):
            return strategy.compare(result_value, filter_val, operator)

    # Fallback: version comparison strategy (always tries, falls back to string)
    return VersionComparisonStrategy().compare(result_value, filter_val, operator)


# ============================================================================
# Main Filter Application
# ============================================================================


def apply_filter(
    source: Any,
    filter_key: str,
    filter_val: FilterValue,
    tags_mode: str = "all",
) -> bool:
    """
    Apply a filter to a source object.

    This is the main entry point for filtering logic. Works for both
    submission-level and result-entry filtering.

    Parameters
    ----------
    source : Any
        Source to filter (BenchmarkResults or result dict)
    filter_key : str
        Filter key (may have exclude_ prefix, may have __ for nested fields)
    filter_val : FilterValue
        Filter value (single value, list, or comparison string)
    tags_mode : str, default='all'
        Mode for tags filtering: 'all' (AND) or 'any' (OR)

    Returns
    -------
    bool
        True if source matches filter

    Examples
    --------
    >>> from dataclasses import dataclass
    >>> @dataclass
    ... class Submission:
    ...     tags: list[str]
    ...     date: str
    >>> submission = Submission(tags=["rbfe", "jacs"], date="2026-01-15")

    >>> # Direct field matching
    >>> apply_filter(submission, "tags", ["rbfe"], tags_mode="any")
    True

    >>> # Comparison operator
    >>> apply_filter(submission, "date", ">=2026-01-01")
    True

    >>> # Exclusion filter
    >>> apply_filter(submission, "exclude_tags", ["test"])
    True
    """
    # Parse filter specification
    parsed = parse_filter(filter_key, filter_val)

    # Extract value from source
    value = extract_value(source, parsed.field_name, parsed.is_nested)

    # If field not found, return False (with negation if applicable)
    if value is None:
        return parsed.negate

    # Apply comparison or match logic
    if parsed.comparison_op is not None:
        matched = compare_values(
            value, parsed.filter_value, parsed.comparison_op, parsed.field_name
        )
    else:
        matched = match_value(value, parsed.filter_value, parsed.field_name, tags_mode)

    # Apply negation if needed
    return not matched if parsed.negate else matched
