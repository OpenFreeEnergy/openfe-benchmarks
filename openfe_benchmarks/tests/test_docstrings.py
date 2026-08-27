"""
Test docstring examples using Python's built-in doctest module.

Most tests validate docstring syntax without execution, avoiding data dependencies
while ensuring examples remain valid Python.

The test_docstring_examples_with_real_data() test executes all docstring examples
with real data files and runs unconditionally with the full test suite.
"""

import ast
import doctest
import pytest


def test_benchmark_results_docstring_syntax():
    """
    Test that all docstring code examples are syntactically valid Python.

    Uses doctest.DocTestFinder to extract examples, then validates syntax
    with ast.parse (without execution).
    """
    import openfe_benchmarks.results._benchmark_results as br_module

    finder = doctest.DocTestFinder()
    total_examples = 0
    invalid_examples = []
    tests_found = []

    # Find all doctests in the module
    for test in finder.find(br_module, module=br_module):
        test_examples = len(test.examples)
        if test_examples > 0:
            tests_found.append(f"{test.name}: {test_examples} examples")

        for example in test.examples:
            total_examples += 1
            try:
                # Validate syntax without executing
                ast.parse(example.source)
            except SyntaxError as e:
                invalid_examples.append(
                    {
                        "test": test.name,
                        "line": example.lineno,
                        "source": example.source,
                        "error": str(e),
                    }
                )

    # Always print summary for visibility
    print(f"\n_benchmark_results.py: Validated {total_examples} code examples")
    print(f"Found in: {', '.join(tests_found[:5])}")
    if len(tests_found) > 5:
        print(f"  ... and {len(tests_found) - 5} more")

    # Ensure we found a reasonable number of examples (sanity check)
    assert total_examples >= 50, (
        f"Expected at least 50 examples, found only {total_examples}. "
        "Examples may not be detected properly."
    )

    # Report results
    if invalid_examples:
        error_msg = f"Found {len(invalid_examples)} syntax errors in {total_examples} examples:\n"
        for err in invalid_examples:
            error_msg += f"\n  {err['test']} (line {err['line']}): {err['error']}\n"
            error_msg += f"    {err['source'][:60]}...\n"
        pytest.fail(error_msg)


def test_filtering_docstring_syntax():
    """Test that filtering module docstring examples are syntactically valid."""
    import openfe_benchmarks.results._filtering as filtering_module

    finder = doctest.DocTestFinder()
    total_examples = 0
    invalid_examples = []
    tests_found = []

    # Find all doctests in the module
    for test in finder.find(filtering_module, module=filtering_module):
        test_examples = len(test.examples)
        if test_examples > 0:
            tests_found.append(f"{test.name}: {test_examples} examples")

        for example in test.examples:
            total_examples += 1
            try:
                # Validate syntax without executing
                ast.parse(example.source)
            except SyntaxError as e:
                invalid_examples.append(
                    {
                        "test": test.name,
                        "line": example.lineno,
                        "source": example.source,
                        "error": str(e),
                    }
                )

    # Always print summary for visibility
    print(f"\n_filtering.py: Validated {total_examples} code examples")
    print(f"Found in: {', '.join(tests_found[:5])}")
    if len(tests_found) > 5:
        print(f"  ... and {len(tests_found) - 5} more")

    # Ensure we found a reasonable number of examples (sanity check)
    assert total_examples >= 25, (
        f"Expected at least 25 examples, found only {total_examples}. "
        "Examples may not be detected properly."
    )

    # Report results
    if invalid_examples:
        error_msg = f"Found {len(invalid_examples)} syntax errors in {total_examples} examples:\n"
        for err in invalid_examples:
            error_msg += f"\n  {err['test']} (line {err['line']}): {err['error']}\n"
            error_msg += f"    {err['source'][:60]}...\n"
        pytest.fail(error_msg)


def test_docstring_example_imports():
    """
    Verify that examples use correct public API patterns.

    Ensures examples use get_benchmark_results() instead of direct
    BenchmarkResults() instantiation.
    """
    import openfe_benchmarks.results._benchmark_results as br_module

    finder = doctest.DocTestFinder()
    issues = []

    # Check BenchmarkResults class examples
    for test in finder.find(br_module.BenchmarkResults):
        for example in test.examples:
            source = example.source
            # Should use get_benchmark_results(), not direct instantiation
            if "BenchmarkResults(" in source and "get_benchmark_results" not in source:
                # Exception: allow for specific internal examples
                if "BenchmarkResults(**yaml_data)" not in source:
                    issues.append(
                        f"Line {example.lineno}: Direct instantiation instead of get_benchmark_results()\n"
                        f"  Found: {source.strip()[:60]}..."
                    )

    if issues:
        pytest.fail("API usage issues in docstring examples:\n" + "\n".join(issues))


def test_docstring_completeness():
    """
    Verify that key functions and classes have docstrings with examples.

    Checks that:
    1. Docstring exists
    2. Contains 'Examples' section
    3. Has runnable code examples (>>> markers)
    """
    import openfe_benchmarks.results._benchmark_results as br_module

    required_docstrings = [
        ("filter_results", br_module.filter_results),
        ("get_benchmark_results", br_module.get_benchmark_results),
        ("BenchmarkResults", br_module.BenchmarkResults),
    ]

    finder = doctest.DocTestFinder()
    missing_or_incomplete = []

    for name, obj in required_docstrings:
        if not obj.__doc__:
            missing_or_incomplete.append(f"{name}: Missing docstring")
        elif "Examples" not in obj.__doc__:
            missing_or_incomplete.append(f"{name}: Missing 'Examples' section")
        else:
            # Use doctest to find examples
            tests = finder.find(obj, name=name)
            example_count = sum(len(test.examples) for test in tests)
            if example_count == 0:
                missing_or_incomplete.append(f"{name}: No runnable code examples")

    if missing_or_incomplete:
        pytest.fail(
            "Documentation completeness issues:\n  "
            + "\n  ".join(missing_or_incomplete)
        )


def test_docstring_validation_summary():
    """
    Report total docstring examples validated across both modules.

    Prints a summary of validation coverage and asserts the combined total
    meets the sum of individual module minimums (50 + 25 = 75).
    """
    import openfe_benchmarks.results._benchmark_results as br_module
    import openfe_benchmarks.results._filtering as filtering_module

    finder = doctest.DocTestFinder()

    # Count examples in both modules
    br_examples = sum(
        len(test.examples) for test in finder.find(br_module, module=br_module)
    )

    filtering_examples = sum(
        len(test.examples)
        for test in finder.find(filtering_module, module=filtering_module)
    )

    total = br_examples + filtering_examples

    # Always print summary
    print(f"\n{'=' * 60}")
    print("DOCSTRING VALIDATION SUMMARY")
    print(f"{'=' * 60}")
    print(f"_benchmark_results.py: {br_examples} examples")
    print(f"_filtering.py:         {filtering_examples} examples")
    print(f"{'=' * 60}")
    print(f"TOTAL VALIDATED:       {total} code examples")
    print(f"{'=' * 60}")

    # Verify combined total meets sum of individual module minimums
    assert total >= 75, (
        f"Expected at least 75 total examples (50 + 25 from individual modules), "
        f"found {total}. Examples may not be detected properly."
    )


def test_docstring_examples_with_real_data():
    """
    Execute docstring examples with real data files.

    Requires actual submission directories and results files.
    All doctests must pass - failures indicate code/API changes that broke examples.
    """
    import openfe_benchmarks.results._benchmark_results as br_module
    import openfe_benchmarks.results._filtering as filtering_module

    failures = {}

    # Run full doctest on both modules
    for module_name, module in [
        ("_benchmark_results", br_module),
        ("_filtering", filtering_module),
    ]:
        results = doctest.testmod(
            module,
            optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE,
            verbose=False,  # Only show failures, not all examples
        )

        if results.failed > 0:
            failures[module_name] = (results.attempted, results.failed)

    # Assert all doctests pass
    if failures:
        error_msg = "\n" + "=" * 70 + "\n"
        error_msg += "DOCTEST FAILURES\n"
        error_msg += "=" * 70 + "\n"
        for module_name, (attempted, failed) in failures.items():
            error_msg += f"\n{module_name}: {failed}/{attempted} examples failed\n"
        error_msg += "\nDoctests must pass to ensure examples stay in sync with code.\n"
        error_msg += "Run with pytest -v to see detailed failure output."
        pytest.fail(error_msg)


if __name__ == "__main__":
    """
    Direct execution for quick validation.

    Usage:
        python -m openfe_benchmarks.tests.test_docstrings
    """
    import openfe_benchmarks.results._benchmark_results as br_module
    import openfe_benchmarks.results._filtering as filtering_module

    print("Docstring Example Validation")
    print("=" * 50)

    # Use doctest's DocTestFinder to count examples
    finder = doctest.DocTestFinder()

    br_tests = finder.find(br_module, module=br_module)
    br_examples = sum(len(test.examples) for test in br_tests)

    filtering_tests = finder.find(filtering_module, module=filtering_module)
    filtering_examples = sum(len(test.examples) for test in filtering_tests)

    print(f"\n_benchmark_results.py: {br_examples} code examples")
    print(f"_filtering.py: {filtering_examples} code examples")

    print("\nTo run tests:")
    print("  pytest openfe_benchmarks/tests/test_docstrings.py -v")
