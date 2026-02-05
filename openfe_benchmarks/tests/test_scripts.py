"""Run example scripts and validate their outputs.

This test module runs every script matching ``_example*.py`` in the
``openfe_benchmarks/scripts`` directory and validates their outputs using
checks provided by the example script itself.

Contract for an example script to be testable:
- Define ``__test_config__`` (dict) at module level. Required keys:
    - ``outputs``: list[str] — paths (relative to the script) the test runner
        should check and (optionally) remove after successful validation.
    - ``custom_validate``: callable(script_dir: Path, config: dict) ->
        - a ``list[str]`` of error messages (empty list means success), or
        - ``True``/``False`` where ``False`` indicates failure, or
        - ``None`` for no-op. The runner treats non-empty lists or ``False`` as
        validation failures.

Runner flow (generic, calculation-agnostic):
1. Import the example module to read ``__test_config__``
2. Execute the script with ``runpy.run_path(..., run_name='__main__')``.
3. Check that files listed in ``outputs`` exist and are non-empty.
4. Invoke ``custom_validate(script_dir, config)`` if provided and collect
     any returned errors.
5. If validation passed and ``DEBUG`` is false, remove the output files.

For a concrete RBFE example and a richer validator that uses
``BenchmarkData``, see ``openfe_benchmarks/scripts/_example_plan_rbfe.py``.
"""

from pathlib import Path
import runpy
import multiprocessing as mp
import pytest
import shutil

TEST_DEBUG = False  # set to True to avoid deleting output files from scripts
# (just for testing the tests)

SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"


def collect_example_scripts():
    """Discover all example scripts."""
    return list(SCRIPTS_DIR.glob("_example*.py"))


def _run_and_validate(script_path: Path, queue: mp.Queue):
    """Run script, perform validations, and cleanup."""
    script_dir = script_path.parent
    errors = []

    # Import as module to get __test_config__ and optional DEBUG flag
    try:
        import importlib.util
        import os

        spec = importlib.util.spec_from_file_location("script_module", script_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        config = getattr(module, "__test_config__", {})
        debug_flag = getattr(module, "DEBUG", TEST_DEBUG)
    except Exception as e:
        queue.put((False, f"Failed to import script: {e}"))
        return

    # Run the script as __main__ from the script's directory
    try:
        import os

        original_cwd = os.getcwd()
        os.chdir(script_dir)
        try:
            runpy.run_path(str(script_path), run_name="__main__")
        finally:
            os.chdir(original_cwd)
    except Exception as e:
        queue.put((False, f"Script execution failed: {e}"))
        return

    # Validate outputs exist and are non-empty
    outputs = config.get("outputs", [])
    output_paths = []
    for out in outputs:
        path = Path(out) if Path(out).is_absolute() else script_dir / out
        output_paths.append(path)
        if not path.exists():
            errors.append(f"Missing output: {path.name}")
        elif path.is_file() and path.stat().st_size == 0:
            errors.append(f"Empty output: {path.name}")

    # Note: no calculation-specific declarative validations here. Individual
    # example scripts must provide a `custom_validate(script_dir, config)`
    # function via `__test_config__` for domain-specific checks.
    # Run custom validation if provided
    custom_validate = config.get("custom_validate")
    if custom_validate:
        try:
            result = custom_validate(script_dir, config)
            # Accept either a list of error strings, an empty list, or a boolean
            if isinstance(result, list):
                errors.extend(result)
            elif isinstance(result, bool):
                if not result:
                    errors.append("Custom validation returned False")
            elif result is not None:
                # anything else treated as a single error message
                errors.append(str(result))
        except Exception as e:
            errors.append(f"Custom validation failed: {e}")

    # Cleanup: only remove outputs when validation passed and debugging is not enabled
    if not errors and not debug_flag:
        for path in output_paths:
            try:
                if path.exists():
                    if path.is_file():
                        path.unlink()
                    elif path.is_dir():
                        shutil.rmtree(path)
            except Exception:
                pass  # Best effort cleanup

    queue.put((not errors, "\n".join(errors) if errors else "ok"))


@pytest.mark.parametrize("script_path", collect_example_scripts(), ids=lambda p: p.name)
def test_example_script(script_path: Path):
    """Run and validate a single example script.

    Each script runs in isolation with a 30s timeout.
    """
    queue: mp.Queue = mp.Queue()
    process = mp.Process(target=_run_and_validate, args=(script_path, queue))
    process.start()
    process.join(timeout=30)

    if process.is_alive():
        process.terminate()
        process.join()
        pytest.fail("Script timed out after 30s")

    try:
        success, message = queue.get_nowait()
    except Exception:
        pytest.fail("No result returned from subprocess")

    if not success:
        pytest.fail(message)
