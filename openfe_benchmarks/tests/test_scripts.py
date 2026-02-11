"""Run each `_example*.py` as a standalone script.

This test discovers example files in `openfe_benchmarks/scripts`, executes
each with `runpy.run_path(..., run_name="__main__")` in an isolated
process (120s timeout), and fails if the script raises or times out. The
test performs no output validation or cleanup.
"""

from pathlib import Path
import runpy
import multiprocessing as mp
import pytest
import os

TEST_DEBUG = False  # set to True to avoid leaving outputs behind
SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"


def collect_example_scripts():
    """Discover all example scripts matching ``_example*.py``."""
    return list(SCRIPTS_DIR.glob("_example*.py"))


def _run_script_as_main(script_path: Path, queue: mp.Queue):
    """Execute the script with runpy as ``__main__`` and report exceptions.

    This intentionally performs no validation or cleanup — the test only
    asserts that the example script runs to completion without raising.
    """
    try:
        original_cwd = os.getcwd()
        os.chdir(script_path.parent)
        try:
            runpy.run_path(str(script_path), run_name="__main__")
        finally:
            os.chdir(original_cwd)
    except Exception as exc:  # catch and return the exception to the parent
        queue.put((False, repr(exc)))
        return

    queue.put((True, "ok"))


@pytest.mark.parametrize("script_path", collect_example_scripts(), ids=lambda p: p.name)
def test_example_scripts_run_as_main(script_path: Path):
    """Run each `_example*.py` as a script (no validation logic).

    The test runs each example in an isolated process with a 120s timeout and
    fails if the script raises or times out.
    """
    queue: mp.Queue = mp.Queue()
    proc = mp.Process(target=_run_script_as_main, args=(script_path, queue))
    proc.start()
    proc.join(timeout=120)

    if proc.is_alive():
        proc.terminate()
        proc.join()
        pytest.fail("Script timed out after 120s")

    try:
        success, msg = queue.get_nowait()
    except Exception:
        pytest.fail("No result returned from subprocess")

    if not success:
        pytest.fail(msg)
