"""Tests that execute example scripts to ensure they run as `__main__`.

These tests locate any files matching `scripts/_example*.py` and execute
them using `runpy.run_path` in a separate process to avoid state leakage
and to allow a timeout for scripts that might hang.
"""
from pathlib import Path
import runpy
import multiprocessing as mp
import pytest


SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"


def _run_script(path: Path, queue: mp.Queue):
    """Helper to run a script and report exceptions via a queue."""
    try:
        runpy.run_path(str(path), run_name="__main__")
    except Exception as e:
        queue.put((False, repr(e)))
        return
    queue.put((True, "ok"))


def test_example_scripts_run_as_main():
    """Run each `scripts/_example*.py` file as `__main__` in a subprocess.

    Fails if any script raises an exception or times out.
    """
    patterns = list(SCRIPTS_DIR.glob("_example*.py"))
    if not patterns:
        pytest.skip("No example scripts found to test")

    errors = []
    for script in patterns:
        queue: mp.Queue = mp.Queue()
        p = mp.Process(target=_run_script, args=(script, queue))
        p.start()
        p.join(timeout=30)
        if p.is_alive():
            p.terminate()
            errors.append(f"{script}: timeout")
            continue

        try:
            ok, msg = queue.get_nowait()
        except Exception:
            errors.append(f"{script}: no result returned")
            continue

        if not ok:
            errors.append(f"{script}: error: {msg}")

    if errors:
        pytest.fail("\n".join(errors))
