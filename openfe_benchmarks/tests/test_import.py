"""
Import tests
"""

import sys
from pathlib import Path


def test_benchmark_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "openfe_benchmarks" in sys.modules


def test_yaml_files_contain_no_todos():
    """Ensure no YAML files in the repository still contain TODO markers."""
    repo_root = Path(__file__).resolve().parents[2]
    yaml_files = list(repo_root.rglob("*.yaml"))
    todos: list[str] = []
    for path in yaml_files:
        text = path.read_text(encoding="utf-8", errors="ignore")
        for line_number, line in enumerate(text.splitlines(), start=1):
            if "TODO" in line:
                todos.append(
                    f"{path.relative_to(repo_root)}:{line_number}: {line.strip()}"
                )

    assert not todos, "Found TODO markers in YAML files:\n" + "\n".join(todos)
