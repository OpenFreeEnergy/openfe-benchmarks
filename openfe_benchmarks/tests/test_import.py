"""
Import tests
"""

import sys


def test_benchmark_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "openfe_benchmarks" in sys.modules
