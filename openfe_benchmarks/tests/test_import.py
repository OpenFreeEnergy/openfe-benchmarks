"""
Import tests
"""

import openfe_benchmarks
import pytest
import sys


def test_arsenic_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "openfe_benchmarks" in sys.modules
