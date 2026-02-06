"""
Benchmark Data Systems for OpenFE

This module provides access to remediated benchmark system inputs ready for use
with the OpenFE toolkit.
"""

from ._benchmark_systems import BenchmarkData as BenchmarkData
from ._benchmark_systems import BenchmarkIndex as BenchmarkIndex
from ._benchmark_systems import get_benchmark_data_system as get_benchmark_data_system
from ._benchmark_systems import (
    get_benchmark_set_data_systems as get_benchmark_set_data_systems,
)
from ._benchmark_systems import PARTIAL_CHARGE_TYPES as PARTIAL_CHARGE_TYPES

from ._benchmark_systems import __all__ as _benchmark_systems_all

__all__ = _benchmark_systems_all

"""
.. include:: README.md
"""
