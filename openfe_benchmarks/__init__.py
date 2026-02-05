"""OpenFE Benchmarks: Free Energy Benchmark Data Repository"""

from importlib.metadata import version
import logging
import sys

__all__ = ()

__version__ = version("openfe_benchmarks")

# Configure logging for the package
# Only configure if no handlers are already set up (to avoid duplicate output)
_logger = logging.getLogger(__name__)
if not _logger.handlers:
    _handler = logging.StreamHandler(sys.stderr)
    _formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    _handler.setFormatter(_formatter)
    _logger.addHandler(_handler)
    _logger.setLevel(logging.INFO)
