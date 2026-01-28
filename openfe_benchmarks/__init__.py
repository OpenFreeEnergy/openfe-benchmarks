"""OpenFE Benchmarks: Benchmarking systems for OpenFE."""

from importlib.metadata import version
from loguru import logger
import sys

__all__ = ()

__version__ = version("openfe_benchmarks")

# Configure logger for the package
# Remove default handler and add a single configured handler
logger.remove()
logger.add(
    sys.stderr,
    format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",
    level="INFO"
)