"""
Industry Benchmark Systems for OpenFE

This module provides access to remediated benchmark system inputs ready for use 
with the OpenFE toolkit.
"""

from pathlib import Path
from typing import Optional, List
from dataclasses import dataclass
from loguru import logger


# Supported partial charge types
PARTIAL_CHARGE_TYPES = ["antechamber_am1bcc", "openeye_elf10"]

# Base directory for industry benchmark systems
_BASE_DIR = Path(__file__).parent


@dataclass
class BenchmarkSystem:
    """
    Represents a benchmark system with protein, ligands, optional cofactors, 
    and network mappings.
    
    Attributes:
        name: Name of the benchmark system
        benchmark_set: Name of the benchmark set this system belongs to
        protein: Path to the protein PDB file
        ligands: Dictionary mapping charge type to ligand SDF file path
        cofactors: Dictionary mapping charge type to cofactor SDF file path (optional)
        mappings: List of paths to network mapping files (e.g., .graphml files)
    """
    name: str
    benchmark_set: str
    protein: Path
    ligands: dict[str, Path]
    cofactors: dict[str, Path]
    mappings: List[Path]
    
    def __repr__(self):
        return (f"BenchmarkSystem(name='{self.name}', "
                f"benchmark_set='{self.benchmark_set}', "
                f"protein={self.protein.name}, "
                f"ligands={list(self.ligands.keys())}, "
                f"cofactors={list(self.cofactors.keys())}, "
                f"mappings={[m.name for m in self.mappings]})")


def _discover_benchmark_sets() -> dict[str, List[str]]:
    """
    Discover all available benchmark sets and their systems.
    
    Returns:
        Dictionary mapping benchmark set names to lists of system names
    """
    benchmark_sets = {}
    
    for item in _BASE_DIR.iterdir():
        if item.is_dir() and not item.name.startswith('_') and not item.name.startswith('.'):
            systems = []
            for subitem in item.iterdir():
                if subitem.is_dir() and not subitem.name.startswith('_') and not subitem.name.startswith('.'):
                    systems.append(subitem.name)
            if systems:
                benchmark_sets[item.name] = sorted(systems)
    
    return benchmark_sets


def _validate_and_load_system(system_path: Path, system_name: str, 
                               benchmark_set: str) -> BenchmarkSystem:
    """
    Validate and load a benchmark system from a directory.
    
    Args:
        system_path: Path to the system directory
        system_name: Name of the system
        benchmark_set: Name of the benchmark set
        
    Returns:
        BenchmarkSystem object
        
    Raises:
        ValueError: If required files are missing or improperly named
    """
    protein_path = None
    ligands = {}
    cofactors = {}
    mappings = []
    
    # Track all files for validation
    all_files = list(system_path.glob('*'))
    categorized_files = set()
    
    for file_path in all_files:
        if not file_path.is_file():
            continue
            
        filename = file_path.name
        
        # Skip PREPARATION_DETAILS.md
        if filename == 'PREPARATION_DETAILS.md':
            categorized_files.add(file_path)
            continue
        
        # Check for protein PDB
        if filename == 'protein.pdb':
            protein_path = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found protein: {filename}")
            continue
        
        # Check for ligand files
        if filename.startswith('ligands_') and filename.endswith('.sdf'):
            charge_type = filename[len('ligands_'):-len('.sdf')]
            if charge_type not in PARTIAL_CHARGE_TYPES:
                raise ValueError(
                    f"Unsupported partial charge type '{charge_type}' in file '{filename}' "
                    f"for system '{system_name}' in benchmark set '{benchmark_set}'. "
                    f"Supported types: {PARTIAL_CHARGE_TYPES}. "
                    f"If this is a new charge type, please contribute it by adding to "
                    f"PARTIAL_CHARGE_TYPES in the module."
                )
            ligands[charge_type] = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found ligands with {charge_type} charges: {filename}")
            continue
        
        # Check for cofactor files
        if filename.startswith('cofactors_') and filename.endswith('.sdf'):
            charge_type = filename[len('cofactors_'):-len('.sdf')]
            if charge_type not in PARTIAL_CHARGE_TYPES:
                raise ValueError(
                    f"Unsupported partial charge type '{charge_type}' in file '{filename}' "
                    f"for system '{system_name}' in benchmark set '{benchmark_set}'. "
                    f"Supported types: {PARTIAL_CHARGE_TYPES}. "
                    f"If this is a new charge type, please contribute it by adding to "
                    f"PARTIAL_CHARGE_TYPES in the module."
                )
            cofactors[charge_type] = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found cofactors with {charge_type} charges: {filename}")
            continue
        
        # Check for network mapping files (graphml, etc.)
        if filename.endswith('.graphml') or 'network' in filename.lower():
            mappings.append(file_path)
            categorized_files.add(file_path)
            logger.debug(f"Found network mapping: {filename}")
            continue
    
    # Check for uncategorized files
    uncategorized = set(all_files) - categorized_files
    for file_path in uncategorized:
        if not file_path.is_file():
            continue
        
        filename = file_path.name
        
        # Error on uncategorized PDB or SDF files
        if filename.endswith('.pdb'):
            raise ValueError(
                f"Uncategorized PDB file '{filename}' found in system '{system_name}' "
                f"in benchmark set '{benchmark_set}'. Expected 'protein.pdb'."
            )
        
        if filename.endswith('.sdf'):
            raise ValueError(
                f"Uncategorized SDF file '{filename}' found in system '{system_name}' "
                f"in benchmark set '{benchmark_set}'. Expected format: "
                f"'ligands_<charge_type>.sdf' or 'cofactors_<charge_type>.sdf' "
                f"where <charge_type> is one of {PARTIAL_CHARGE_TYPES}."
            )
        
        # Warning for other uncategorized files
        logger.warning(
            f"Uncategorized file '{filename}' found in system '{system_name}' "
            f"in benchmark set '{benchmark_set}'. This file will be ignored."
        )
    
    # Validate required files
    if protein_path is None:
        raise ValueError(
            f"Missing required 'protein.pdb' file in system '{system_name}' "
            f"in benchmark set '{benchmark_set}'."
        )
    
    if not ligands:
        raise ValueError(
            f"No ligand files found in system '{system_name}' "
            f"in benchmark set '{benchmark_set}'. Expected files named "
            f"'ligands_<charge_type>.sdf' where <charge_type> is one of "
            f"{PARTIAL_CHARGE_TYPES}."
        )
    
    logger.info(
        f"Loaded system '{system_name}' from benchmark set '{benchmark_set}' "
        f"with {len(ligands)} ligand file(s), {len(cofactors)} cofactor file(s), "
        f"and {len(mappings)} mapping file(s)."
    )
    
    return BenchmarkSystem(
        name=system_name,
        benchmark_set=benchmark_set,
        protein=protein_path,
        ligands=ligands,
        cofactors=cofactors,
        mappings=mappings
    )


def get_benchmark_system(benchmark_set: str, system_name: str) -> BenchmarkSystem:
    """
    Factory method to retrieve a benchmark system from a given benchmark set.
    
    Args:
        benchmark_set: Name of the benchmark set (e.g., 'jacs_set', 'fragments')
        system_name: Name of the system within the benchmark set (e.g., 'p38', 'tyk2')
        
    Returns:
        BenchmarkSystem object with paths to all relevant files
        
    Raises:
        ValueError: If the benchmark set or system does not exist, or if files 
                   are improperly formatted
        
    Examples:
        >>> system = get_benchmark_system('jacs_set', 'p38')
        >>> print(system.protein)
        >>> print(system.ligands['antechamber_am1bcc'])
    """
    # Discover available benchmark sets and systems
    available_sets = _discover_benchmark_sets()
    
    # Check if benchmark set exists
    if benchmark_set not in available_sets:
        raise ValueError(
            f"Benchmark set '{benchmark_set}' not found. "
            f"Available benchmark sets: {sorted(available_sets.keys())}"
        )
    
    # Check if system exists in the benchmark set
    if system_name not in available_sets[benchmark_set]:
        raise ValueError(
            f"System '{system_name}' not found in benchmark set '{benchmark_set}'. "
            f"Available systems in '{benchmark_set}': {available_sets[benchmark_set]}"
        )
    
    # Load and validate the system
    system_path = _BASE_DIR / benchmark_set / system_name
    
    logger.debug(f"Loading benchmark system '{system_name}' from '{benchmark_set}'...")
    
    return _validate_and_load_system(system_path, system_name, benchmark_set)


def list_benchmark_sets() -> List[str]:
    """
    List all available benchmark sets.
    
    Returns:
        List of benchmark set names
    """
    return sorted(_discover_benchmark_sets().keys())


def list_systems(benchmark_set: str) -> List[str]:
    """
    List all systems in a given benchmark set.
    
    Args:
        benchmark_set: Name of the benchmark set
        
    Returns:
        List of system names in the benchmark set
        
    Raises:
        ValueError: If the benchmark set does not exist
    """
    available_sets = _discover_benchmark_sets()
    
    if benchmark_set not in available_sets:
        raise ValueError(
            f"Benchmark set '{benchmark_set}' not found. "
            f"Available benchmark sets: {sorted(available_sets.keys())}"
        )
    
    return available_sets[benchmark_set]


__all__ = [
    'BenchmarkSystem',
    'get_benchmark_system',
    'list_benchmark_sets',
    'list_systems',
    'PARTIAL_CHARGE_TYPES',
]
