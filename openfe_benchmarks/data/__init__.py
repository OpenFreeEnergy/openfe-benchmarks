"""
Industry Benchmark Systems for OpenFE

This module provides access to remediated benchmark system inputs ready for use 
with the OpenFE toolkit.
"""

from pathlib import Path
from dataclasses import dataclass
from loguru import logger


# Supported partial charge types
PARTIAL_CHARGE_TYPES = ["antechamber_am1bcc", "nagl_openff-gnn-am1bcc-1.0.0.pt", "openeye_am1bcc", "openeye_am1bccelf10"]

# Base directory for industry benchmark systems
_BASE_DIR = Path(__file__).parent


@dataclass
class BenchmarkSystem:
    """
    Represents a benchmark system with protein, ligands, optional cofactors, 
    and networks.
    
    Attributes:
        name: Name of the benchmark system
        benchmark_set: Fully qualified name of the benchmark set this system belongs to
                      (e.g., 'industry_benchmark_systems.charge_annihilation_set')
        protein: Path to the protein PDB file
        ligands: Dictionary mapping charge type to ligand SDF file path.
                Always includes 'no_charges' key for the base ligands.sdf file.
                May include charge type keys from PARTIAL_CHARGE_TYPES for charged versions.
        cofactors: Dictionary mapping charge type to cofactor SDF file path (optional).
                  May include 'no_charges' key for the base cofactors.sdf file.
                  May include charge type keys from PARTIAL_CHARGE_TYPES for charged versions.
        networks: list of paths to network files (e.g., *network.json files)
    """
    name: str
    benchmark_set: str
    protein: Path
    ligands: dict[str, Path]
    cofactors: dict[str, Path]
    networks: list[Path]
    
    def __repr__(self):
        return (f"BenchmarkSystem(name='{self.name}', "
                f"benchmark_set='{self.benchmark_set}', "
                f"protein={self.protein.name}, "
                f"ligands={list(self.ligands.keys())}, "
                f"cofactors={list(self.cofactors.keys())}, "
                f"networks={[n.name for n in self.networks]})")


def _discover_benchmark_sets() -> dict[str, list[str]]:
    """
    Discover all available benchmark sets and their systems.
    
    A benchmark system is identified by the presence of a 'PREPARATION_DETAILS.md' file
    (case-insensitive) in its directory. Benchmark sets are represented as hierarchical
    dot-separated paths (e.g., 'industry_benchmark_systems.charge_annihilation_set').
    
    Returns:
        Dictionary mapping fully qualified benchmark set names to lists of system names
    """
    def _get_qualified_name(path: Path) -> str:
        """Convert a path to a dot-separated qualified name."""
        return ".".join(path.relative_to(_BASE_DIR).parts)
    
    def _has_preparation_details(directory: Path) -> bool:
        """Check if directory contains PREPARATION_DETAILS.md (case-insensitive)."""
        for file in directory.glob('*'):
            if file.is_file() and file.name.lower() == 'preparation_details.md':
                return True
        return False
    
    benchmark_sets = {}
    
    # Recursively traverse directory structure to find benchmark systems
    def _traverse(current_path: Path, depth: int = 0, max_depth: int = 5):
        """Recursively traverse directories to find benchmark systems."""
        if depth > max_depth:
            return
        
        for item in current_path.iterdir():
            if not item.is_dir() or item.name.startswith('_') or item.name.startswith('.'):
                continue
            
            # Check if this directory contains systems with PREPARATION_DETAILS.md
            systems = []
            for potential_system in item.iterdir():
                if potential_system.is_dir() and _has_preparation_details(potential_system):
                    systems.append(potential_system.name)
            
            if systems:
                # This is a benchmark set
                set_name = _get_qualified_name(item)
                benchmark_sets[set_name] = sorted(systems)
            else:
                # Continue traversing deeper
                _traverse(item, depth + 1, max_depth)
    
    _traverse(_BASE_DIR)
    
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
    networks = []
    
    # Track all files for validation
    all_files = list(system_path.glob('*'))
    categorized_files = set()
    
    for file_path in all_files:
        if not file_path.is_file():
            continue
            
        filename = file_path.name
        
        # Skip PREPARATION_DETAILS.md (case-insensitive)
        if filename.lower() == 'preparation_details.md':
            categorized_files.add(file_path)
            continue
        
        # Check for protein PDB
        if filename == 'protein.pdb':
            protein_path = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found protein: {filename}")
            continue
        
        # Check for base ligands file (without charges)
        if filename == 'ligands.sdf':
            ligands['no_charges'] = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found ligands without charges: {filename}")
            continue
        
        # Check for ligand files with specific charge types
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
        
        # Check for base cofactors file (without charges)
        if filename == 'cofactors.sdf':
            cofactors['no_charges'] = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found cofactors without charges: {filename}")
            continue
        
        # Check for cofactor files with specific charge types
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
        
        # Check for network files (*network.json)
        if filename.endswith('network.json'):
            networks.append(file_path)
            categorized_files.add(file_path)
            logger.debug(f"Found network: {filename}")
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
    
    if 'no_charges' not in ligands:
        raise ValueError(
            f"Missing required 'ligands.sdf' file in system '{system_name}' "
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
        f"and {len(networks)} network file(s)."
    )
    
    return BenchmarkSystem(
        name=system_name,
        benchmark_set=benchmark_set,
        protein=protein_path,
        ligands=ligands,
        cofactors=cofactors,
        networks=networks
    )


def get_benchmark_system(benchmark_set: str, system_name: str) -> BenchmarkSystem:
    """
    Factory method to retrieve a benchmark system from a given benchmark set.
    
    Args:
        benchmark_set: Fully qualified name of the benchmark set 
                      (e.g., 'industry_benchmark_systems.charge_annihilation_set')
        system_name: Name of the system within the benchmark set (e.g., 'cdk2', 'tyk2')
        
    Returns:
        BenchmarkSystem object with paths to all relevant files
        
    Raises:
        ValueError: If the benchmark set or system does not exist, or if files 
                   are improperly formatted
        
    Examples:
        >>> system = get_benchmark_system('industry_benchmark_systems.jacs_set', 'p38')
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
    
    # Load and validate the system - convert dot-separated path to filesystem path
    system_path = _BASE_DIR / benchmark_set.replace('.', '/') / system_name
    
    logger.debug(f"Loading benchmark system '{system_name}' from '{benchmark_set}'...")
    
    return _validate_and_load_system(system_path, system_name, benchmark_set)


def list_benchmark_sets() -> list[str]:
    """
    list all available benchmark sets.
    
    Returns:
        list of benchmark set names
    """
    return sorted(_discover_benchmark_sets().keys())


def list_systems(benchmark_set: str) -> list[str]:
    """
    list all systems in a given benchmark set.
    
    Args:
        benchmark_set: Name of the benchmark set
        
    Returns:
        list of system names in the benchmark set
        
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