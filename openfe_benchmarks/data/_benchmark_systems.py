"""
Industry Benchmark Systems for OpenFE

This module provides access to remediated benchmark system inputs ready for use 
with the OpenFE toolkit.
"""

from pathlib import Path
from dataclasses import dataclass
import yaml

from loguru import logger

__all__ = [
    'BenchmarkData',
    'BenchmarkIndex',
    'get_benchmark_data_system',
    'list_benchmark_sets',
    'list_data_systems',
    'PARTIAL_CHARGE_TYPES',
    'get_benchmark_set_data_systems',
    'get_system_index',
]

# Supported partial charge types
PARTIAL_CHARGE_TYPES = ["antechamber_am1bcc", "nagl_openff-gnn-am1bcc-1.0.0.pt", "openeye_am1bcc", "openeye_am1bccelf10"]

_BASE_DIR = Path(__file__).resolve().parent / "benchmark_systems"
_INDEX_FILE = Path(__file__).resolve().parent / "benchmark_system_indexing.yml"


class BenchmarkIndex:
    """
    Singleton class that manages the benchmark system index.
    
    Loads and caches the benchmark system index from YAML on first access.
    Provides methods to query systems by calculation type and validate
    the index against the filesystem.
    """
    _instance = None
    _data = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    def __init__(self):
        # Only load once
        if self._data is None:
            self.reload()
    
    def reload(self):
        """Reload the index from the YAML file."""
        try:
            with open(_INDEX_FILE, 'r') as f:
                self._data = yaml.safe_load(f)
            logger.debug(f"Loaded benchmark index from {_INDEX_FILE}")
        except FileNotFoundError:
            raise ValueError(
                f"Benchmark system index file not found at {_INDEX_FILE}"
            )
        except Exception as e:
            raise ValueError(
                f"Error reading benchmark system index file: {e}"
            )
        
        if 'calculation_types' not in self._data:
            raise ValueError(
                "Invalid index file format: missing 'calculation_types' key"
            )
    
    def get_all_systems(self) -> dict[str, list[str]]:
        """
        Get all benchmark systems from all calculation types.
        
        Returns
        -------
        dict[str, list[str]]
            Dictionary mapping benchmark set names to lists of system names.
        """
        all_systems = {}
        
        for calc_type in self._data['calculation_types'].values():
            if 'benchmark_sets' in calc_type and calc_type['benchmark_sets']:
                for set_name, systems in calc_type['benchmark_sets'].items():
                    if set_name in all_systems:
                        # Merge and deduplicate
                        all_systems[set_name] = sorted(list(set(all_systems[set_name] + systems)))
                    else:
                        all_systems[set_name] = systems
        
        return all_systems
    
    def get_systems_for_calculation(self, calculation_type: str) -> dict[str, list[str]]:
        """
        Get all systems that can be used for a specific calculation type.
        
        Applies logic for system reuse:
        - RBFE systems can be used for any calculation type
        - RSFE and ABFE systems can be used for ASFE
        
        Parameters
        ----------
        calculation_type : str
            The type of calculation: 'rbfe', 'abfe', 'asfe', or 'rsfe'.
        
        Returns
        -------
        dict[str, list[str]]
            Dictionary mapping benchmark set names to lists of system names.
        
        Raises
        ------
        ValueError
            If the calculation_type is not valid.
        """
        valid_types = ['rbfe', 'abfe', 'asfe', 'rsfe']
        
        if calculation_type not in valid_types:
            raise ValueError(
                f"Invalid calculation type '{calculation_type}'. "
                f"Valid types are: {valid_types}"
            )
        
        calc_types = self._data['calculation_types']
        
        # Helper function to merge systems
        def _merge_systems(target_dict: dict, source_dict: dict):
            for set_name, systems in source_dict.items():
                if set_name in target_dict:
                    target_dict[set_name] = sorted(list(set(target_dict[set_name] + systems)))
                else:
                    target_dict[set_name] = systems
        
        # Get systems for the requested calculation type
        result = {}
        
        # Add systems specifically categorized for this calculation type
        if calculation_type in calc_types:
            calc_data = calc_types[calculation_type]
            if 'benchmark_sets' in calc_data and calc_data['benchmark_sets']:
                result = dict(calc_data['benchmark_sets'])
        
        # RBFE systems can be used for any calculation type
        if calculation_type != 'rbfe' and 'rbfe' in calc_types:
            rbfe_data = calc_types['rbfe']
            if 'benchmark_sets' in rbfe_data and rbfe_data['benchmark_sets']:
                _merge_systems(result, rbfe_data['benchmark_sets'])
        
        # For ASFE, also include RSFE and ABFE systems
        if calculation_type == 'asfe':
            if 'rsfe' in calc_types:
                rsfe_data = calc_types['rsfe']
                if 'benchmark_sets' in rsfe_data and rsfe_data['benchmark_sets']:
                    _merge_systems(result, rsfe_data['benchmark_sets'])
            
            if 'abfe' in calc_types:
                abfe_data = calc_types['abfe']
                if 'benchmark_sets' in abfe_data and abfe_data['benchmark_sets']:
                    _merge_systems(result, abfe_data['benchmark_sets'])
        
        return result

# Module-level instance
_benchmark_index = BenchmarkIndex()


def get_system_index(calculation_type: str) -> dict[str, list[str]]:
    """
    Get all benchmark systems that can be used for a specific calculation type.
    
    Systems categorized as 'rbfe' can be used for any calculation type since they
    have all required inputs (protein, ligands, network). Additionally:
    - RSFE and ABFE systems can be used for ASFE calculations
    - Other calculation types return their specifically categorized systems plus
      any systems that can be reused.
    
    Parameters
    ----------
    calculation_type : str
        The type of calculation: 'rbfe', 'abfe', 'asfe', or 'rsfe'.
    
    Returns
    -------
    dict[str, list[str]]
        Dictionary mapping benchmark set names to lists of system names that can
        be used for the specified calculation type.
    
    Raises
    ------
    ValueError
        If the calculation_type is not valid or if the index file cannot be read.
    
    Examples
    --------
    >>> systems = get_system_index('rbfe')
    >>> print(systems['jacs_set'])
    ['bace', 'cdk2', 'jnk1', 'mcl1', 'p38', 'ptp1b', 'thrombin', 'tyk2']
    
    >>> systems = get_system_index('asfe')
    >>> # Returns rbfe, abfe, and rsfe systems since they can be used for asfe
    """
    return _benchmark_index.get_systems_for_calculation(calculation_type)
    

@dataclass
class BenchmarkData:
    """
    Represents a benchmark data system with protein, ligands, optional cofactors, 
    and network.

    Attributes
    ----------
    name : str
        Name of the benchmark data system
    benchmark_set : str
        Fully qualified name of the benchmark set this data system belongs to
        (e.g., 'charge_annihilation_set')
    protein : Path | None
        Path to the protein PDB file (optional, can be None)
    ligands : dict[str, Path]
        Dictionary mapping charge type to ligand SDF file path.
        Always includes 'no_charges' key for the base ligands.sdf file.
        May include charge type keys from PARTIAL_CHARGE_TYPES for charged versions.
    cofactors : dict[str, Path]
        Dictionary mapping charge type to cofactor SDF file path (optional).
        May include 'no_charges' key for the base cofactors.sdf file.
        May include charge type keys from PARTIAL_CHARGE_TYPES for charged versions.
    network : Path
        Paths to network file '*network.json'
    """
    name: str
    benchmark_set: str
    protein: Path | None  # Updated typing to allow None
    ligands: dict[str, Path]
    cofactors: dict[str, Path]
    network: Path
    
    def __repr__(self):
        return (f"BenchmarkData(name='{self.name}', "
                f"benchmark_set='{self.benchmark_set}', "
                f"protein={self.protein.name if self.protein else 'None'}, "
                f"ligands={list(self.ligands.keys())}, "
                f"cofactors={list(self.cofactors.keys())}, "
                f"network={self.network}")


def _discover_benchmark_sets() -> dict[str, list[str]]:
    """
    Discover all available benchmark sets and their data systems from the index.
    
    Uses the YAML index file instead of filesystem traversal for better
    performance. Systems are loaded from benchmark_system_indexing.yml.
    
    Returns
    -------
    dict[str, list[str]]
        Dictionary mapping fully qualified benchmark set names to lists of data system names.
    """
    return _benchmark_index.get_all_systems()


def _validate_and_load_data_system(system_path: Path, system_name: str, 
                               benchmark_set: str) -> BenchmarkData:
    """
    Validate and load a benchmark system from a directory.

    Parameters
    ----------
    system_path : Path
        Path to the system directory.
    system_name : str
        Name of the system.
    benchmark_set : str
        Name of the benchmark set.

    Returns
    -------
    BenchmarkData
        BenchmarkData object.

    Raises
    ------
    ValueError
        If required files are missing or improperly named.
    """
    protein_path = None
    ligands = {}
    cofactors = {}
    network = None
    
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
        
        # Check for network file (network.json)
        if filename.endswith('network.json'):
            if network is not None:
                raise ValueError(f"Multiple network files have been detected. The following has already been saved: {network}")
            network = file_path
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

        logger.error(
            f"Uncategorized file '{filename}' found in system '{system_name}' "
            f"in benchmark set '{benchmark_set}'."
        )
    
    # Validate required files
    if not ligands:
        raise ValueError(
            f"No ligand files found in system '{system_name}' "
            f"in benchmark set '{benchmark_set}'. Expected files named "
            f"'ligands_<charge_type>.sdf' where <charge_type> is one of "
            f"{PARTIAL_CHARGE_TYPES}."
        )

    if 'no_charges' not in ligands:
        raise ValueError(
            f"Missing required 'ligands.sdf' file in system '{system_name}' "
            f"in benchmark set '{benchmark_set}'."
        )
        
    if network is None:
        raise ValueError(f"Missing '*network.json' file in system '{system_name}' "
                         f"in benchmark set '{benchmark_set}'.")
    
    logger.info(
        f"Loaded system '{system_name}' from benchmark set '{benchmark_set}' "
        f"with {len(ligands)} ligand file(s), and {len(cofactors)} cofactor file(s)."
        f" Found protein file: {protein_path is not None}."
    )
    
    return BenchmarkData(
        name=system_name,
        benchmark_set=benchmark_set,
        protein=protein_path,
        ligands=ligands,
        cofactors=cofactors,
        network=network
    )


def get_benchmark_data_system(benchmark_set: str, system_name: str) -> BenchmarkData:
    """
    Factory method to retrieve a benchmark system from a given benchmark set.

    Parameters
    ----------
    benchmark_set : str
        Fully qualified name of the benchmark set (e.g., 'charge_annihilation_set').
    system_name : str
        Name of the system within the benchmark set (e.g., 'cdk2', 'tyk2').

    Returns
    -------
    BenchmarkData
        BenchmarkData object with paths to all relevant files.

    Raises
    ------
    ValueError
        If the benchmark set or system does not exist, or if files are improperly formatted.

    Examples
    --------
    >>> system = get_benchmark_data_system('jacs_set', 'p38')
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
    
    return _validate_and_load_data_system(system_path, system_name, benchmark_set)


def list_benchmark_sets() -> list[str]:
    """
    List all available benchmark sets.

    Returns
    -------
    list[str]
        List of benchmark set names.
    """
    return sorted(_discover_benchmark_sets().keys())


def list_data_systems(benchmark_set: str) -> list[str]:
    """
    List all systems in a given benchmark set.

    Parameters
    ----------
    benchmark_set : str
        Name of the benchmark set.

    Returns
    -------
    list[str]
        List of system names in the benchmark set.

    Raises
    ------
    ValueError
        If the benchmark set does not exist.
    """
    available_sets = _discover_benchmark_sets()
    
    if benchmark_set not in available_sets:
        raise ValueError(
            f"Benchmark set '{benchmark_set}' not found. "
            f"Available benchmark sets: {sorted(available_sets.keys())}"
        )
    
    return available_sets[benchmark_set]


def get_benchmark_set_data_systems(benchmark_set: str) -> dict[str, BenchmarkData]:
    """
    Return all systems in a given benchmark set.

    Parameters
    ----------
    benchmark_set : str
        Name of the benchmark set.

    Returns
    -------
    dict[str, BenchmarkData]
        Dictionary of system names and BenchmarkData objects.

    Raises
    ------
    ValueError
        If the benchmark set does not exist.
    """
    available_sets = _discover_benchmark_sets()
    
    if benchmark_set not in available_sets:
        raise ValueError(
            f"Benchmark set '{benchmark_set}' not found. "
            f"Available benchmark sets: {sorted(available_sets.keys())}"
        )
        
    benchmark_systems = {}
    for system_name in available_sets[benchmark_set]:
        system_path = _BASE_DIR / benchmark_set.replace('.', '/') / system_name
        logger.debug(f"Loading benchmark system '{system_name}' from '{benchmark_set}'...")
        benchmark_systems[system_name] = _validate_and_load_data_system(system_path, system_name, benchmark_set)
    
    return benchmark_systems