"""This module provides access to remediated benchmark system inputs ready for use
with the OpenFE toolkit.
"""

from pathlib import Path
from dataclasses import dataclass
import yaml
import logging

# Get logger for this module - will inherit configuration from parent
logger = logging.getLogger(__name__)

__all__ = [
    "BenchmarkData",
    "BenchmarkIndex",
    "get_data_by_system_name",
    "get_data_by_system_group",
    "PARTIAL_CHARGE_TYPES",
]

# Supported partial charge types
PARTIAL_CHARGE_TYPES = [
    "antechamber_am1bcc",
    "nagl_openff-gnn-am1bcc-1.0.0.pt",
    "openeye_am1bcc",
    "openeye_am1bccelf10",
]

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
            with open(_INDEX_FILE, "r") as f:
                raw_data = yaml.safe_load(f) or {}
            logger.debug(f"Loaded benchmark index from {_INDEX_FILE}")
        except FileNotFoundError:
            raise ValueError(f"Benchmark system index file not found at {_INDEX_FILE}")
        except Exception as e:
            raise ValueError(f"Error reading benchmark system index file: {e}")

        # Adapt the structure to ensure 'systems' key exists
        systems_data = {
            key: value for key, value in raw_data.items() if key != "metadata"
        }

        if not systems_data:
            raise ValueError("Invalid index file format: no benchmark systems found.")

        self._data = {"systems": systems_data}

        logger.debug("Benchmark index successfully reloaded and validated.")

    def list_system_names_by_tag(self, tags: list[str] = []) -> list[tuple[str, str]]:
        """
        Get all systems that match **all** of the provided tags.

        Parameters
        ----------
        tags : list[str]
            List of tags to filter by (e.g., ['protein', 'cofactors']).
            Only systems containing every tag in this list will be returned.
            Defaults to ``[]``.

        Returns
        -------
        list[tuple[str, str]]
            List of tuples containing (system_group, system_name).

        Examples
        --------
        >>> index = BenchmarkIndex()
        >>> systems = index.get_systems_by_tag(['protein'])
        >>> # Returns all systems tagged with 'protein'
        >>> systems_with_cofactors = index.get_systems_by_tag(['cofactors'])
        >>> # Returns all systems tagged with 'cofactors'
        >>> protein_and_cofactors = index.get_systems_by_tag(['protein', 'cofactors'])
        >>> # Returns systems tagged with both 'protein' and 'cofactors'
        """
        if not self._data or not self._data.get("systems"):
            logger.error("Benchmark index data is not loaded or is invalid.")
            return []

        matching_systems = []

        for system_group, systems in self._data["systems"].items():
            for system_name, system_data in systems.items():
                # Check if all of the requested tags match
                if all(tag in system_data for tag in tags):
                    matching_systems.append((system_group, system_name))

        return matching_systems

    def list_available_tags(self) -> set[str]:
        """
        Get all unique tags used across all systems.

        Returns
        -------
        set[str]
            Set of all available tags.
        """
        if not self._data or not self._data.get("systems"):
            logger.error("Benchmark index data is not loaded or is invalid.")
            return set()

        tags = set()
        for systems in self._data["systems"].values():
            for system_tags in systems.values():
                tags.update(system_tags)

        return tags

    def list_system_names_by_group(self, system_group: str) -> list[str]:
        """
        Get all system names in a specific system group.

        Parameters
        ----------
        system_group : str
            The system group name.

        Returns
        -------
        list[str]
            List of system names in the system group.

        Raises
        ------
        ValueError
            If the system group doesn't exist.
        """
        if not self._data or not self._data.get("systems"):
            logger.error("Benchmark index data is not loaded or is invalid.")
            return []

        if system_group not in self._data["systems"]:
            available = list(self._data["systems"].keys())
            raise ValueError(
                f"System group '{system_group}' not found. Available groups: {available}"
            )

        return list(self._data["systems"][system_group].keys())

    def list_system_groups(self) -> list[str]:
        """
        List all available system groups.

        Returns
        -------
        list[str]
            Sorted list of system group names.
        """
        if not self._data or not self._data.get("systems"):
            logger.error("Benchmark index data is not loaded or is invalid.")
            return []

        return sorted(self._data["systems"].keys())


# Module-level instance
_benchmark_index = BenchmarkIndex()


@dataclass
class BenchmarkData:
    """
    Represents a benchmark data system with protein, ligands, optional cofactors,
    and ligand networks.

    Attributes
    ----------
    system_name : str
        Name of the benchmark data system
    system_group : str
        Fully qualified name of the system group this data system belongs to
        (e.g., 'charge_annihilation_set')
    protein : Path | None
        Path to the protein PDB file (optional, can be None)
    ligands : dict[str, Path]
        Dictionary mapping charge type to ligand SDF file path.
        Always includes 'no_charges' key for the base ligands.sdf file.
        May include charge type keys from PARTIAL_CHARGE_TYPES for charged versions.
    cofactors : dict[str, Path] | None
        Dictionary mapping charge type to cofactor SDF file path (optional).
        May include 'no_charges' key for the base cofactors.sdf file.
        May include charge type keys from PARTIAL_CHARGE_TYPES for charged versions.
    ligand_networks : dict[str, Path] | None
        Dictionary of available ligand networks where the key is the filename, '*network.json',
        and the value is a Path to a ligand network file.
    reference_data : dict[str, Path] | None
        Dictionary of available reference data files where the key is the filename 'experimental*data.json',
        and the value is a Path to a reference data file.
    subset_data : dict[str, Path] | None
        Dictionary of available subset files where the key is the filename stem 'subset*.json',
        and the value is a Path to a subset JSON file defining a named subset of ligands.
    details : str
        Information available in the preparation_details.md file
    """

    system_name: str
    system_group: str
    protein: Path | None  # Updated typing to allow None
    ligands: dict[str, Path]
    cofactors: dict[str, Path] | None
    ligand_networks: dict[str, Path] | None
    reference_data: dict[str, Path] | None
    subset_data: dict[str, Path] | None
    details: str

    def __repr__(self):
        return (
            f"BenchmarkData(system_name='{self.system_name}', "
            f"system_group='{self.system_group}', "
            f"protein={self.protein.name if self.protein else 'None'}, "
            f"ligands={list(self.ligands.keys())}, "
            f"cofactors={list(self.cofactors.keys()) if self.cofactors is not None else 'None'}, "
            f"ligand_network={list(self.ligand_networks.keys()) if self.ligand_networks is not None else 'None'}, "
            f"reference_data={list(self.reference_data.keys()) if self.reference_data is not None else 'None'}, "
            f"subset_data={list(self.subset_data.keys()) if self.subset_data is not None else 'None'})"
        )


def _validate_and_load_data_system(
    system_path: Path, system_name: str, system_group: str
) -> BenchmarkData:
    """
    Validate and load a benchmark system from a directory.

    Parameters
    ----------
    system_path : Path
        Path to the system directory.
    system_name : str
        Name of the system.
    system_group : str
        Name of the system group.

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
    ligand_networks = {}
    reference_data = {}
    subset_data = {}
    details = None

    # Track all files for validation
    all_files = list(system_path.glob("*"))
    categorized_files = set()

    for file_path in all_files:
        if not file_path.is_file():
            continue

        filename = file_path.name

        # Skip PREPARATION_DETAILS.md (case-insensitive)
        if filename.lower() == "preparation_details.md":
            details = file_path.read_text(encoding="utf-8")
            categorized_files.add(file_path)
            continue

        # Check for protein PDB
        if filename == "protein.pdb":
            protein_path = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found protein: {filename}")
            continue

        # Check for base ligands file (without charges)
        if filename == "ligands.sdf":
            ligands["no_charges"] = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found ligands without charges: {filename}")
            continue

        # Check for ligand files with specific charge types
        if filename.startswith("ligands_") and filename.endswith(".sdf"):
            charge_type = filename[len("ligands_") : -len(".sdf")]
            if charge_type not in PARTIAL_CHARGE_TYPES:
                raise ValueError(
                    f"Unsupported partial charge type '{charge_type}' in file '{filename}' "
                    f"for system '{system_name}' in system group '{system_group}'. "
                    f"Supported types: {PARTIAL_CHARGE_TYPES}. "
                    f"If this is a new charge type, please contribute it by adding to "
                    f"PARTIAL_CHARGE_TYPES in the module."
                )
            ligands[charge_type] = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found ligands with {charge_type} charges: {filename}")
            continue

        # Check for base cofactors file (without charges)
        if filename == "cofactors.sdf":
            cofactors["no_charges"] = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found cofactors without charges: {filename}")
            continue

        # Check for cofactor files with specific charge types
        if filename.startswith("cofactors_") and filename.endswith(".sdf"):
            charge_type = filename[len("cofactors_") : -len(".sdf")]
            if charge_type not in PARTIAL_CHARGE_TYPES:
                raise ValueError(
                    f"Unsupported partial charge type '{charge_type}' in file '{filename}' "
                    f"for system '{system_name}' in system group '{system_group}'. "
                    f"Supported types: {PARTIAL_CHARGE_TYPES}. "
                    f"If this is a new charge type, please contribute it by adding to "
                    f"PARTIAL_CHARGE_TYPES in the module."
                )
            cofactors[charge_type] = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found cofactors with {charge_type} charges: {filename}")
            continue

        # Check for ligand network file (network.json)
        if "network" in filename and filename.endswith(".json"):
            ligand_networks[file_path.stem] = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found ligand network: {filename}")
            continue

        # check for reference data file (experimental*data.json)
        if filename.startswith("experimental") and filename.endswith("data.json"):
            reference_data[file_path.stem] = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found reference data: {filename}")
            continue

        # check for subset data file (subset*.json)
        if filename.startswith("subset") and filename.endswith(".json"):
            subset_data[file_path.stem] = file_path
            categorized_files.add(file_path)
            logger.debug(f"Found subset: {filename}")
            continue

    # Check for uncategorized files
    uncategorized = set(all_files) - categorized_files
    for file_path in uncategorized:
        if not file_path.is_file():
            continue

        filename = file_path.name

        # Error on uncategorized PDB or SDF files
        if filename.endswith(".pdb"):
            raise ValueError(
                f"Uncategorized PDB file '{filename}' found in system '{system_name}' "
                f"in system group '{system_group}'. Expected 'protein.pdb'."
            )

        if filename.endswith(".sdf"):
            raise ValueError(
                f"Uncategorized SDF file '{filename}' found in system '{system_name}' "
                f"in system group '{system_group}'. Expected format: "
                f"'ligands_<charge_type>.sdf' or 'cofactors_<charge_type>.sdf' "
                f"where <charge_type> is one of {PARTIAL_CHARGE_TYPES}."
            )

        if filename.endswith(".json"):
            raise ValueError(
                f"Uncategorized JSON file '{filename}' found in system '{system_name}' "
                f"in system group '{system_group}'. Expected format: "
                f"'*network*.json' for ligand networks, 'experimental*data.json', or 'subset*.json'."
            )

        raise ValueError(
            f"Uncategorized file '{filename}' found in system '{system_name}' "
            f"in system group '{system_group}'."
        )

    # Validate required files
    if not ligands:
        raise ValueError(
            f"No ligand files found in system '{system_name}' "
            f"in system group '{system_group}'. Expected files named "
            f"'ligands_<charge_type>.sdf' where <charge_type> is one of "
            f"{PARTIAL_CHARGE_TYPES}."
        )

    if "no_charges" not in ligands:
        raise ValueError(
            f"Missing required 'ligands.sdf' file in system '{system_name}' "
            f"in system group '{system_group}'."
        )

    if details is None:
        raise ValueError(
            "The file 'preparation_details.md' if missing from"
            f"system '{system_name}' in system group '{system_group}'."
        )

    logger.info(
        f"Loaded system '{system_name}' from system group '{system_group}' with:\n"
        f"  {len(ligands)} ligand file(s), and {len(cofactors)} cofactor file(s).\n"
        f"  Found protein file: {protein_path is not None}.\n"
        f"  Found {len(ligand_networks)} ligand network files with keys: {', '.join(list(ligand_networks.keys()))}"
    )

    return BenchmarkData(
        system_name=system_name,
        system_group=system_group,
        protein=protein_path,
        ligands=ligands,
        cofactors=cofactors or None,
        ligand_networks=ligand_networks or None,
        reference_data=reference_data or None,
        subset_data=subset_data or None,
        details=details,
    )


def get_data_by_system_name(system_group: str, system_name: str) -> BenchmarkData:
    """
    Factory method to retrieve a benchmark system from a given system group.

    Parameters
    ----------
    system_group : str
        Fully qualified name of the system group (e.g., 'charge_annihilation_set').
    system_name : str
        Name of the system within the system group (e.g., 'cdk2', 'tyk2').

    Returns
    -------
    BenchmarkData
        BenchmarkData object with paths to all relevant files.

    Raises
    ------
    ValueError
        If the system group or system does not exist, or if files are improperly formatted.

    Examples
    --------
    >>> system = get_data_by_system_name('jacs_set', 'p38')
    >>> print(system.protein)
    >>> print(system.ligands['antechamber_am1bcc'])
    """
    # Check if system group exists
    available_groups = _benchmark_index.list_system_groups()
    if system_group not in available_groups:
        raise ValueError(
            f"System group '{system_group}' not found. "
            f"Available system groups: {available_groups}"
        )

    # Check if system exists in the system group
    available_systems = _benchmark_index.list_system_names_by_group(system_group)
    if system_name not in available_systems:
        raise ValueError(
            f"System '{system_name}' not found in system group '{system_group}'. "
            f"Available systems in '{system_group}': {available_systems}"
        )

    # Load and validate the system - convert dot-separated path to filesystem path
    system_path = _BASE_DIR / system_group.replace(".", "/") / system_name

    logger.debug(f"Loading benchmark system '{system_name}' from '{system_group}'...")

    return _validate_and_load_data_system(system_path, system_name, system_group)


def get_data_by_system_group(system_group: str) -> dict[str, BenchmarkData]:
    """
    Retrieve all benchmark systems in a given system group.

    Parameters
    ----------
    system_group : str
        The name of the system group.

    Returns
    -------
    dict[str, BenchmarkData]
        A dictionary of BenchmarkData objects for all systems in the system group.

    Raises
    ------
    ValueError
        If the system group does not exist.
    """
    # Check if system group exists
    available_groups = _benchmark_index.list_system_groups()
    if system_group not in available_groups:
        raise ValueError(
            f"System group '{system_group}' not found. "
            f"Available system groups: {available_groups}"
        )

    # Retrieve all systems in the system group
    system_names = _benchmark_index.list_system_names_by_group(system_group)

    # Load and return all systems as BenchmarkData objects
    return {
        system_name: get_data_by_system_name(system_group, system_name)
        for system_name in system_names
    }
