"""
Tests for the industry benchmark systems module.
"""
import pytest
from pathlib import Path

from rdkit import Chem

from openfe import ProteinComponent, SmallMoleculeComponent, LigandNetwork

from openfe_benchmarks.data import (
    get_benchmark_system,
    list_benchmark_sets,
    list_systems,
    PARTIAL_CHARGE_TYPES,
)

def test_list_benchmark_sets_output():
    """Test that list_benchmark_sets() outputs a list of strings."""
    benchmark_sets = list_benchmark_sets()
    assert isinstance(benchmark_sets, list)
    assert all(isinstance(item, str) for item in benchmark_sets)

@pytest.fixture(scope="module")
def benchmark_sets():
    """Fixture to provide all available benchmark sets."""
    return list_benchmark_sets()

def test_list_benchmark_sets_tuple_output(benchmark_sets):
    """Test that list_benchmark_sets() outputs a list of tuples of length 2 containing strings."""
    for benchmark_set in benchmark_sets:
        benchmark_systems = list_systems(benchmark_set)
        assert isinstance(benchmark_systems, list)
        assert all(isinstance(item, str) for item in benchmark_systems)

@pytest.mark.parametrize("benchmark_set, system_name", [
    (benchmark_set, system_name)
    for benchmark_set in list_benchmark_sets()
    for system_name in list_systems(benchmark_set)
])
def test_benchmark_system_initialization(benchmark_set, system_name):
    """Test initialization of all BenchmarkSystem objects."""
    system = get_benchmark_system(benchmark_set, system_name)
    assert system is not None
    assert system.name == system_name
    assert system.benchmark_set == benchmark_set

    missing_ligand_charges = [charge for charge in PARTIAL_CHARGE_TYPES if charge not in system.ligands]
    if missing_ligand_charges:
        print(f"Warning: Missing ligand charge types for {system_name}: {missing_ligand_charges}")

    missing_cofactor_charges = [charge for charge in PARTIAL_CHARGE_TYPES if charge not in system.cofactors]
    if missing_cofactor_charges:
        print(f"Warning: Missing cofactor charge types for {system_name}: {missing_cofactor_charges}")

@pytest.mark.parametrize("benchmark_set, system_name", [
    (benchmark_set, system_name)
    for benchmark_set in list_benchmark_sets()
    for system_name in list_systems(benchmark_set)
])
def test_benchmark_system_components(benchmark_set, system_name):
    """Test loading and validation of BenchmarkSystem components."""
    system = get_benchmark_system(benchmark_set, system_name)

    # Validate protein
    print(system.protein)
    protein = ProteinComponent.from_pdb_file(str(system.protein))
    assert protein.to_rdkit().GetNumAtoms() > 0

    # Validate ligands
    for charge_type, ligand_path in system.ligands.items():
        ligand_supplier = Chem.SDMolSupplier(str(ligand_path), removeHs=False)
        ligands = [SmallMoleculeComponent(mol) for mol in ligand_supplier if mol is not None]
        assert len(ligands) > 0

    # Validate cofactors
    for _, cofactor_path in system.cofactors.items():
        cofactor_supplier = Chem.SDMolSupplier(str(cofactor_path), removeHs=False)
        cofactors = [SmallMoleculeComponent(mol) for mol in cofactor_supplier if mol is not None]
        assert len(cofactors) > 0

    # Validate network
    if system.networks:
        network = LigandNetwork.from_json(file=system.networks[0])
        assert hasattr(network, 'edges')
        assert len(network.edges) > 0
