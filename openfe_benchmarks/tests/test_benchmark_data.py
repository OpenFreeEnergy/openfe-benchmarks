"""
Tests for the benchmark systems module.
"""

import pytest
import json

from gufe.tokenization import JSON_HANDLER
from rdkit import Chem
from openfe import ProteinComponent, SmallMoleculeComponent, LigandNetwork

from openfe_benchmarks.data import (
    BenchmarkData,
    BenchmarkIndex,
    get_benchmark_data_system,
    get_benchmark_set_data_systems,
    PARTIAL_CHARGE_TYPES,
)

from openff.units import unit


def test_list_benchmark_sets_output():
    """Test that BenchmarkIndex lists benchmark sets."""
    index = BenchmarkIndex()
    benchmark_sets = index.list_benchmark_sets()
    assert isinstance(benchmark_sets, list)
    assert all(isinstance(item, str) for item in benchmark_sets)


@pytest.fixture(scope="module")
def benchmark_sets():
    """Fixture to provide all available benchmark sets."""
    index = BenchmarkIndex()
    return index.list_benchmark_sets()


@pytest.mark.parametrize(
    "benchmark_set",
    [benchmark_set for benchmark_set in BenchmarkIndex().list_benchmark_sets()],
)
def test_benchmark_set_systems(benchmark_set):
    """Test that systems in a benchmark set can be retrieved."""
    systems = get_benchmark_set_data_systems(benchmark_set)
    assert isinstance(systems, dict)
    assert all(isinstance(system, str) for system, _ in systems.items())
    assert all(isinstance(system, BenchmarkData) for _, system in systems.items())


@pytest.mark.parametrize(
    "benchmark_set, system_name",
    [
        (benchmark_set, system_name)
        for benchmark_set in BenchmarkIndex().list_benchmark_sets()
        for system_name in get_benchmark_set_data_systems(benchmark_set)
    ],
)
def test_benchmark_system_initialization(benchmark_set, system_name):
    """Test initialization of all BenchmarkData objects."""
    system = get_benchmark_data_system(benchmark_set, system_name)
    assert system is not None
    assert system.name == system_name
    assert system.benchmark_set == benchmark_set

    missing_ligand_charges = [
        charge for charge in PARTIAL_CHARGE_TYPES if charge not in system.ligands
    ]
    if missing_ligand_charges:
        print(
            f"Warning: Missing ligand charge types for {system_name}: {missing_ligand_charges}"
        )

    if system.cofactors:
        missing_cofactor_charges = [
            charge for charge in PARTIAL_CHARGE_TYPES if charge not in system.cofactors
        ]
        if missing_cofactor_charges:
            print(
                f"Warning: Missing cofactor charge types for {system_name}: {missing_cofactor_charges}"
            )


@pytest.mark.parametrize(
    "benchmark_set, system_name",
    [
        (benchmark_set, system_name)
        for benchmark_set in BenchmarkIndex().list_benchmark_sets()
        for system_name in get_benchmark_set_data_systems(benchmark_set)
    ],
)
def test_benchmark_system_components_with_openfe(benchmark_set, system_name):
    """Test loading and validation of BenchmarkData components through OpenFE."""
    system = get_benchmark_data_system(benchmark_set, system_name)

    # Validate protein can be loaded with OpenFE
    if system.protein:
        protein = ProteinComponent.from_pdb_file(str(system.protein))
        assert protein.to_rdkit().GetNumAtoms() > 0, (
            f"Protein for {benchmark_set}/{system_name} has no atoms"
        )

    # Validate ligands can be loaded with RDKit and converted to OpenFE components
    for charge_type, ligand_path in system.ligands.items():
        ligand_supplier = Chem.SDMolSupplier(str(ligand_path), removeHs=False)
        ligands = [
            SmallMoleculeComponent(mol) for mol in ligand_supplier if mol is not None
        ]
        assert len(ligands) > 0, (
            f"No valid ligands loaded from {charge_type} for {benchmark_set}/{system_name}"
        )

        # Verify each ligand has atoms
        for i, ligand in enumerate(ligands):
            rdkit_mol = ligand.to_rdkit()
            assert rdkit_mol.GetNumAtoms() > 0, (
                f"Ligand {i + 1} ({charge_type}) for {benchmark_set}/{system_name} has no atoms"
            )

    # Validate cofactors can be loaded with RDKit and converted to OpenFE components
    if system.cofactors:
        for charge_type, cofactor_path in system.cofactors.items():
            cofactor_supplier = Chem.SDMolSupplier(str(cofactor_path), removeHs=False)
            cofactors = [
                SmallMoleculeComponent(mol)
                for mol in cofactor_supplier
                if mol is not None
            ]
            assert len(cofactors) > 0, (
                f"No valid cofactors loaded from {charge_type} for {benchmark_set}/{system_name}"
            )

            # Verify each cofactor has atoms
            for i, cofactor in enumerate(cofactors):
                rdkit_mol = cofactor.to_rdkit()
                assert rdkit_mol.GetNumAtoms() > 0, (
                    f"Cofactor {i + 1} ({charge_type}) for {benchmark_set}/{system_name} has no atoms"
                )

    # Validate network can be loaded with OpenFE
    if system.ligand_networks:
        for network_name, network_path in system.ligand_networks.items():
            network = LigandNetwork.from_json(file=str(network_path))
            assert hasattr(network, "edges"), (
                f"Network {network_name} for {benchmark_set}/{system_name} has no edges attribute"
            )
            assert len(network.edges) > 0, (
                f"Network {network_name} for {benchmark_set}/{system_name} has no edges"
            )
            assert hasattr(network, "nodes"), (
                f"Network {network_name} for {benchmark_set}/{system_name} has no nodes attribute"
            )
            assert len(network.nodes) > 0, (
                f"Network {network_name} for {benchmark_set}/{system_name} has no nodes"
            )

    # make sure the reference data is present and can be loaded and has dg
    if system.reference_data:
        for ref_path in system.reference_data.values():
            with open(ref_path, "r") as f:
                ref_data = json.load(f, cls=JSON_HANDLER.decoder)
            assert isinstance(ref_data, dict), (
                f"Reference data for {benchmark_set}/{system_name} is not a dict"
            )
            for ref_entry in ref_data.values():
                assert "dg" in ref_entry, (
                    f"Reference data for {benchmark_set}/{system_name} missing 'dg' key"
                )
                # make sure the units are compatible with kilocalories per mole
                dg = ref_entry["dg"]
                assert isinstance(dg, unit.Quantity), (
                    f"'dg' value for {benchmark_set}/{system_name} is not a Quantity"
                )
                # convert to kcal/mol to check compatibility
                dg.to("kcal/mol")


class TestErrorHandling:
    """Tests for error handling."""

    def test_nonexistent_benchmark_set(self):
        """Test accessing nonexistent benchmark set."""
        with pytest.raises(ValueError, match="not found"):
            get_benchmark_data_system("nonexistent_set", "system")

    def test_nonexistent_system(self):
        """Test accessing nonexistent system in valid set."""
        index = BenchmarkIndex()
        sets = index.list_benchmark_sets()
        if sets:
            with pytest.raises(ValueError, match="not found"):
                get_benchmark_data_system(sets[0], "nonexistent_system")

    def test_invalid_set_name_in_get_benchmark_set_data_systems(self):
        """Test get_benchmark_set_data_systems with invalid set."""
        with pytest.raises(ValueError, match="not found"):
            get_benchmark_set_data_systems("nonexistent_set")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
