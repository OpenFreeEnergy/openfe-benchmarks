"""
Tests for the industry benchmark systems module.
"""
import pytest
from pathlib import Path
from openfe_benchmarks.data import (
    get_benchmark_system,
    list_benchmark_sets,
    list_systems,
    BenchmarkSystem,
    PARTIAL_CHARGE_TYPES,
)


class TestListBenchmarkSets:
    """Tests for listing available benchmark sets."""
    
    def test_list_benchmark_sets_returns_list(self):
        """Test that list_benchmark_sets returns a list."""
        result = list_benchmark_sets()
        assert isinstance(result, list)
    
    def test_list_benchmark_sets_not_empty(self):
        """Test that there is at least one benchmark set."""
        result = list_benchmark_sets()
        assert len(result) > 0
    
    def test_list_benchmark_sets_contains_known_sets(self):
        """Test that known benchmark sets are present."""
        result = list_benchmark_sets()
        # These should exist based on the directory structure
        expected_sets = ['fragments', 'jacs_set', 'janssen_bace', 'mcs_docking_set']
        for expected in expected_sets:
            assert expected in result
    
    def test_list_benchmark_sets_sorted(self):
        """Test that benchmark sets are returned in sorted order."""
        result = list_benchmark_sets()
        assert result == sorted(result)


class TestListSystems:
    """Tests for listing systems within a benchmark set."""
    
    def test_list_systems_jacs_set(self):
        """Test listing systems in jacs_set."""
        result = list_systems('jacs_set')
        assert isinstance(result, list)
        assert len(result) > 0
        # Known systems in jacs_set
        assert 'p38' in result
        assert 'tyk2' in result
    
    def test_list_systems_fragments(self):
        """Test listing systems in fragments."""
        result = list_systems('fragments')
        assert isinstance(result, list)
        assert len(result) > 0
    
    def test_list_systems_sorted(self):
        """Test that systems are returned in sorted order."""
        result = list_systems('jacs_set')
        assert result == sorted(result)
    
    def test_list_systems_invalid_benchmark_set(self):
        """Test that listing systems for invalid benchmark set raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            list_systems('nonexistent_set')
        
        assert "not found" in str(exc_info.value)
        assert "Available benchmark sets:" in str(exc_info.value)


class TestGetBenchmarkSystem:
    """Tests for the get_benchmark_system factory method."""
    
    def test_get_benchmark_system_returns_benchmark_system(self):
        """Test that get_benchmark_system returns a BenchmarkSystem object."""
        system = get_benchmark_system('jacs_set', 'p38')
        assert isinstance(system, BenchmarkSystem)
    
    def test_get_benchmark_system_attributes(self):
        """Test that BenchmarkSystem has all required attributes."""
        system = get_benchmark_system('jacs_set', 'p38')
        assert hasattr(system, 'name')
        assert hasattr(system, 'benchmark_set')
        assert hasattr(system, 'protein')
        assert hasattr(system, 'ligands')
        assert hasattr(system, 'cofactors')
        assert hasattr(system, 'mappings')
    
    def test_get_benchmark_system_name_and_set(self):
        """Test that system name and benchmark set are correctly set."""
        system = get_benchmark_system('jacs_set', 'p38')
        assert system.name == 'p38'
        assert system.benchmark_set == 'jacs_set'
    
    def test_get_benchmark_system_protein_path(self):
        """Test that protein path is valid and points to protein.pdb."""
        system = get_benchmark_system('jacs_set', 'p38')
        assert isinstance(system.protein, Path)
        assert system.protein.exists()
        assert system.protein.name == 'protein.pdb'
        assert system.protein.suffix == '.pdb'
    
    def test_get_benchmark_system_ligands_dict(self):
        """Test that ligands is a dict with charge types as keys."""
        system = get_benchmark_system('jacs_set', 'p38')
        assert isinstance(system.ligands, dict)
        assert len(system.ligands) > 0
        
        # Check that all keys are valid charge types
        for charge_type in system.ligands.keys():
            assert charge_type in PARTIAL_CHARGE_TYPES
        
        # Check that all values are valid paths
        for ligand_path in system.ligands.values():
            assert isinstance(ligand_path, Path)
            assert ligand_path.exists()
            assert ligand_path.suffix == '.sdf'
    
    def test_get_benchmark_system_ligands_naming(self):
        """Test that ligand files follow the correct naming convention."""
        system = get_benchmark_system('jacs_set', 'p38')
        
        for charge_type, ligand_path in system.ligands.items():
            expected_name = f'ligands_{charge_type}.sdf'
            assert ligand_path.name == expected_name
    
    def test_get_benchmark_system_cofactors_dict(self):
        """Test that cofactors is a dict (may be empty)."""
        system = get_benchmark_system('jacs_set', 'p38')
        assert isinstance(system.cofactors, dict)
        
        # If cofactors exist, validate them
        for charge_type, cofactor_path in system.cofactors.items():
            assert charge_type in PARTIAL_CHARGE_TYPES
            assert isinstance(cofactor_path, Path)
            assert cofactor_path.exists()
            assert cofactor_path.suffix == '.sdf'
            assert cofactor_path.name.startswith('cofactors_')
    
    def test_get_benchmark_system_mappings_list(self):
        """Test that mappings is a list of paths."""
        system = get_benchmark_system('jacs_set', 'p38')
        assert isinstance(system.mappings, list)
        
        # If mappings exist, validate them
        for mapping_path in system.mappings:
            assert isinstance(mapping_path, Path)
            assert mapping_path.exists()
    
    def test_get_benchmark_system_invalid_set(self):
        """Test that invalid benchmark set raises ValueError with available sets."""
        with pytest.raises(ValueError) as exc_info:
            get_benchmark_system('invalid_set', 'p38')
        
        error_msg = str(exc_info.value)
        assert "Benchmark set 'invalid_set' not found" in error_msg
        assert "Available benchmark sets:" in error_msg
    
    def test_get_benchmark_system_invalid_system(self):
        """Test that invalid system name raises ValueError with available systems."""
        with pytest.raises(ValueError) as exc_info:
            get_benchmark_system('jacs_set', 'nonexistent_system')
        
        error_msg = str(exc_info.value)
        assert "System 'nonexistent_system' not found in benchmark set 'jacs_set'" in error_msg
        assert "Available systems in 'jacs_set':" in error_msg
    
    def test_get_benchmark_system_multiple_charge_types(self):
        """Test that systems can have multiple charge type files."""
        system = get_benchmark_system('jacs_set', 'p38')
        
        # Should have at least antechamber_am1bcc
        assert 'antechamber_am1bcc' in system.ligands
        
        # If multiple charge types exist, they should point to different files
        if len(system.ligands) > 1:
            assert 'openeye_elf10' in system.ligands
            assert system.ligands['antechamber_am1bcc'] != system.ligands['openeye_elf10']
    
    def test_get_benchmark_system_repr(self):
        """Test that BenchmarkSystem has a useful string representation."""
        system = get_benchmark_system('jacs_set', 'p38')
        repr_str = repr(system)
        
        assert 'BenchmarkSystem' in repr_str
        assert 'p38' in repr_str
        assert 'jacs_set' in repr_str


class TestBenchmarkSystemValidation:
    """Tests for validation of benchmark system files and structure."""
    
    def test_all_systems_have_protein(self):
        """Test that all systems in all benchmark sets have a protein.pdb file."""
        for benchmark_set in list_benchmark_sets():
            for system_name in list_systems(benchmark_set):
                system = get_benchmark_system(benchmark_set, system_name)
                assert system.protein is not None
                assert system.protein.exists()
                assert system.protein.name == 'protein.pdb'
    
    def test_all_systems_have_ligands(self):
        """Test that all systems have at least one ligand file."""
        for benchmark_set in list_benchmark_sets():
            for system_name in list_systems(benchmark_set):
                system = get_benchmark_system(benchmark_set, system_name)
                assert len(system.ligands) > 0
    
    def test_ligand_charge_types_are_valid(self):
        """Test that all ligand charge types are in PARTIAL_CHARGE_TYPES."""
        for benchmark_set in list_benchmark_sets():
            for system_name in list_systems(benchmark_set):
                system = get_benchmark_system(benchmark_set, system_name)
                for charge_type in system.ligands.keys():
                    assert charge_type in PARTIAL_CHARGE_TYPES
    
    def test_cofactor_charge_types_are_valid(self):
        """Test that all cofactor charge types are in PARTIAL_CHARGE_TYPES."""
        for benchmark_set in list_benchmark_sets():
            for system_name in list_systems(benchmark_set):
                system = get_benchmark_system(benchmark_set, system_name)
                for charge_type in system.cofactors.keys():
                    assert charge_type in PARTIAL_CHARGE_TYPES


class TestPartialChargeTypes:
    """Tests for PARTIAL_CHARGE_TYPES constant."""
    
    def test_partial_charge_types_is_list(self):
        """Test that PARTIAL_CHARGE_TYPES is a list."""
        assert isinstance(PARTIAL_CHARGE_TYPES, list)
    
    def test_partial_charge_types_not_empty(self):
        """Test that PARTIAL_CHARGE_TYPES is not empty."""
        assert len(PARTIAL_CHARGE_TYPES) > 0
    
    def test_partial_charge_types_contains_expected(self):
        """Test that PARTIAL_CHARGE_TYPES contains expected charge types."""
        assert 'antechamber_am1bcc' in PARTIAL_CHARGE_TYPES
        assert 'openeye_elf10' in PARTIAL_CHARGE_TYPES


class TestSystemWithCofactors:
    """Tests for systems that include cofactors."""
    
    def test_cofactor_system_exists(self):
        """Test that we can detect systems with cofactors if they exist."""
        found_cofactors = False
        
        for benchmark_set in list_benchmark_sets():
            for system_name in list_systems(benchmark_set):
                system = get_benchmark_system(benchmark_set, system_name)
                if system.cofactors:
                    found_cofactors = True
                    # If we find cofactors, verify they're properly structured
                    for charge_type, cofactor_path in system.cofactors.items():
                        assert charge_type in PARTIAL_CHARGE_TYPES
                        assert cofactor_path.exists()
                    break
            if found_cofactors:
                break
        
        # Note: Not all benchmark sets may have cofactors, so this test
        # only validates the structure if cofactors are present
    
    def test_cofactor_naming_convention(self):
        """Test that cofactor files follow the naming convention."""
        for benchmark_set in list_benchmark_sets():
            for system_name in list_systems(benchmark_set):
                system = get_benchmark_system(benchmark_set, system_name)
                
                for charge_type, cofactor_path in system.cofactors.items():
                    expected_name = f'cofactors_{charge_type}.sdf'
                    assert cofactor_path.name == expected_name


class TestKnownSystems:
    """Tests for specific known systems."""
    
    def test_jacs_set_p38_system(self):
        """Test the jacs_set/p38 system specifically."""
        system = get_benchmark_system('jacs_set', 'p38')
        
        assert system.name == 'p38'
        assert system.benchmark_set == 'jacs_set'
        assert system.protein.exists()
        assert len(system.ligands) >= 1
        assert len(system.mappings) >= 1
    
    def test_jacs_set_tyk2_system(self):
        """Test the jacs_set/tyk2 system specifically."""
        system = get_benchmark_system('jacs_set', 'tyk2')
        
        assert system.name == 'tyk2'
        assert system.benchmark_set == 'jacs_set'
        assert system.protein.exists()
        assert len(system.ligands) >= 1
    
    def test_fragments_system(self):
        """Test a system from the fragments benchmark set."""
        # Get first available system from fragments
        fragment_systems = list_systems('fragments')
        assert len(fragment_systems) > 0, "No systems found in fragments benchmark set"
        
        system = get_benchmark_system('fragments', fragment_systems[0])
        
        assert system.name == fragment_systems[0]
        assert system.benchmark_set == 'fragments'
        assert system.protein.exists()
        assert len(system.ligands) >= 1
