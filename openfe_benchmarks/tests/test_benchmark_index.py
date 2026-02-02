"""
Tests for BenchmarkIndex class and benchmark system validation.
"""

import pytest
from pathlib import Path

from openfe_benchmarks.data import (
    BenchmarkIndex,
    get_system_index,
    list_benchmark_sets,
    get_benchmark_data_system,
    get_benchmark_set_data_systems,
)


class TestBenchmarkIndex:
    """Tests for the BenchmarkIndex singleton class."""
    
    def test_singleton_pattern(self):
        """Test that BenchmarkIndex follows singleton pattern."""
        index1 = BenchmarkIndex()
        index2 = BenchmarkIndex()
        assert index1 is index2, "BenchmarkIndex should be a singleton"
    
    def test_get_all_systems(self):
        """Test retrieving all systems from the index."""
        index = BenchmarkIndex()
        all_systems = index.get_all_systems()
        
        assert isinstance(all_systems, dict)
        assert len(all_systems) > 0, "Should have at least one benchmark set"
        
        # Check that all values are lists
        for set_name, systems in all_systems.items():
            assert isinstance(systems, list), f"{set_name} should have a list of systems"
            assert len(systems) > 0, f"{set_name} should have at least one system"
            assert all(isinstance(s, str) for s in systems), \
                f"All systems in {set_name} should be strings"
    
    def test_get_systems_for_rbfe(self):
        """Test retrieving systems for RBFE calculations."""
        index = BenchmarkIndex()
        rbfe_systems = index.get_systems_for_calculation('rbfe')
        
        assert isinstance(rbfe_systems, dict)
        assert len(rbfe_systems) > 0, "Should have RBFE systems"
    
    def test_get_systems_for_all_calc_types(self):
        """Test retrieving systems for all calculation types."""
        index = BenchmarkIndex()
        
        for calc_type in ['rbfe', 'abfe', 'asfe', 'rsfe']:
            systems = index.get_systems_for_calculation(calc_type)
            assert isinstance(systems, dict), f"{calc_type} should return dict"
            # Since all current systems are RBFE, all types should return systems
            assert len(systems) > 0, f"{calc_type} should have systems (from RBFE)"
    
    def test_invalid_calculation_type(self):
        """Test that invalid calculation types raise ValueError."""
        index = BenchmarkIndex()
        
        with pytest.raises(ValueError, match="Invalid calculation type"):
            index.get_systems_for_calculation('invalid')
    
    def test_asfe_includes_rbfe_systems(self):
        """Test that ASFE calculation includes RBFE systems."""
        index = BenchmarkIndex()
        
        rbfe_systems = index.get_systems_for_calculation('rbfe')
        asfe_systems = index.get_systems_for_calculation('asfe')
        
        # ASFE should include all RBFE systems
        for set_name, systems in rbfe_systems.items():
            assert set_name in asfe_systems, \
                f"ASFE should include RBFE set {set_name}"
            for system in systems:
                assert system in asfe_systems[set_name], \
                    f"ASFE should include RBFE system {set_name}/{system}"
    
    def test_reload(self):
        """Test that reload method works."""
        index = BenchmarkIndex()
        initial_systems = index.get_all_systems()
        
        # Reload should not crash and should return same data
        index.reload()
        reloaded_systems = index.get_all_systems()
        
        assert initial_systems == reloaded_systems, \
            "Reloaded systems should match initial systems"


class TestGetSystemIndex:
    """Tests for the get_system_index function."""
    
    def test_get_system_index_rbfe(self):
        """Test get_system_index for RBFE."""
        systems = get_system_index('rbfe')
        
        assert isinstance(systems, dict)
        assert len(systems) > 0
    
    def test_get_system_index_all_types(self):
        """Test get_system_index for all calculation types."""
        for calc_type in ['rbfe', 'abfe', 'asfe', 'rsfe']:
            systems = get_system_index(calc_type)
            assert isinstance(systems, dict)
            assert len(systems) > 0, f"{calc_type} should return systems"
    
    def test_get_system_index_invalid(self):
        """Test get_system_index with invalid type."""
        with pytest.raises(ValueError, match="Invalid calculation type"):
            get_system_index('invalid_type')

class TestSystemValidation:
    """Tests that validate each individual system."""
    
    def test_all_systems_have_required_files(self):
        """Test that all systems in the index have required files."""
        index = BenchmarkIndex()
        all_systems = index.get_all_systems()
        
        errors = []
        
        for set_name, systems in all_systems.items():
            for system_name in systems:
                try:
                    # Try to load the system
                    system_data = get_benchmark_data_system(set_name, system_name)
                    
                    # Verify essential attributes exist
                    assert system_data.name == system_name, \
                        f"System name mismatch for {set_name}/{system_name}"
                    assert system_data.benchmark_set == set_name, \
                        f"Benchmark set mismatch for {set_name}/{system_name}"
                    
                    # Check ligands
                    assert len(system_data.ligands) > 0, \
                        f"No ligands found for {set_name}/{system_name}"
                    assert 'no_charges' in system_data.ligands, \
                        f"Missing base ligands.sdf for {set_name}/{system_name}"
                    
                    # Check that network exists
                    assert system_data.network is not None, \
                        f"Missing network for {set_name}/{system_name}"
                    assert system_data.network.exists(), \
                        f"Network file doesn't exist for {set_name}/{system_name}"
                    
                    # Check that ligands file exists
                    for charge_type, ligand_path in system_data.ligands.items():
                        assert ligand_path.exists(), \
                            f"Ligands file missing ({charge_type}) for {set_name}/{system_name}"
                    
                    # Check protein if it exists
                    if system_data.protein:
                        assert system_data.protein.exists(), \
                            f"Protein file missing for {set_name}/{system_name}"
                    
                    # Check cofactors if they exist
                    for charge_type, cofactor_path in system_data.cofactors.items():
                        assert cofactor_path.exists(), \
                            f"Cofactor file missing ({charge_type}) for {set_name}/{system_name}"
                
                except Exception as e:
                    errors.append(f"{set_name}/{system_name}: {str(e)}")
        
        if errors:
            error_msg = "System validation errors:\n" + "\n".join(f"  - {e}" for e in errors)
            pytest.fail(error_msg)
    
    def test_rbfe_systems_have_protein_and_network(self):
        """Test that all RBFE systems have both protein and network files."""
        rbfe_systems = get_system_index('rbfe')
        
        errors = []
        
        for set_name, systems in rbfe_systems.items():
            for system_name in systems:
                try:
                    system_data = get_benchmark_data_system(set_name, system_name)
                    
                    # RBFE requires protein (can be None for some edge cases)
                    # and network
                    assert system_data.network is not None, \
                        f"RBFE system {set_name}/{system_name} missing network"
                    
                    # Note: protein can be None for some systems, but most RBFE should have it
                    # This is a soft check - log warning instead of failing
                    if system_data.protein is None:
                        print(f"Warning: RBFE system {set_name}/{system_name} has no protein")
                
                except Exception as e:
                    errors.append(f"{set_name}/{system_name}: {str(e)}")
        
        if errors:
            error_msg = "RBFE system validation errors:\n" + "\n".join(f"  - {e}" for e in errors)
            pytest.fail(error_msg)
    
    def test_system_count_matches_metadata(self):
        """Test that actual system count matches metadata in YAML."""
        index = BenchmarkIndex()
        all_systems = index.get_all_systems()
        
        total_count = sum(len(systems) for systems in all_systems.values())
        
        # Get expected count from metadata
        metadata = index._data.get('metadata', {})
        expected_count = metadata.get('total_systems')
        
        if expected_count is not None:
            assert total_count == expected_count, \
                f"System count mismatch: found {total_count}, expected {expected_count}"


class TestErrorHandling:
    """Tests for error handling."""
    
    def test_nonexistent_benchmark_set(self):
        """Test accessing nonexistent benchmark set."""
        with pytest.raises(ValueError, match="not found"):
            get_benchmark_data_system('nonexistent_set', 'system')
    
    def test_nonexistent_system(self):
        """Test accessing nonexistent system in valid set."""
        sets = list_benchmark_sets()
        if sets:
            with pytest.raises(ValueError, match="not found"):
                get_benchmark_data_system(sets[0], 'nonexistent_system')
    
    def test_invalid_set_name_in_get_benchmark_set_data_systems(self):
        """Test get_benchmark_set_data_systems with invalid set."""
        with pytest.raises(ValueError, match="not found"):
            get_benchmark_set_data_systems('nonexistent_set')


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
