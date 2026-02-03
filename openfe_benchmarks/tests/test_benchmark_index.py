"""
Tests for BenchmarkIndex class and benchmark system validation.
"""

import pytest
from pathlib import Path

from openfe_benchmarks.data import (
    BenchmarkIndex,
)
_BASE_DIR = Path(__file__).resolve().parent.parent / "data" / "benchmark_systems"

class TestBenchmarkIndex:
    """Tests for the BenchmarkIndex singleton class."""
    
    def test_singleton_pattern(self):
        """Test that BenchmarkIndex follows singleton pattern."""
        index1 = BenchmarkIndex()
        index2 = BenchmarkIndex()
        assert index1 is index2, "BenchmarkIndex should be a singleton"

    def test_reload(self):
        """Test that reload method works."""
        index = BenchmarkIndex()
        initial_systems = index.list_benchmark_sets()

        # Reload should not crash and should return same data
        index.reload()
        reloaded_systems = index.list_benchmark_sets()

        assert initial_systems == reloaded_systems, \
            "Reloaded systems should match initial systems"

    def test_list_systems_by_tag(self):
        """Test that systems can be filtered by tags."""
        index = BenchmarkIndex()
        tags = index.list_available_tags()

        for tag in tags:
            systems = index.list_systems_by_tag([tag])
            assert isinstance(systems, list), f"Systems for tag '{tag}' should be a list"
            for benchmark_set, system_name in systems:
                assert isinstance(benchmark_set, str)
                assert isinstance(system_name, str)

    def test_list_systems_by_multiple_tags(self):
        """Test that systems can be filtered by multiple tags."""
        index = BenchmarkIndex()
        tags = list(index.list_available_tags())

        if len(tags) > 1:
            tag1, tag2 = tags[0], tags[1]
            systems_tag1 = set(index.list_systems_by_tag([tag1]))
            systems_tag2 = set(index.list_systems_by_tag([tag2]))
            systems_both_tags = set(index.list_systems_by_tag([tag1, tag2]))

            # The systems for both tags should be the intersection of the two sets
            assert systems_both_tags == systems_tag1 & systems_tag2, \
                f"Systems for tags '{tag1}' and '{tag2}' should be the intersection of their individual systems"

    def test_list_systems_by_nonexistent_tag(self):
        """Test that listing systems by a nonexistent tag returns an empty list."""
        index = BenchmarkIndex()
        nonexistent_tag = "nonexistent_tag"
        systems = index.list_systems_by_tag([nonexistent_tag])

        assert systems == [], f"Systems for nonexistent tag '{nonexistent_tag}' should be an empty list"

    def test_list_systems_by_empty_tag_list(self):
        """Test that listing systems by an empty tag list returns all systems."""
        index = BenchmarkIndex()
        all_systems = [
            (benchmark_set, system_name)
            for benchmark_set, systems in index._data['systems'].items()
            for system_name in systems
        ]
        systems = index.list_systems_by_tag([])

        assert set(systems) == set(all_systems), \
            "Systems for an empty tag list should include all systems in the index"
        

    def test_protein_files_for_bfe(self):
        """Check that systems with 'bfe' tag have protein.pdb files."""
        index = BenchmarkIndex()
        bfe_systems = index.list_systems_by_tag(['bfe'])

        for benchmark_set, system_name in bfe_systems:
            system_path = _BASE_DIR / benchmark_set / system_name
            assert (system_path / 'protein.pdb').exists(), \
                f"System '{system_name}' in '{benchmark_set}' should have protein.pdb"

        non_bfe_systems = [
            (benchmark_set, system_name)
            for benchmark_set, systems in index._data['systems'].items()
            for system_name in systems
            if 'bfe' not in systems[system_name]
        ]

        for benchmark_set, system_name in non_bfe_systems:
            system_path = _BASE_DIR / benchmark_set / system_name
            assert not (system_path / 'protein.pdb').exists(), \
                f"System '{system_name}' in '{benchmark_set}' should have 'bfe' tag"

    def test_cofactor_files_for_cofactor_tag(self):
        """Check that systems with 'cofactor' tag have cofactors.sdf files."""
        index = BenchmarkIndex()
        cofactor_systems = index.list_systems_by_tag(['cofactor'])

        for benchmark_set, system_name in cofactor_systems:
            system_path = _BASE_DIR / benchmark_set / system_name
            assert (system_path / 'cofactors.sdf').exists(), \
                f"System '{system_name}' in '{benchmark_set}' should have cofactors.sdf"

        non_cofactor_systems = [
            (benchmark_set, system_name)
            for benchmark_set, systems in index._data['systems'].items()
            for system_name in systems
            if 'cofactor' not in systems[system_name]
        ]

        for benchmark_set, system_name in non_cofactor_systems:
            system_path = _BASE_DIR / benchmark_set / system_name
            assert not (system_path / 'cofactors.sdf').exists(), \
                f"System '{system_name}' in '{benchmark_set}' should have 'cofactor' tag"

    def test_ligands_files_for_sfe_tag(self):
        """Check that systems with 'ligands.sdf' have 'sfe' tag and vice versa."""
        index = BenchmarkIndex()
        sfe_systems = index.list_systems_by_tag(['sfe'])

        for benchmark_set, system_name in sfe_systems:
            system_path = _BASE_DIR / benchmark_set / system_name
            assert (system_path / 'ligands.sdf').exists(), \
                f"System '{system_name}' in '{benchmark_set}' should have ligands.sdf"

        all_systems = [
            (benchmark_set, system_name)
            for benchmark_set, systems in index._data['systems'].items()
            for system_name in systems
        ]

        for benchmark_set, system_name in all_systems:
            system_path = _BASE_DIR / benchmark_set / system_name
            has_ligands = (system_path / 'ligands.sdf').exists()
            has_sfe_tag = 'sfe' in index._data['systems'][benchmark_set][system_name]

            assert has_ligands == has_sfe_tag, \
                f"Mismatch for system '{system_name}' in '{benchmark_set}': " \
                f"has_ligands={has_ligands}, has_sfe_tag={has_sfe_tag}"

    def test_all_ligands_sdf_files_in_index(self):
        """Test that all ligands.sdf files in the repo are associated with systems in the index."""
        index = BenchmarkIndex()
        
        # Find all ligands.sdf files in the repository
        all_ligands_files = list(_BASE_DIR.rglob('ligands.sdf'))
        
        # Build a set of expected paths from the index
        indexed_systems = set()
        for benchmark_set, systems in index._data['systems'].items():
            for system_name in systems:
                system_path = _BASE_DIR / benchmark_set / system_name / 'ligands.sdf'
                indexed_systems.add(system_path)
        
        # Check that all found ligands.sdf files are in the index
        nonindexed_files = []
        for ligands_file in all_ligands_files:
            if ligands_file not in indexed_systems:
                # Get relative path for better error message
                rel_path = ligands_file.relative_to(_BASE_DIR)
                nonindexed_files.append(str(rel_path))
        
        assert len(nonindexed_files) == 0, \
            pytest.fail(
                f"Found {len(nonindexed_files)} ligands.sdf file(s) not indexed in benchmark_system_indexing.yml:"
                ", ".join(f"{f}" for f in nonindexed_files)
            )
        
        # Also verify the reverse: that indexed systems with ligands.sdf actually have the file
        missing_files = []
        for system_path in indexed_systems:
            if not system_path.exists():
                rel_path = system_path.relative_to(_BASE_DIR)
                missing_files.append(str(rel_path))
        
        assert len(missing_files) == 0, \
            pytest.fail(
                f"Found {len(missing_files)} ligands.sdf file(s) indexed but not present on disk:"
                ", ".join(f"{f}" for f in nonindexed_files)
            )



if __name__ == '__main__':
    pytest.main([__file__, '-v'])
