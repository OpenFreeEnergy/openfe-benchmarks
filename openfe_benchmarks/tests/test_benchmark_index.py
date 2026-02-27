"""
Tests for BenchmarkIndex class and benchmark system validation.
"""

import pytest
from pathlib import Path

from openfe_benchmarks.data import (
    BenchmarkIndex,
)

_BASE_DIR = Path(__file__).resolve().parent.parent / "data" / "benchmark_systems"

# Each entry is (tag, [candidate_files]).
# For tags with files: systems WITH the tag must contain at least one of the listed files;
# systems WITHOUT the tag must contain none of them.
# Use an empty list if no files are associated with the tag.
TAG_CHECKS = [
    ("protein", ["protein.pdb"]),
    ("cofactor", ["cofactors.sdf"]),
    ("bfe", ["experimental_binding_data.json"]),
    ("sfe", ["experimental_solvation_free_energy_data.json", "systems_data.json"]),
]


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

        assert initial_systems == reloaded_systems, (
            "Reloaded systems should match initial systems"
        )

    def test_list_systems_by_tag(self):
        """Test that systems can be filtered by tags."""
        index = BenchmarkIndex()
        tags = index.list_available_tags()

        for tag in tags:
            systems = index.list_systems_by_tag([tag])
            assert isinstance(systems, list), (
                f"Systems for tag '{tag}' should be a list"
            )
            for benchmark_set, system_name in systems:
                assert isinstance(benchmark_set, str)
                assert isinstance(system_name, str)

    def test_pytest_tag_defined_criteria(self):
        """Test that systems can be filtered by tags."""
        index = BenchmarkIndex()
        tags = index.list_available_tags()
        checked_tags = [x[0] for x in TAG_CHECKS]
        for tag in tags:
            assert tag in checked_tags

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
            assert systems_both_tags == systems_tag1 & systems_tag2, (
                f"Systems for tags '{tag1}' and '{tag2}' should be the intersection of their individual systems"
            )

    def test_list_systems_by_nonexistent_tag(self):
        """Test that listing systems by a nonexistent tag returns an empty list."""
        index = BenchmarkIndex()
        nonexistent_tag = "nonexistent_tag"
        systems = index.list_systems_by_tag([nonexistent_tag])

        assert systems == [], (
            f"Systems for nonexistent tag '{nonexistent_tag}' should be an empty list"
        )

    def test_list_systems_by_empty_tag_list(self):
        """Test that listing systems by an empty tag list returns all systems."""
        index = BenchmarkIndex()
        all_systems = [
            (benchmark_set, system_name)
            for benchmark_set, systems in index._data["systems"].items()
            for system_name in systems
        ]
        systems = index.list_systems_by_tag([])

        assert set(systems) == set(all_systems), (
            "Systems for an empty tag list should include all systems in the index"
        )

    @pytest.mark.parametrize("tag, candidate_files", TAG_CHECKS)
    def test_files_for_tags(self, tag, candidate_files):
        """Check that systems with a tag have at least one of its candidate files,
        and systems without the tag have none of them."""
        index = BenchmarkIndex()

        if not candidate_files:
            return

        all_systems = [
            (benchmark_set, system_name)
            for benchmark_set, systems in index._data["systems"].items()
            for system_name in systems
        ]

        for benchmark_set, system_name in all_systems:
            system_path = _BASE_DIR / benchmark_set / system_name
            has_file = any(
                any(system_path.glob(pattern)) for pattern in candidate_files
            )
            has_tag = tag in index._data["systems"][benchmark_set][system_name]

            if has_tag:
                assert has_file, (
                    f"System '{system_name}' in '{benchmark_set}' is tagged '{tag}' but "
                    f"has none of {candidate_files}"
                )
            else:
                assert not has_file, (
                    f"System '{system_name}' in '{benchmark_set}' has a '{tag}' candidate "
                    f"file but is not tagged '{tag}'"
                )

    def test_all_systems_in_index(self):
        """Test that all ligands.sdf and PREPARATION_DETAILS.md files in the repo are associated with systems in the index."""
        index = BenchmarkIndex()

        # Find all relevant files in the repository
        all_relevant_files = list(_BASE_DIR.rglob("ligands.sdf")) + list(
            _BASE_DIR.rglob("preparation_details.md")
        )

        # Build a set of expected paths from the index
        indexed_systems = set()
        for benchmark_set, systems in index._data["systems"].items():
            for system_name in systems:
                system_path = _BASE_DIR / benchmark_set / system_name
                indexed_systems.add(system_path)

        # Check that all found relevant files are in the index
        nonindexed_files = []
        for relevant_file in all_relevant_files:
            parent_dir = relevant_file.parent
            if parent_dir not in indexed_systems:
                # Get relative path for better error message
                rel_path = relevant_file.relative_to(_BASE_DIR)
                nonindexed_files.append(str(rel_path))

        assert len(nonindexed_files) == 0, pytest.fail(
            f"Found {len(nonindexed_files)} relevant file(s) not indexed in benchmark_system_indexing.yml: "
            + ", ".join(nonindexed_files)
        )

        # Also verify the reverse: that indexed systems actually have both of the relevant files
        missing_files = []
        for system_path in indexed_systems:
            if not system_path.exists():
                has_relevant_file = False
            else:
                present = {p.name.lower() for p in system_path.iterdir() if p.is_file()}
                has_relevant_file = all(
                    fname.lower() in present
                    for fname in ["ligands.sdf", "preparation_details.md"]
                )
            if not has_relevant_file:
                rel_path = system_path.relative_to(_BASE_DIR)
                missing_files.append(str(rel_path))

        assert len(missing_files) == 0, pytest.fail(
            f"Found {len(missing_files)} indexed system(s) missing relevant files: "
            + ", ".join(missing_files)
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
