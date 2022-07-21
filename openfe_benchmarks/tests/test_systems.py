import pytest

from openfe_benchmarks import tyk2, ptp1b


class TestTyk2System:
    """
    Tests to check that the tyk2 system is properly created
    """
    @pytest.fixture()
    def system(self):
        return tyk2.get_system()

    def test_ligand_components(self, system):
        assert len(system.ligand_components) == 13

    def test_edges(self, system):
        assert len(system.ligand_network.edges) == 12

    def test_protein_component(self, system):
        assert system.protein_component.name == "tyk2"

    def test_solvent_component(self, system):
        # TODO use name once we bump up to the next gufe release
        assert system.solvent_component.smiles == "O"


class TestPTP1BSystem:
    """
    Tests to check that the ptp1b system is properly created
    """
    @pytest.fixture()
    def system(self):
        return ptp1b.get_system()

    def test_ligand_components(self, system):
        assert len(system.ligand_components) == 23

    def test_edges(self, system):
        assert len(system.ligand_network.edges) == 49

    def test_protein_component(self, system):
        assert system.protein_component.name == "ptp1b"

    def test_solvent_component(self, system):
        # TODO use name once we bump up to the next gufe release
        assert system.solvent_component.smiles == "O"
