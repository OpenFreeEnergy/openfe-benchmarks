import pytest

from openfe_benchmarks import tyk2


class TestTyk2System:
    """
    Tests to check that the tyk2 system is properly created
    """
    @pytest.fixture()
    def system(self):
        return tyk2.get_system()

    def test_ligand_components(self, system):
        assert len(system.ligand_components) == 16

    def test_edges(self, system):
        assert len(system.ligand_network.edges) == 24

    def test_protein_component(self, system):
        assert system.protein_component.name == "tyk2"

    def test_solvent_component(self, system):
        # TODO use name once we bump up to the next gufe release
        assert system.solvent_component.smiles == "O"
