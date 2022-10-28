import pytest

from openfe_benchmarks import (
        tyk2, ptp1b, cmet, tnsk2, mcl1, p38, thrombin, hif2a, syk
)


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
        assert len(system.ligand_components) == 22

    def test_edges(self, system):
        assert len(system.ligand_network.edges) == 21

    def test_protein_component(self, system):
        assert system.protein_component.name == "ptp1b"

    def test_solvent_component(self, system):
        # TODO use name once we bump up to the next gufe release
        assert system.solvent_component.smiles == "O"


class TestCMETSystem:
    """
    Tests to check that the cmet system is properly created
    """
    @pytest.fixture()
    def system(self):
        return cmet.get_system()

    def test_ligand_components(self, system):
        assert len(system.ligand_components) == 5

    def test_edges(self, system):
        assert len(system.ligand_network.edges) == 4

    def test_protein_component(self, system):
        assert system.protein_component.name == "cmet"

    def test_solvent_component(self, system):
        # TODO use name once we bump up to the next gufe release
        assert system.solvent_component.smiles == "O"


class TestTNSK2System:
    """
    Tests to check that the tnsk2 system is properly created
    """
    @pytest.fixture()
    def system(self):
        return tnsk2.get_system()

    def test_ligand_components(self, system):
        assert len(system.ligand_components) == 27

    def test_edges(self, system):
        assert len(system.ligand_network.edges) == 26

    def test_protein_component(self, system):
        assert system.protein_component.name == "tnsk2"

    def test_solvent_component(self, system):
        # TODO use name once we bump up to the next gufe release
        assert system.solvent_component.smiles == "O"


class TestMCL1System:
    """
    Tests to check that the mcl1 system is properly created
    """
    @pytest.fixture()
    def system(self):
        return mcl1.get_system()

    def test_ligand_components(self, system):
        assert len(system.ligand_components) == 24

    def test_edges(self, system):
        assert len(system.ligand_network.edges) == 23

    def test_protein_component(self, system):
        assert system.protein_component.name == "mcl1"

    def test_solvent_component(self, system):
        # TODO use name once we bump up to the next gufe release
        assert system.solvent_component.smiles == "O"


class TestP38System:
    """
    Tests to check that the p38 system is properly created
    """
    @pytest.fixture()
    def system(self):
        return p38.get_system()

    def test_ligand_components(self, system):
        assert len(system.ligand_components) == 29

    def test_edges(self, system):
        assert len(system.ligand_network.edges) == 28

    def test_protein_component(self, system):
        assert system.protein_component.name == "p38"

    def test_solvent_component(self, system):
        # TODO use name once we bump up to the next gufe release
        assert system.solvent_component.smiles == "O"


class TestP38OldSystem:
    """
    Tests to check that the p38 system is properly created
    """
    @pytest.fixture()
    def system(self):
        return p38.get_old_system()

    def test_ligand_components(self, system):
        assert len(system.ligand_components) == 34

    def test_edges(self, system):
        assert len(system.ligand_network.edges) == 56

    def test_protein_component(self, system):
        assert system.protein_component.name == "p38"

    def test_solvent_component(self, system):
        # TODO use name once we bump up to the next gufe release
        assert system.solvent_component.smiles == "O"


class TestThrombinSystem:
    """
    Tests to check that the thrombin system is properly created
    """
    @pytest.fixture()
    def system(self):
        return thrombin.get_system()

    def test_ligand_components(self, system):
        assert len(system.ligand_components) == 11

    def test_edges(self, system):
        assert len(system.ligand_network.edges) == 10

    def test_protein_component(self, system):
        assert system.protein_component.name == "thrombin"

    def test_solvent_component(self, system):
        # TODO use name once we bump up to the next gufe release
        assert system.solvent_component.smiles == "O"


class TestHIF2ASystem:
    """
    Tests to check that the hif2a system is properly created
    """
    @pytest.fixture()
    def system(self):
        return hif2a.get_system()

    def test_ligand_components(self, system):
        assert len(system.ligand_components) == 37

    def test_edges(self, system):
        assert len(system.ligand_network.edges) == 36

    def test_protein_component(self, system):
        assert system.protein_component.name == "hif2a"

    def test_solvent_component(self, system):
        # TODO use name once we bump up to the next gufe release
        assert system.solvent_component.smiles == "O"


class TestSykSystem:
    """
    Tests to check that the syk system is properly created
    """
    @pytest.fixture()
    def system(self):
        return syk.get_system()

    def test_ligand_components(self, system):
        assert len(system.ligand_components) == 43

    def test_edges(self, system):
        assert len(system.ligand_network.edges) == 42

    def test_protein_component(self, system):
        assert system.protein_component.name == "syk"

    def test_solvent_component(self, system):
        # TODO use name once we bump up to the next gufe release
        assert system.solvent_component.smiles == "O"
