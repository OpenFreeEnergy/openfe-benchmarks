import itertools
from importlib import resources
from typing import Iterable, List, Tuple

from rdkit import Chem
from openfe.setup.lomap_mapper import LomapAtomMapper
from openfe.setup import (
    Network, SmallMoleculeComponent, SolventComponent, ProteinComponent,
    LigandAtomMapper, LigandAtomMapping,
)
from openff.units import unit


def generate_relative_network_from_names(ligands: Iterable[SmallMoleculeComponent],
                                         connections: List[Tuple[str, str]],
                                         mappers: Iterable[LigandAtomMapper],
                                         scorer=None):
    """
    Generate a network from an input list of tuples each containing
    a pair of ligand names to connect with an edge.

    Parameters
    ----------
    ligands : Iterable[SmallMoleculeComponent]
        The ligands to create the network from.
    connections : List[Tuple[str, str]]
        The list of edges to create with each node identified by its
        SmallMoleculeComponent name.
    mapper : Iterable[LigandAtomMapper]
        Mappers to use. At least 1 is required.
    scorer : scoring function, optional
        A callable which returns a float for any LigandAtomMapping. Used to
        assign score to potential mappings, higher scores indicate worse
        mappings.

    Raises
    ------
    ValueError
        If no mapping can be found for a supplied edge.

    Returns
    -------
    network : Network
        Network of SmallMoleculeComponent transformations.
    """
    edges = []

    for entry in connections:
        nodes = [idx for idx, lig in enumerate(ligands) if lig.name in entry]

        for mapping in itertools.chain.from_iterable(
            mapper.suggest_mappings(ligands[nodes[0]], ligands[nodes[1]])
            for mapper in mappers
        ):
            if not scorer:
                best_mapping = mapping
                break

            score = scorer(mapping)
            mapping = mapping.with_annotations({"score": score})

            if score < best_score:
                best_mapping = mapping
                best_score = score

        if best_mapping is None:
            raise ValueError(f"No mapping for pair {entry}")
        edges.append(best_mapping)

    return Network(edges)


class RBFEBenchmarkSystem:
    """
    Class defining the components and alchemical network of a relative free
    energy benchmark system.

    Parameters
    ----------
    system_name : str
        The name of the benchmark system
    connections : List[Tuple[str, str]]
        The list of edges to create with each node identified by its
        SmallMoleculeComponent name.
    mappers : Iterable[LigandAtomMapper]
        Mappers to use. At least 1 is required.
    scorer : scoring function, optional
        A callable which returns a float for any LigandAtomMapping. Used to
        assign score to potential mappings, higher scores indicate worse
        mappings.

    Attributes
    ----------
    system_name : str
        The name / identifier of the benchmark system
    mappers : Iterable[LigandAtomMapper]
        Mappers used to create the ligand network
    scorer : Union[Callable, None]
        Scorer (if provided) used to create the ligand network
    ligand_components : List[SmallMoleculeComponent]
        List of SmallMoleculeComponent objects for each ligand in the benchmark
        system.
    ligand_network : Network
        Network of SmallMoleculeComponent transformations.
    protein_component : ProteinComponent
        ProteinComponent defining the host molecule of the benchmark system
    solvent_component : SolventComponent
        SolventComponent defining the solvent used for the benchmark system
    """
    def __init__(self, system_name: str, connections: List[Tuple[str, str]],
                 mappers=Iterable[LigandAtomMapper], scorer=None):
        self.system_name = system_name
        self.mappers = mappers
        self.scorer = scorer

        # Extract ligands
        self.ligand_components = self.extract_ligands(self.system_name)
        self.ligand_network = generate_relative_network_from_names(
            self.ligand_components, connections=connections,
            mappers=mappers, scorer=scorer)

        # Extract protein
        self.protein_component = self.extract_protein(self.system_name)

        # Create solvent component
        self.solvent_component = SolventComponent(
            positive_ion='Na', negative_ion='Cl',
            neutralize=True, ion_concentration=0.15*unit.molar)

    @staticmethod
    def extract_ligands(systemname: str):
        with resources.path('openfe_benchmarks.data',
                            f'{systemname}_ligands.sdf') as fn:
            ligands_sdf = Chem.SDMolSupplier(str(fn), removeHs=False)
        return [SmallMoleculeComponent(sdf) for sdf in ligands_sdf]

    @staticmethod
    def extract_protein(systemname: str):
        with resources.path('openfe_benchmarks.data',
                            f'{systemname}_protein.pdb') as fn:
            protein = ProteinComponent.from_pdbfile(str(fn), name=systemname)
        return protein
