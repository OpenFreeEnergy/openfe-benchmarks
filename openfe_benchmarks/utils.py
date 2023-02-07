import itertools
from importlib import resources
from typing import Iterable, List, Tuple

from rdkit import Chem
from rdkit.Geometry.rdGeometry import Point3D
from openfe.setup.atom_mapping.lomap_mapper import LomapAtomMapper
from openfe.setup.atom_mapping import LigandAtomMapper, LigandAtomMapping
from openfe.setup import (
    LigandNetwork, SmallMoleculeComponent, SolventComponent, ProteinComponent,
)
from openff.units import unit

import py3Dmol
from matplotlib import pyplot as plt
from matplotlib.colors import rgb2hex


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
    network : LigandNetwork
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

    return LigandNetwork(edges)


class RHFEBenchmarkSystem:
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
    ligand_network : LigandNetwork
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
        self.connections = connections

        # Extract ligands
        self.ligand_components = self.extract_ligands(self.system_name)
        self.ligand_network = generate_relative_network_from_names(
            self.ligand_components, connections=connections,
            mappers=mappers, scorer=scorer)

        # Create solvent component
        self.solvent_component = SolventComponent(
            positive_ion='Na', negative_ion='Cl',
            neutralize=True, ion_concentration=0.15*unit.molar)

    @staticmethod
    def extract_ligands(systemname: str):
        with resources.path('openfe_benchmarks.data',
                            f'{systemname}_RHFE.sdf') as fn:
            ligands_sdf = Chem.SDMolSupplier(str(fn), removeHs=False)
        return [SmallMoleculeComponent(sdf) for sdf in ligands_sdf]
    
    
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
    ligand_network : LigandNetwork
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
        self.connections = connections

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
            protein = ProteinComponent.from_pdb_file(str(fn), name=systemname)
        return protein


def show_edge_3D(edge: LigandAtomMapping, spheres: bool=True,
                 style: str='stick',
                 shift: Tuple[float, float, float]=(10, 0, 0)):
    """
    Render relative transformation edge in 3D using py3Dmol.

    By default matching atoms will be annotated using colored spheres.

    Parameters
    ----------
    edge : LigandAtomMapping
        The ligand transformation edge to visualize.
    spheres : bool
        Whether or not to show matching atoms as spheres.
    style : str
        Style in which to represent the molecules in py3Dmol.
    shift : Tuple of floats
        Amount to shift molB by in order to visualise the two ligands.
        By default molB is shifted by 10 A in the z direction.

    Returns
    -------
    view : py3Dmol.view
        View of the system containing both molecules in the edge.
    """
    def translate(mol, shift):
        conf = mol.GetConformer()
        for i, atom in enumerate(mol.GetAtoms()):
            x, y, z = conf.GetAtomPosition(i)
            point = Point3D(x+shift[0], y+shift[1], z+shift[2])
            conf.SetAtomPosition(i, point)
        return mol

    def add_spheres(view, mol1, mol2, mapping):
        # Get colourmap of size mapping
        cmap = plt.cm.get_cmap('hsv', len(mapping))
        for i, pair in enumerate(mapping.items()):
            p1 = mol1.GetConformer().GetAtomPosition(pair[0])
            p2 = mol2.GetConformer().GetAtomPosition(pair[1])
            color = rgb2hex(cmap(i))
            view.addSphere({"center": {"x":p1.x, "y": p1.y, "z": p1.z},
                            "radius": 0.6, "color": color, "alpha": 0.8})
            view.addSphere({"center": {"x":p2.x, "y": p2.y, "z": p2.z},
                            "radius": 0.6, "color": color, "alpha": 0.8})

    molA = edge.componentA.to_rdkit()
    molB = edge.componentB.to_rdkit()

    mblock1 = Chem.MolToMolBlock(molA)
    mblock2 = Chem.MolToMolBlock(translate(molB, shift))

    view = py3Dmol.view(width=600, height=600)
    view.addModel(mblock1, 'molA')
    view.addModel(mblock2, 'molB')

    if spheres:
        add_spheres(view, molA, molB, edge.componentA_to_componentB)

    view.setStyle({style:{}})
    view.zoomTo()
    return view
