from rdkit import Chem
from openfe.setup.atom_mapping.lomap_mapper import LomapAtomMapper

from openfe_benchmarks.utils import RBFEBenchmarkSystem


def get_system(
    mappers=[LomapAtomMapper(time=20, threed=True, element_change=False, max3d=1),],
    scorer=None
):
    """
    Returns a RBFEBenchmarkSystem describing the tnsk2 benchmark system
    with a network of ligand atom mappings defined by the input mappers
    and scorer.

    Parameters
    ----------
    mappers : Iterable[LigandAtomMapper]
        Mappers to use to create the ligand transformation network.
    scorer : Union[Callable, None]
        Scoring function for potential mappings, higher scores indicate
        worse mappings.

    Returns
    -------
    system : RBFEBenchmarkSystem
        RBFEBenchmarkSystem defining the various components and the
        ligand network of the tnsk2 benchmark system
    """

    connections = [("lig_1a", "lig_1b"),
                   ("lig_3a", "lig_5b"),
                   ("lig_1a", "lig_3a"),
                   ("lig_5f", "lig_5l"),
                   ("lig_1b", "lig_8b"),
                   ("lig_5a", "lig_5i"),
                   ("lig_5a", "lig_5g"),
                   ("lig_5d", "lig_5m"),
                   ("lig_8b", "lig_8d"),
                   ("lig_5d", "lig_5p"),
                   ("lig_5a", "lig_5c"),
                   ("lig_8c", "lig_8e"),
                   ("lig_8e", "lig_8f"),
                   ("lig_5d", "lig_5j"),
                   ("lig_5n", "lig_5o"),
                   ("lig_8a", "lig_8b"),
                   ("lig_3a", "lig_5a"),
                   ("lig_8d", "lig_8f"),
                   ("lig_1b", "lig_7"),
                   ("lig_5a", "lig_5f"),
                   ("lig_5a", "lig_5h"),
                   ("lig_5b", "lig_5e"),
                   ("lig_5f", "lig_5k"),
                   ("lig_3a", "lig_3b"),
                   ("lig_5m", "lig_5o"),
                   ("lig_5d", "lig_5k"),]

    return RBFEBenchmarkSystem(system_name="tnsk2", connections=connections,
                               mappers=mappers, scorer=scorer)
