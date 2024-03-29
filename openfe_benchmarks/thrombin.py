from rdkit import Chem
from openfe.setup.atom_mapping.lomap_mapper import LomapAtomMapper

from openfe_benchmarks.utils import RBFEBenchmarkSystem


def get_system(
    mappers=[LomapAtomMapper(time=20, threed=True, element_change=False, max3d=1),],
    scorer=None
):
    """
    Returns a RBFEBenchmarkSystem describing the thrombin benchmark system
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
        ligand network of the thrombin benchmark system
    """

    connections = [("lig_1d", "lig_5"),
                   ("lig_1b", "lig_6b"),
                   ("lig_1b", "lig_6e"),
                   ("lig_1b", "lig_6a"),
                   ("lig_3a", "lig_3b"),
                   ("lig_6b", "lig_7a"),
                   ("lig_1b", "lig_5"),
                   ("lig_1a", "lig_5"),
                   ("lig_1c", "lig_5"),
                   ("lig_3a", "lig_5"),]

    return RBFEBenchmarkSystem(system_name="thrombin", connections=connections,
                               mappers=mappers, scorer=scorer)
