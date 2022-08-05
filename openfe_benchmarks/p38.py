from rdkit import Chem
from openfe.setup.atom_mapping.lomap_mapper import LomapAtomMapper

from openfe_benchmarks.utils import RBFEBenchmarkSystem


def get_system(
    mappers=[LomapAtomMapper(time=20, threed=True, element_change=False, max3d=1),],
    scorer=None
):
    """
    Returns a RBFEBenchmarkSystem describing the p38 benchmark system
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
        ligand network of the p38 benchmark system
    """

    connections = [("lig_p38a_2f", "lig_p38a_2i"),
                   ("lig_p38a_2l", "lig_p38a_3fln"),
                   ("lig_p38a_2f", "lig_p38a_3flz"),
                   ("lig_p38a_2g", "lig_p38a_3flz"),
                   ("lig_p38a_2l", "lig_p38a_2o"),
                   ("lig_p38a_2bb", "lig_p38a_2v"),
                   ("lig_p38a_2e", "lig_p38a_3fln"),
                   ("lig_p38a_2i", "lig_p38a_2j"),
                   ("lig_p38a_2f", "lig_p38a_2t"),
                   ("lig_p38a_2z", "lig_p38a_3fly"),
                   ("lig_p38a_2g", "lig_p38a_3fln"),
                   ("lig_p38a_2n", "lig_p38a_2t"),
                   ("lig_p38a_2y", "lig_p38a_3fly"),
                   ("lig_p38a_3flw", "lig_p38a_3fly"),
                   ("lig_p38a_2v", "lig_p38a_3fly"),
                   ("lig_p38a_2m", "lig_p38a_2t"),
                   ("lig_p38a_2v", "lig_p38a_3fmh"),
                   ("lig_p38a_2v", "lig_p38a_2x"),
                   ("lig_p38a_2l", "lig_p38a_2r"),
                   ("lig_p38a_2gg", "lig_p38a_3fln"),
                   ("lig_p38a_2l", "lig_p38a_2p"),
                   ("lig_p38a_2n", "lig_p38a_2s"),
                   ("lig_p38a_2v", "lig_p38a_3fln"),
                   ("lig_p38a_2ff", "lig_p38a_3fln"),
                   ("lig_p38a_2aa", "lig_p38a_2z"),
                   ("lig_p38a_2h", "lig_p38a_3flz"),
                   ("lig_p38a_2ee", "lig_p38a_3fln"),
                   ("lig_p38a_2k", "lig_p38a_3fln"),]

    return RBFEBenchmarkSystem(system_name="p38", connections=connections,
                               mappers=mappers, scorer=scorer)
