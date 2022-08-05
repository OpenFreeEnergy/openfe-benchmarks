from rdkit import Chem
from openfe.setup.atom_mapping.lomap_mapper import LomapAtomMapper

from openfe_benchmarks.utils import RBFEBenchmarkSystem


def get_system(
    mappers=[LomapAtomMapper(time=20, threed=True, element_change=False, max3d=1),],
    scorer=None
):
    """
    Returns a RBFEBenchmarkSystem describing the MCL1 benchmark system
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
        ligand network of the MCL1 benchmark system
    """
    connections = [("lig_27", "lig_48"),
                   ("lig_27", "lig_46"),
                   ("lig_27", "lig_30"),
                   ("lig_37", "lig_60"),
                   ("lig_28", "lig_33"),
                   ("lig_30", "lig_35"),
                   ("lig_35", "lig_36"),
                   ("lig_35", "lig_50"),
                   ("lig_27", "lig_47"),
                   ("lig_37", "lig_67"),
                   ("lig_27", "lig_32"),
                   ("lig_60", "lig_65"),
                   ("lig_30", "lig_31"),
                   ("lig_35", "lig_52"),
                   ("lig_60", "lig_61"),
                   ("lig_32", "lig_34"),
                   ("lig_49", "lig_52"),
                   ("lig_33", "lig_52"),
                   ("lig_37", "lig_52"),
                   ("lig_60", "lig_63"),
                   ("lig_37", "lig_56"),
                   ("lig_53", "lig_58"),
                   ("lig_35", "lig_53"),]

    return RBFEBenchmarkSystem(system_name="mcl1", connections=connections,
                               mappers=mappers, scorer=scorer)
