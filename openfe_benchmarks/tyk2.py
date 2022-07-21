from rdkit import Chem
from openfe.setup.lomap_mapper import LomapAtomMapper

from openfe_benchmarks.utils import RBFEBenchmarkSystem


def get_system(
    mappers=[LomapAtomMapper(time=20, threed=True, element_changes=False, max3d=1),],
    scorer=None
):
    """
    Returns a RBFEBenchmarkSystem describing the TYK2 benchmark system
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
        ligand network of the TY2K benchmark system
    """

    connections = [("lig_ejm_31", "lig_ejm_50"),
                   ("lig_ejm_46", "lig_jmc_23"),
                   ("lig_ejm_31", "lig_ejm_55"),
                   ("lig_ejm_31", "lig_ejm_48"),
                   ("lig_ejm_31", "lig_ejm_54"),
                   ("lig_ejm_31", "lig_ejm_47"),
                   ("lig_ejm_31", "lig_ejm_46"),
                   ("lig_ejm_46", "lig_jmc_27"),
                   ("lig_ejm_46", "lig_jmc_28"),
                   ("lig_ejm_42", "lig_ejm_43"),
                   ("lig_ejm_31", "lig_ejm_42"),
                   ("lig_ejm_45", "lig_ejm_55"),]

    return RBFEBenchmarkSystem(system_name="tyk2", connections=connections,
                               mappers=mappers, scorer=scorer)
