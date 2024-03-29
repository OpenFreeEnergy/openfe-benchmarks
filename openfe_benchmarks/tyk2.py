from rdkit import Chem
from openfe.setup.atom_mapping.lomap_mapper import LomapAtomMapper

from openfe_benchmarks.utils import RBFEBenchmarkSystem


def get_system(
    mappers=[LomapAtomMapper(time=20, threed=True, element_change=False, max3d=1),],
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
        ligand network of the TYK2 benchmark system
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


def get_old_system(
    mappers=[LomapAtomMapper(time=20, threed=True, element_change=False, max3d=1),],
    scorer=None
):
    """
    Returns a RBFEBenchmarkSystem describing the "old" TYK2 benchmark system
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
        ligand network of the TYK2 benchmark system
    """

    connections = [("lig_ejm_31", "lig_ejm_43"),
                   ("lig_ejm_31", "lig_ejm_45"),
                   ("lig_ejm_31", "lig_ejm_46"),
                   ("lig_ejm_31", "lig_ejm_48"),
                   ("lig_ejm_31", "lig_jmc_28"),
                   ("lig_ejm_42", "lig_ejm_48"),
                   ("lig_ejm_42", "lig_ejm_54"),
                   ("lig_ejm_42", "lig_ejm_55"),
                   ("lig_ejm_43", "lig_ejm_55"),
                   ("lig_ejm_44", "lig_ejm_42"),
                   ("lig_ejm_44", "lig_ejm_55"),
                   ("lig_ejm_45", "lig_ejm_42"),
                   ("lig_ejm_47", "lig_ejm_31"),
                   ("lig_ejm_47", "lig_ejm_55"),
                   ("lig_ejm_49", "lig_ejm_31"),
                   ("lig_ejm_49", "lig_ejm_50"),
                   ("lig_ejm_50", "lig_ejm_42"),
                   ("lig_ejm_55", "lig_ejm_54"),
                   ("lig_jmc_23", "lig_ejm_46"),
                   ("lig_jmc_23", "lig_ejm_55"),
                   ("lig_jmc_23", "lig_jmc_27"),
                   ("lig_jmc_23", "lig_jmc_30"),
                   ("lig_jmc_28", "lig_jmc_27"),
                   ("lig_jmc_28", "lig_jmc_30"),]

    return RBFEBenchmarkSystem(system_name="tyk2_old", connections=connections,
                               mappers=mappers, scorer=scorer)


def get_new_system(
    mappers=[LomapAtomMapper(time=20, threed=True, element_change=False, max3d=1),],
    scorer=None
):
    """
    Returns a RBFEBenchmarkSystem describing the "new" TYK2 benchmark system
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
        ligand network of the TYK2 benchmark system
    """

    connections = [("lig_ejm_45", "lig_ejm_48"),
                   ("lig_ejm_48", "lig_ejm_43"),
                   ("lig_jmc_30", "lig_ejm_48"),
                   ("lig_ejm_46", "lig_ejm_48"),
                   ("lig_ejm_31", "lig_ejm_48"),
                   ("lig_ejm_50", "lig_ejm_48"),
                   ("lig_ejm_54", "lig_ejm_48"),
                   ("lig_jmc_27", "lig_ejm_48"),
                   ("lig_jmc_23", "lig_ejm_48"),
                   ("lig_ejm_55", "lig_ejm_48"),
                   ("lig_ejm_47", "lig_ejm_48"),
                   ("lig_ejm_42", "lig_ejm_48"),
                   ("lig_jmc_28", "lig_ejm_48"),]

    return RBFEBenchmarkSystem(system_name="tyk2_new", connections=connections,
                               mappers=mappers, scorer=scorer)
