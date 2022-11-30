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


def get_old_system(
    mappers=[LomapAtomMapper(time=20, threed=True, element_change=False, max3d=1),],
    scorer=None,
):
    """
    Returns a RBFEBenchmarkSystem describing the "old" mcl1 benchmark system
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
        ligand network of the mcl1 benchmark system
    """

    connections = [("lig_26", "lig_44"),
                   ("lig_26", "lig_57"),
                   ("lig_26", "lig_64"),
                   ("lig_27", "lig_23"),
                   ("lig_27", "lig_45"),
                   ("lig_27", "lig_46"),
                   ("lig_28", "lig_27"),
                   ("lig_28", "lig_35"),
                   ("lig_28", "lig_47"),
                   ("lig_29", "lig_27"),
                   ("lig_29", "lig_35"),
                   ("lig_29", "lig_40"),
                   ("lig_30", "lig_27"),
                   ("lig_30", "lig_35"),
                   ("lig_30", "lig_40"),
                   ("lig_30", "lig_48"),
                   ("lig_31", "lig_35"),
                   ("lig_32", "lig_34"),
                   ("lig_32", "lig_46"),
                   ("lig_33", "lig_27"),
                   ("lig_35", "lig_33"),
                   ("lig_35", "lig_34"),
                   ("lig_35", "lig_36"),
                   ("lig_35", "lig_37"),
                   ("lig_35", "lig_39"),
                   ("lig_35", "lig_53"),
                   ("lig_35", "lig_60"),
                   ("lig_38", "lig_35"),
                   ("lig_38", "lig_60"),
                   ("lig_39", "lig_32"),
                   ("lig_41", "lig_32"),
                   ("lig_41", "lig_35"),
                   ("lig_42", "lig_51"),
                   ("lig_42", "lig_64"),
                   ("lig_43", "lig_27"),
                   ("lig_43", "lig_47"),
                   ("lig_44", "lig_23"),
                   ("lig_48", "lig_27"),
                   ("lig_49", "lig_35"),
                   ("lig_49", "lig_67"),
                   ("lig_50", "lig_60"),
                   ("lig_51", "lig_45"),
                   ("lig_52", "lig_60"),
                   ("lig_54", "lig_23"),
                   ("lig_54", "lig_42"),
                   ("lig_56", "lig_35"),
                   ("lig_56", "lig_60"),
                   ("lig_57", "lig_23"),
                   ("lig_58", "lig_60"),
                   ("lig_60", "lig_36"),
                   ("lig_61", "lig_60"),
                   ("lig_62", "lig_26"),
                   ("lig_62", "lig_45"),
                   ("lig_63", "lig_60"),
                   ("lig_65", "lig_60"),
                   ("lig_65", "lig_67"),
                   ("lig_66", "lig_23"),
                   ("lig_66", "lig_42"),
                   ("lig_67", "lig_27"),
                   ("lig_67", "lig_31"),
                   ("lig_67", "lig_32"),
                   ("lig_67", "lig_35"),
                   ("lig_67", "lig_37"),
                   ("lig_67", "lig_50"),
                   ("lig_67", "lig_52"),
                   ("lig_67", "lig_53"),
                   ("lig_67", "lig_58"),
                   ("lig_67", "lig_61"),
                   ("lig_67", "lig_63"),
                   ("lig_68", "lig_23"),
                   ("lig_68", "lig_45"),]

    return RBFEBenchmarkSystem(system_name="mcl1_old", connections=connections,
                               mappers=mappers, scorer=scorer)
