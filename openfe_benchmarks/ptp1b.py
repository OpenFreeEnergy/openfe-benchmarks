from rdkit import Chem
from openfe.setup.atom_mapping.lomap_mapper import LomapAtomMapper

from openfe_benchmarks.utils import RBFEBenchmarkSystem


def get_system(mappers=[LomapAtomMapper(threed=True),], scorer=None):
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

    connections = [('lig_20667_2qbp', 'lig_23479'),
                   ('lig_20667_2qbp', 'lig_23482'),
                   ('lig_20667_2qbp', 'lig_23484'),
                   ('lig_20667_2qbp', 'lig_23485'),
                   ('lig_20667_2qbp', 'lig_23486'),
                   ('lig_20669_2qbr', 'lig_23466'),
                   ('lig_20669_2qbr', 'lig_23472'),
                   ('lig_20670_2qbs', 'lig_23330_2qbq'),
                   ('lig_20670_2qbs', 'lig_23466'),
                   ('lig_20670_2qbs', 'lig_23477'),
                   ('lig_20670_2qbs', 'lig_23479'),
                   ('lig_20670_2qbs', 'lig_23482'),
                   ('lig_20670_2qbs', 'lig_23483'),
                   ('lig_23466', 'lig_23475'),
                   ('lig_23467', 'lig_23466'),
                   ('lig_23467', 'lig_23468'),
                   ('lig_23467', 'lig_23469'),
                   ('lig_23467', 'lig_23470'),
                   ('lig_23467',  'lig_23473'),
                   ('lig_23467', 'lig_23474'),
                   ('lig_23467', 'lig_23475'),
                   ('lig_23467', 'lig_23476'),
                   ('lig_23469', 'lig_20669_2qbr'),
                   ('lig_23469', 'lig_23472'),
                   ('lig_23471', 'lig_23466'),
                   ('lig_23471', 'lig_23468'),
                   ('lig_23471', 'lig_23470'),
                   ('lig_23473', 'lig_20669_2qbr'),
                   ('lig_23474', 'lig_23466'),
                   ('lig_23476', 'lig_23466'),
                   ('lig_23477', 'lig_23330_2qbq'),
                   ('lig_23477', 'lig_23466'),
                   ('lig_23477', 'lig_23467'),
                   ('lig_23477', 'lig_23479'),
                   ('lig_23477',  'lig_23482'),
                   ('lig_23477', 'lig_23483'),
                   ('lig_23480', 'lig_23479'),
                   ('lig_23480', 'lig_23482'),
                   ('lig_23482', 'lig_23479'),
                   ('lig_23482', 'lig_23485'),
                   ('lig_23482', 'lig_23486'),
                   ('lig_23483', 'lig_23479'),
                   ('lig_23484', 'lig_23479'),
                   ('lig_23484', 'lig_23482'),
                   ('lig_23484', 'lig_23485'),
                   ('lig_23484', 'lig_23486'),
                   ('lig_23485', 'lig_23479'),
                   ('lig_23486', 'lig_23479'),
                   ('lig_23486', 'lig_23485'),]

    return RBFEBenchmarkSystem(system_name="ptp1b", connections=connections,
                               mappers=mappers, scorer=scorer)
