from rdkit import Chem
from openfe.setup.atom_mapping.lomap_mapper import LomapAtomMapper

from openfe_benchmarks.utils import RHFEBenchmarkSystem


def get_system(
    mappers=[LomapAtomMapper(time=20, threed=True, element_change=False, max3d=1),],
    scorer=None
):
    """
    Returns a RBFEBenchmarkSystem describing the CMET benchmark system
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
        ligand network of the CMET benchmark system
    """

    connections = [('lig_1', 'lig_12'),
                    ('lig_2', 'lig_12'),
                    ('lig_3', 'lig_12'),
                    ('lig_4', 'lig_12'),
                    ('lig_5', 'lig_12'),
                    ('lig_6', 'lig_12'),
                    ('lig_7', 'lig_12'),
                    ('lig_8', 'lig_12'),
                    ('lig_9', 'lig_12'),
                    ('lig_10', 'lig_12'),
                    ('lig_11', 'lig_12'),
                    ('lig_13', 'lig_12'),
                    ('lig_14', 'lig_12'),
                    ('lig_15', 'lig_12'),
                    ('lig_16', 'lig_12')]

    return RHFEBenchmarkSystem(system_name="benzenes", connections=connections,
                               mappers=mappers, scorer=scorer)
