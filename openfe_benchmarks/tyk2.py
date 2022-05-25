from rdkit import Chem
from openfe.setup.lomap_mapper import LomapAtomMapper

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

    connections = [('ligand_23', 'ligand_46'),
                   ('ligand_23', 'ligand_55'),
                   ('ligand_23', 'ligand_27'),
                   ('ligand_23', 'ligand_30'),
                   ('ligand_28', 'ligand_27'),
                   ('ligand_28', 'ligand_30'),
                   ('ligand_31', 'ligand_43'),
                   ('ligand_31', 'ligand_45'),
                   ('ligand_31', 'ligand_46'),
                   ('ligand_31', 'ligand_48'),
                   ('ligand_31', 'ligand_28'),
                   ('ligand_42', 'ligand_48'),
                   ('ligand_42', 'ligand_54'),
                   ('ligand_42', 'ligand_55'),
                   ('ligand_43', 'ligand_55'),
                   ('ligand_44', 'ligand_42'),
                   ('ligand_44', 'ligand_55'),
                   ('ligand_45', 'ligand_42'),
                   ('ligand_47', 'ligand_31'),
                   ('ligand_47', 'ligand_55'),
                   ('ligand_49', 'ligand_31'),
                   ('ligand_49', 'ligand_50'),
                   ('ligand_50', 'ligand_42'),
                   ('ligand_55', 'ligand_54'),]

    return RBFEBenchmarkSystem(system_name="tyk2", connections=connections,
                               mappers=mappers, scorer=scorer)
