from rdkit import Chem
from openfe.setup.atom_mapping.lomap_mapper import LomapAtomMapper

from openfe_benchmarks.utils import RBFEBenchmarkSystem


def get_system(
    mappers=[LomapAtomMapper(time=20, threed=True, element_change=False, max3d=1),],
    scorer=None
):
    """
    Returns a RBFEBenchmarkSystem describing the syk benchmark system
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
        ligand network of the syk benchmark system
    """

    connections = [("lig_CHEMBL3264994", "lig_CHEMBL3264995"),
                   ("lig_CHEMBL3265016", "lig_CHEMBL3265019"),
                   ("lig_CHEMBL3264999", "lig_CHEMBL3264996"),
                   ("lig_CHEMBL3265020", "lig_CHEMBL3265021"),
                   ("lig_CHEMBL3265025", "lig_CHEMBL3265026"),
                   ("lig_CHEMBL3265005", "lig_CHEMBL3265019"),
                   ("lig_CHEMBL3265015", "lig_CHEMBL3259820"),
                   ("lig_CHEMBL3265003", "lig_CHEMBL3264999"),
                   ("lig_CHEMBL3265023", "lig_CHEMBL3265005"),
                   ("lig_CHEMBL3265035", "lig_CHEMBL3264999"),
                   ("lig_CHEMBL3265013", "lig_CHEMBL3259820"),
                   ("lig_CHEMBL3265005", "lig_CHEMBL3265008"),
                   ("lig_CHEMBL3265020", "lig_CHEMBL3265018"),
                   ("lig_CHEMBL3265005", "lig_CHEMBL3265011"),
                   ("lig_CHEMBL3264995", "lig_CHEMBL3264998"),
                   ("lig_CHEMBL3265005", "lig_CHEMBL3259820"),
                   ("lig_CHEMBL3265037", "lig_CHEMBL3265032"),
                   ("lig_CHEMBL3265027", "lig_CHEMBL3265031"),
                   ("lig_CHEMBL3265028", "lig_CHEMBL3265027"),
                   ("lig_CHEMBL3264995", "lig_CHEMBL3264997"),
                   ("lig_CHEMBL3265024", "lig_CHEMBL3265026"),
                   ("lig_CHEMBL3265027", "lig_CHEMBL3265033"),
                   ("lig_CHEMBL3265012", "lig_CHEMBL3259820"),
                   ("lig_CHEMBL3265017", "lig_CHEMBL3265029"),
                   ("lig_CHEMBL3265010", "lig_CHEMBL3265005"),
                   ("lig_CHEMBL3265027", "lig_CHEMBL3265036"),
                   ("lig_CHEMBL3265027", "lig_CHEMBL3265034"),
                   ("lig_CHEMBL3265002", "lig_CHEMBL3265004"),
                   ("lig_CHEMBL3265020", "lig_CHEMBL3265019"),
                   ("lig_CHEMBL3265009", "lig_CHEMBL3265005"),
                   ("lig_CHEMBL3264999", "lig_CHEMBL3259820"),
                   ("lig_CHEMBL3265030_n", "lig_CHEMBL3265035"),
                   ("lig_CHEMBL3264999", "lig_CHEMBL3265004"),
                   ("lig_CHEMBL3265020", "lig_CHEMBL3265017"),
                   ("lig_CHEMBL3265027", "lig_CHEMBL3265029"),
                   ("lig_CHEMBL3265005", "lig_CHEMBL3265026"),
                   ("lig_CHEMBL3265014", "lig_CHEMBL3259820"),
                   ("lig_CHEMBL3265027", "lig_CHEMBL3265032"),
                   ("lig_CHEMBL3265017", "lig_CHEMBL3265022"),
                   ("lig_CHEMBL3265006", "lig_CHEMBL3265005"),
                   ("lig_CHEMBL3264999", "lig_CHEMBL3265001"),
                   ("lig_CHEMBL3264995", "lig_CHEMBL3265001"),]

    return RBFEBenchmarkSystem(system_name="syk", connections=connections,
                               mappers=mappers, scorer=scorer)
