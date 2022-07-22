from rdkit import Chem
from openfe.setup.atom_mapping.lomap_mapper import LomapAtomMapper

from openfe_benchmarks.utils import RBFEBenchmarkSystem


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

    connections = [
            ("lig_CHEMBL3402744_300_4", "lig_CHEMBL3402745_200_5"),
            ("lig_CHEMBL3402745_200_5", "lig_CHEMBL3402749_500_9"),
            ("lig_CHEMBL3402745_200_5", "lig_CHEMBL3402754_40_14"),
            ("lig_CHEMBL3402754_40_14", "lig_CHEMBL3402761_1_21"),]

    return RBFEBenchmarkSystem(system_name="cmet", connections=connections,
                               mappers=mappers, scorer=scorer)
