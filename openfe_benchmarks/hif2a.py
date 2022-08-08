from rdkit import Chem
from openfe.setup.atom_mapping.lomap_mapper import LomapAtomMapper

from openfe_benchmarks.utils import RBFEBenchmarkSystem


def get_system(
    mappers=[LomapAtomMapper(time=20, threed=True, element_change=False, max3d=1),],
    scorer=None
):
    """
    Returns a RBFEBenchmarkSystem describing the hif2a benchmark system
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
        ligand network of the hif2a benchmark system
    """


    connections = [("lig_167", "lig_25"),
                   ("lig_235", "lig_25"),
                   ("lig_155", "lig_57"),
                   ("lig_163", "lig_67"),
                   ("lig_15", "lig_155"),
                   ("lig_15", "lig_163"),
                   ("lig_163", "lig_215"),
                   ("lig_231", "lig_289"),
                   ("lig_227", "lig_57"),
                   ("lig_1",   "lig_155"),
                   ("lig_251", "lig_256"),
                   ("lig_336", "lig_338"),
                   ("lig_163", "lig_231"),
                   ("lig_163", "lig_206"),
                   ("lig_42", "lig_43"),
                   ("lig_35", "lig_41"),
                   ("lig_35", "lig_50"),
                   ("lig_1", "lig_61"),
                   ("lig_124", "lig_54"),
                   ("lig_163", "lig_234"),
                   ("lig_252", "lig_256"),
                   ("lig_167", "lig_234"),
                   ("lig_252", "lig_338"),
                   ("lig_156", "lig_57"),
                   ("lig_61", "lig_84"),
                   ("lig_155", "lig_42"),
                   ("lig_163", "lig_165"),
                   ("lig_31", "lig_43"),
                   ("lig_224", "lig_235"),
                   ("lig_163", "lig_252"),
                   ("lig_227", "lig_290"),
                   ("lig_30", "lig_31"),
                   ("lig_252", "lig_254"),
                   ("lig_50", "lig_54"),
                   ("lig_1", "lig_23"),
                   ("lig_41", "lig_42"),]

    return RBFEBenchmarkSystem(system_name="hif2a", connections=connections,
                               mappers=mappers, scorer=scorer)
