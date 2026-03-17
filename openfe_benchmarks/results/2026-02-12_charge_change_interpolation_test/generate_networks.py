"""
Create alchemical networks for all of the charge change transformations in the benchmarks sets.
"""
import logging
from openfe_benchmarks.data import BenchmarkIndex, get_benchmark_set_data_systems
from openfe_benchmarks.scripts import utils as ofebu
from openfe import SolventComponent, ProteinComponent
from openfe.protocols.openmm_rfe.equil_rfe_methods import RelativeHybridTopologyProtocol
from gufe import LigandAtomMapping, ChemicalSystem, Transformation, AlchemicalNetwork, LigandNetwork
import pathlib


logger = logging.getLogger(__name__)

LIG_NETWORK_FILE = "industry_benchmarks_network"
PARTIAL_CHARGE = "nagl_openff-gnn-am1bcc-1.0.0.pt"
FORCEFIELD = "openff-2.3.0"
OUTPUT_DIR = pathlib.Path("networks")

def _configure_example_logging(level=logging.INFO):
    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(name)s %(levelname)s: %(message)s")
    )
    # Attach to this module's logger so output appears when running the example
    logger.addHandler(handler)
    logger.setLevel(level)
    # Optionally enable package-wide logs:
    logging.getLogger("openfe_benchmarks").setLevel(level)


def main():
    """Create alchemical networks using all charge change transformations available in the bfe sets"""
    _configure_example_logging()
    logger.info(f"Creating charge change transformations networks using {PARTIAL_CHARGE} partial charges and {FORCEFIELD} small molecule force field")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    bench_index = BenchmarkIndex()

    # default solvent component used for all runs
    solvent = SolventComponent()
    # set up the default settings for the protocol, these will be adapted for each edge as needed
    default_settings = RelativeHybridTopologyProtocol.default_settings()
    default_settings.forcefield_settings.small_molecule_forcefield = FORCEFIELD
    default_settings.alchemical_settings.turn_off_core_unique_exceptions = True

    # find all charge change edges across all benchmark sets
    for system_group in bench_index.list_benchmark_sets():
        for system_name, system_data in get_benchmark_set_data_systems(system_group).items():
            # get the network for this system
            network = LigandNetwork.from_json(system_data.ligand_networks[LIG_NETWORK_FILE])
            # get the ligands cofactors and protein components
            ligands_by_name = ofebu.process_sdf(system_data.ligands[PARTIAL_CHARGE], return_dict=True)
            protein = ProteinComponent.from_pdb_file(str(system_data.protein))
            cofactors = None
            if system_data.cofactors is not None:
                cofactors = ofebu.process_sdf(
                    system_data.ligands[PARTIAL_CHARGE], return_dict=False
                )

            transformations = []
            # check if this is a charge change edge
            for edge in network.edges:
                if edge.get_alchemical_charge_difference() != 0.0:
                    logger.info(f"Found charge change edge in {system_group} {system_name}: {edge.componentA.name} -> {edge.componentB.name}")
                    # process this edge
                    # add the system group and system name to the annotations
                    annotations = edge.annotations
                    annotations["system_group"] = system_group
                    annotations["system_name"] = system_name
                    new_edge = LigandAtomMapping(
                        componentA=ligands_by_name[edge.componentA.name],
                        componentB=ligands_by_name[edge.componentB.name],
                        componentA_to_componentB=edge.componentA_to_componentB,
                        annotations=annotations,
                    )
                    # create the transformations for the bound and solvent legs
                    for leg in ["solvent", "complex"]:
                        system_a_dict = {"ligand": new_edge.componentA, "solvent": solvent}
                        system_b_dict = {"ligand": new_edge.componentB, "solvent": solvent}
                        if leg == "complex":
                            system_a_dict["protein"] = protein
                            system_b_dict["protein"] = protein

                            if cofactors is not None:
                                for i, cofactor in enumerate(cofactors):
                                    cofactor_name = f"cofactor_{i}"
                                    system_a_dict[cofactor_name] = cofactor
                                    system_b_dict[cofactor_name] = cofactor

                        system_a = ChemicalSystem(system_a_dict)
                        system_b = ChemicalSystem(system_b_dict)

                        name = f"{leg}_{new_edge.componentA.name}_{new_edge.componentB.name}"

                        transformation_protocol = RelativeHybridTopologyProtocol(
                            settings=RelativeHybridTopologyProtocol._adaptive_settings(
                                stateA=system_a,
                                stateB=system_b,
                                mapping=new_edge,
                                initial_settings=default_settings,
                            ),
                        )

                        # Create transformation
                        transformation = Transformation(
                            stateA=system_a,
                            stateB=system_b,
                            mapping=new_edge,
                            protocol=transformation_protocol,
                            name=name,
                        )
                        transformations.append(transformation)

            if transformations:
                logger.info(f"Found {len(transformations) / 2} charge change transformations for system {system_group} {system_name}")
                # create a network with the transformations
                network = AlchemicalNetwork(edges=transformations, name=f"{system_group}_{system_name}")
                # save the network to disk
                output_name = OUTPUT_DIR / f"{system_group}_{system_name}_alchemicalnetwork.json"
                network.to_json(output_name)
                logger.info(f"Saved network to {output_name}")
            else:
                logger.info(f"No charge change transformations found for system {system_group} {system_name}")

if __name__ == "__main__":
    main()