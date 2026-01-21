import click
import pathlib
import logging
import warnings
import json
from openff.units import unit
import openfe
from openfe.protocols.openmm_rfe.equil_rfe_methods import RelativeHybridTopologyProtocol
from rdkit import Chem
from openff.toolkit import (
    RDKitToolkitWrapper, AmberToolsToolkitWrapper
)
from openff.toolkit.utils.toolkit_registry import ToolkitRegistry
import gufe
from gufe import tokenization


logger = logging.getLogger(__name__)
warnings.filterwarnings("ignore", message="Partial charges have been provided, these will preferentially be used instead of generating new partial charges")
amber_rdkit = ToolkitRegistry([RDKitToolkitWrapper(), AmberToolsToolkitWrapper()])


def get_alchemical_charge_difference(mapping) -> int:
    """
    Checks and returns the difference in formal charge between state A and B.

    Parameters
    ----------
    mapping : dict[str, ComponentMapping]
      Dictionary of mappings between transforming components.

    Returns
    -------
    int
      The formal charge difference between states A and B.
      This is defined as sum(charge state A) - sum(charge state B)
    """
    chg_A = Chem.rdmolops.GetFormalCharge(
        mapping.componentA.to_rdkit()
    )
    chg_B = Chem.rdmolops.GetFormalCharge(
        mapping.componentB.to_rdkit()
    )

    return chg_A - chg_B


def get_complex_settings(forcefield: str):
    """
    Utility method for getting RFEProtocol settings for non charge changing
    transformations.
    """
    # Are there additional settings we should specify here?
    settings = RelativeHybridTopologyProtocol.default_settings()
    settings.engine_settings.compute_platform = "CUDA"
    # Set ligand FF
    settings.forcefield_settings.small_molecule_forcefield = forcefield
    # Fast settings
    settings.simulation_settings.time_per_iteration = 2.5 * unit.picosecond
    settings.simulation_settings.real_time_analysis_interval = 1 * unit.nanosecond
    settings.output_settings.checkpoint_interval = 1 * unit.nanosecond
    settings.solvation_settings.box_shape = 'dodecahedron'
    settings.forcefield_settings.nonbonded_cutoff = 0.9 * unit.nanometer
    settings.solvation_settings.solvent_padding = 1.0 * unit.nanometer
    # Only run one repeat per input json file
    settings.protocol_repeats = 1
    return settings


def get_ligand_settings(forcefield: str):
    """
    Utility method for getting RFEProtocol settings for non charge changing
    transformations.
    """
    # Are there additional settings we should specify here?
    settings = RelativeHybridTopologyProtocol.default_settings()
    settings.engine_settings.compute_platform = "CUDA"
    # Set ligand FF
    settings.forcefield_settings.small_molecule_forcefield = forcefield
    # Fast settings
    settings.simulation_settings.time_per_iteration = 2.5 * unit.picosecond
    settings.simulation_settings.real_time_analysis_interval = 1 * unit.nanosecond
    settings.output_settings.checkpoint_interval = 1 * unit.nanosecond
    settings.solvation_settings.box_shape = 'dodecahedron'
    settings.forcefield_settings.nonbonded_cutoff = 0.9 * unit.nanometer
    settings.solvation_settings.solvent_padding = 1.5 * unit.nanometer
    # Only run one repeat per input json file
    settings.protocol_repeats = 1
    return settings


def get_settings_charge_changes(forcefield: str):
    """
    Utility method for getting RFEProtocol settings for charge changing
    transformations.

    These settings mostly follow defaults but use longer
    simulation times, more lambda windows and an alchemical ion.
    """
    settings = RelativeHybridTopologyProtocol.default_settings()
    settings.engine_settings.compute_platform = "CUDA"
    # Should we use this new OpenFF version or the default?
    settings.forcefield_settings.small_molecule_forcefield = forcefield
    settings.alchemical_settings.explicit_charge_correction = True
    settings.simulation_settings.production_length = 20 * unit.nanosecond
    settings.simulation_settings.n_replicas = 22
    settings.lambda_settings.lambda_windows = 22
    # Only run one repeat per input json file
    settings.protocol_repeats = 1
    return settings


def create_and_store_network(forcefield: str, ligand_network, solv, prot):

    # Create the AlchemicalTransformations, and storing them to an AlchemicalNetwork
    transformations = []
    for mapping in ligand_network.edges:
        # Get different settings depending on whether the transformation
        # involves a change in net charge
        charge_difference = get_alchemical_charge_difference(mapping)
        if abs(charge_difference) > 1e-3:
            # Error out - this shouldn't be happening
            raise RuntimeError("charge change found!")
        
        complex_rfe_settings = get_complex_settings(forcefield)
        ligand_rfe_settings = get_ligand_settings(forcefield)
        for leg in ['solvent', 'complex']:
            # use the solvent and protein created above
            sysA_dict = {'ligand': mapping.componentA,
                         'solvent': solv}
            sysB_dict = {'ligand': mapping.componentB,
                         'solvent': solv}

            if leg == 'complex':
                sysA_dict['protein'] = prot
                sysB_dict['protein'] = prot

            sysA = openfe.ChemicalSystem(sysA_dict)
            sysB = openfe.ChemicalSystem(sysB_dict)

            name = (f"{leg}_{mapping.componentA.name}_"
                    f"{mapping.componentB.name}")

            if leg == 'complex':
                rbfe_protocol = RelativeHybridTopologyProtocol(settings=complex_rfe_settings)
            else:
                rbfe_protocol = RelativeHybridTopologyProtocol(settings=ligand_rfe_settings)

            transformation = openfe.Transformation(
                stateA=sysA,
                stateB=sysB,
                mapping=mapping,
                protocol=rbfe_protocol,
                name=name
            )
            transformations.append(transformation)

    # Create the alchemical network and write it to disk
    alchemical_network = openfe.AlchemicalNetwork(transformations)
    basedir = pathlib.Path('AlchemicalNetworks')
    basedir.mkdir(exist_ok=True)
    alchemical_network_json_fp = basedir / f"{forcefield}_elf10_alchemical_network.json"
    json.dump(
        alchemical_network.to_dict(),
        alchemical_network_json_fp.open(mode="w"),
        cls=tokenization.JSON_HANDLER.encoder
    )


@click.command
@click.option(
    '--pdb',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to the prepared PDB file of the protein",
)
def run_inputs(pdb):
    """
    Generate run json files for RBFE calculations

    Parameters
    ----------
    pdb : pathlib.Path
      A Path to a protein PDB file.
    """
    # Store the ligand network as a graphml file
    with open("ligand_network_elf10.graphml", mode='r') as f:
        ligand_network = gufe.LigandNetwork.from_graphml(f.read())

    # Create the solvent and protein components
    solv = openfe.SolventComponent()
    prot = openfe.ProteinComponent.from_pdb_file(str(pdb))

    for forcefield in ['openff-2.1.1']: # ['openff-2.2.1']:
        create_and_store_network(forcefield, ligand_network, solv, prot)


if __name__ == "__main__":
    run_inputs()
