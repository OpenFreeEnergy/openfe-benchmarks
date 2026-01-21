import click
import pathlib
import logging
import warnings
import openfe
from rdkit import Chem
from openff.toolkit import OpenEyeToolkitWrapper
import gufe


logger = logging.getLogger(__name__)
warnings.filterwarnings(
    "ignore", 
    message="Partial charges have been provided, these will preferentially be used instead of generating new partial charges"
)


def gen_charges(smc):
    """
    Generate AM1BCC partial charges for a SmallMoleculeComponent using
    the input conformer and antechamber as backend.
    """
    print(f"INFO: generating partial charges for ligand {smc.name} -- this may be slow")
    offmol = smc.to_openff()

    offmol.assign_partial_charges(
        partial_charge_method="am1bccelf10",
        toolkit_registry=OpenEyeToolkitWrapper(),
    )
    return openfe.SmallMoleculeComponent.from_openff(offmol)


def gen_ligand_network(smcs):
    """
    Creates the Lomap ligand network using the KartografAtomMapper and
    the Lomap scorer.

    Parameters
    ----------
    smcs : list[SmallMoleculeComponents]
      List of SmallMoleculeComponents of the ligands.

    Returns
    -------
    openfe.LigandNetwork
      The Lomap generated LigandNetwork.
    """
    with open('ligand_network_am1bcc.graphml', 'r') as f:
        old_network = gufe.LigandNetwork.from_graphml(f.read())

    ligands_by_name = {}
    for mol in smcs:
        ligands_by_name[mol.name] = mol

    new_edges = []

    for edge in old_network.edges:
        new_edge = gufe.LigandAtomMapping(
            componentA=ligands_by_name[edge.componentA.name],
            componentB=ligands_by_name[edge.componentB.name],
            componentA_to_componentB=edge.componentA_to_componentB,
        )
        new_edges.append(new_edge)


    ligand_network = gufe.LigandNetwork(new_edges)

    return ligand_network


@click.command
@click.option(
    '--ligands',
    type=click.Path(dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to the prepared SDF file containing the ligands",
)
def run_inputs(ligands):
    """
    Generate run json files for RBFE calculations

    Parameters
    ----------
    ligands : pathlib.Path
      A Path to a ligands SDF.
    """
    # Create the small molecule components of the ligands
    rdmols = [mol for mol in Chem.SDMolSupplier(str(ligands), removeHs=False)]
    smcs = [openfe.SmallMoleculeComponent.from_rdkit(mol) for mol in rdmols]
    # Generate the partial charges
    logger.info("Generating partial charges for ligands")
    smcs = [gen_charges(smc) for smc in smcs]

    # Create ligand network
    ligand_network = gen_ligand_network(smcs)

    # Store the ligand network as a graphml file
    with open("ligand_network_elf10.graphml", mode='w') as f:
        f.write(ligand_network.to_graphml())


if __name__ == "__main__":
    run_inputs()
