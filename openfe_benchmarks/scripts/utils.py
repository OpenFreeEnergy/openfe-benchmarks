import warnings
import itertools

from rdkit import Chem
from loguru import logger

from openff.toolkit.utils.toolkit_registry import ToolkitRegistry

import openfe
from openfe.protocols.openmm_rfe.equil_rfe_methods import RelativeHybridTopologyProtocol
from openff.toolkit import (
    RDKitToolkitWrapper, AmberToolsToolkitWrapper, OpenEyeToolkitWrapper
)
from kartograf import KartografAtomMapper
warnings.filterwarnings(
    "ignore", 
    message="Partial charges have been provided, these will preferentially be used instead of generating new partial charges"
)

def compile_network_transformations(network, protocol, solvent, ligands_by_name, protein, cofactors):
    """
    Compile network transformations for a given alchemical network.

    Parameters
    ----------
    network : LigandNetwork
        The alchemical network containing edges to transform.
    protocol : RelativeHybridTopologyProtocol
        The protocol to use for the transformations.
    solvent : SolventComponent
        The solvent component for the transformations.
    ligands_by_name : dict[str, SmallMoleculeComponent]
        Dictionary mapping ligand names to SmallMoleculeComponent objects.
    protein : ProteinComponent
        The protein component for the transformations.
    cofactors : list[SmallMoleculeComponent] or None
        List of cofactor components, if any.

    Returns
    -------
    list[Transformation]
        List of transformations for the alchemical network.
    """
    
    transformations = []
    for edge in network.edges:
        new_edge = openfe.LigandAtomMapping(
            componentA=ligands_by_name[edge.componentA.name],
            componentB=ligands_by_name[edge.componentB.name], 
            componentA_to_componentB=edge.componentA_to_componentB,
            annotations=edge.annotations,
        )

#        # Get different settings depending on whether the transformation
#        # involves a change in net charge
#        charge_difference = get_alchemical_charge_difference(mapping)
#        if abs(charge_difference) > 1e-3:
#            # Error out - this shouldn't be happening
#            raise RuntimeError("charge change found!")
        complex_rfe_settings = protocol
        ligand_rfe_settings = protocol
    
        # create the transformations for the bound and solvent legs
        for leg in ["solvent", "complex"]:
            system_a_dict = {
                "ligand": new_edge.componentA, "solvent": solvent
            }
            system_b_dict = {
                "ligand": new_edge.componentB, "solvent": solvent
            }
            if leg == "complex":
                system_a_dict["protein"] = protein
                system_b_dict["protein"] = protein
    
                if cofactors is not None:
                    for i, cofactor in enumerate(cofactors):
                        cofactor_name = f"cofactor_{i}"
                        system_a_dict[cofactor_name] = cofactor
                        system_b_dict[cofactor_name] = cofactor
    
            system_a = openfe.ChemicalSystem(system_a_dict)
            system_b = openfe.ChemicalSystem(system_b_dict)
    
            name = f"{leg}_{new_edge.componentA.name}_{new_edge.componentB.name}"
            if leg == 'complex':
                rbfe_protocol = RelativeHybridTopologyProtocol(settings=complex_rfe_settings)
            else:
                rbfe_protocol = RelativeHybridTopologyProtocol(settings=ligand_rfe_settings)

            transformation = openfe.Transformation(
                stateA=system_a,
                stateB=system_b,
                mapping=new_edge,
                protocol=rbfe_protocol,
                name=name
            )
            transformations.append(transformation)

    return transformations


def process_sdf(filename: str, return_dict: bool = False) -> list[openfe.SmallMoleculeComponent] | dict[str, openfe.SmallMoleculeComponent]:
    """
    Process an SDF file and return a list or dictionary of SmallMoleculeComponent objects.

    Parameters
    ----------
    filename : str
        Path to the SDF file to process.
    return_dict : bool, optional
        If True, return a dictionary mapping ligand names to SmallMoleculeComponent objects.
        Defaults to False.

    Returns
    -------
    list[SmallMoleculeComponent] or dict[str, SmallMoleculeComponent]
        A list of SmallMoleculeComponent objects if `return_dict` is False.
        A dictionary mapping ligand names to SmallMoleculeComponent objects if `return_dict` is True.

    Raises
    ------
    ValueError
        If the SDF file cannot be loaded.
    """
    supplier = Chem.SDMolSupplier(str(filename), removeHs=False)
    if supplier is None:
        raise ValueError(f"Failed to load molecules from the provided SDF file: {filename}")
    else:
        molecules = [openfe.SmallMoleculeComponent(mol) for mol in supplier]
        if return_dict:
            molecules = {mol.name: mol for mol in molecules}
    return molecules


def generate_relative_network_from_names(ligands: dict[str, openfe.SmallMoleculeComponent],
                                         connections: list[tuple[str, str]],
                                         mappers: list[openfe.LigandAtomMapper] = [KartografAtomMapper()],
                                         scorer=None):
    """
    Generate a network from an input list of tuples each containing
    a pair of ligand names to connect with an edge.

    Parameters
    ----------
    ligands : dict[str, SmallMoleculeComponent]
        The ligands to create the network from.
    connections : list[tuple[str, str]]
        The list of edges to create with each node identified by its
        SmallMoleculeComponent name.
    mappers : list[LigandAtomMapper], optional
        Mappers to use. Defaults to [KartografAtomMapper()]. At least one is required.
    scorer : callable, optional
        A callable which returns a float for any LigandAtomMapping. Used to
        assign score to potential mappings, higher scores indicate worse
        mappings.

    Raises
    ------
    ValueError
        If no mapping can be found for a supplied edge.

    Returns
    -------
    LigandNetwork
        Network of SmallMoleculeComponent transformations.
    """
    edges = []

    for (lig1_name, lig2_name) in connections:
        print(lig1_name, lig2_name)
        for i, mapping in enumerate(itertools.chain.from_iterable(
            mapper.suggest_mappings(ligands[lig1_name], ligands[lig2_name])
            for mapper in mappers
        )):
            if not scorer:
                best_mapping = mapping
                break

            score = scorer(mapping)
            mapping = mapping.with_annotations({"score": score})

            if i == 0 or score < best_score:
                best_mapping = mapping
                best_score = score

        if best_mapping is None:
            raise ValueError(f"No mapping for pair {entry}")

        edges.append(best_mapping)

    return openfe.LigandNetwork(edges)


def update_partial_charges(smc_list, partial_charge_method, off_wrapper):
    """
    Updates the partial charges of a list of molecules using the specified method and toolkit wrapper.
    Parameters
    ----------
    smc_list : list
        A list of molecules for which partial charges need to be updated.
    partial_charge_method : str
        The name of the partial charge calculation method to use.
    off_wrapper : object
        An Open Force Field (OFF) toolkit wrapper to handle the charge assignment.
    Returns
    -------
    list
        A list of molecules with updated partial charges.
    """
    
    return [assign_partial_charges(mol, partial_charge_method, off_wrapper) for mol in smc_list]


def assign_partial_charges(ofe_smc, partial_charge_method, off_wrapper):
    """
    Generate partial charges for a SmallMoleculeComponent using
    the input conformer and desired method and backend.

    Parameters
    ----------
    ofe_smc : SmallMoleculeComponent
        The small molecule component to assign partial charges to.
    partial_charge_method : str
        The method to use for partial charge assignment.
    off_wrapper : ToolkitRegistry
        The toolkit registry to use for charge assignment.

    Returns
    -------
    SmallMoleculeComponent
        The small molecule component with assigned partial charges.
    """
    logger.debug(f"INFO: generating partial charges for ligand {ofe_smc.name} -- this may be slow")
    offmol = ofe_smc.to_openff()
    offmol.assign_partial_charges(
        partial_charge_method=partial_charge_method,
        toolkit_registry=off_wrapper,
    )
    return openfe.SmallMoleculeComponent.from_openff(offmol)


def get_alchemical_charge_difference(mapping) -> int:
    """
    Checks and returns the difference in formal charge between state A and B.

    Parameters
    ----------
    mapping : dict[str, ComponentMapping]
        Dictionary of mappings between transforming components, i.e., edge.

    Returns
    -------
    int
        The formal charge difference between states A and B.
        This is defined as sum(charge state A) - sum(charge state B).
    """
    chg_A = Chem.rdmolops.GetFormalCharge(
        mapping.componentA.to_rdkit()
    )
    chg_B = Chem.rdmolops.GetFormalCharge(
        mapping.componentB.to_rdkit()
    )

    return chg_A - chg_B