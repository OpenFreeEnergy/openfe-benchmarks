import warnings

from rdkit import Chem
from loguru import logger

import openfe
from openfe.protocols.openmm_rfe.equil_rfe_methods import RelativeHybridTopologyProtocol
warnings.filterwarnings(
    "ignore", 
    message="Partial charges have been provided, these will preferentially be used instead of generating new partial charges"
)

#def compile_network_transformations(network, protocol, solvent, ligands_by_name, protein, cofactors):
#    """
#    Compile network transformations for a given alchemical network.
#
#    Parameters
#    ----------
#    network : LigandNetwork
#        The alchemical network containing edges to transform.
#    protocol : RelativeHybridTopologyProtocol
#        The protocol to use for the transformations.
#    solvent : SolventComponent
#        The solvent component for the transformations.
#    ligands_by_name : dict[str, SmallMoleculeComponent]
#        Dictionary mapping ligand names to SmallMoleculeComponent objects.
#    protein : ProteinComponent
#        The protein component for the transformations.
#    cofactors : list[SmallMoleculeComponent] or None
#        List of cofactor components, if any.
#
#    Returns
#    -------
#    list[Transformation]
#        List of transformations for the alchemical network.
#    """
#    
#    transformations = []
#    for edge in network.edges:
#        new_edge = openfe.LigandAtomMapping(
#            componentA=ligands_by_name[edge.componentA.name],
#            componentB=ligands_by_name[edge.componentB.name], 
#            componentA_to_componentB=edge.componentA_to_componentB,
#            annotations=edge.annotations,
#        )
#
##        # Get different settings depending on whether the transformation
## https://github.com/OpenFreeEnergy/gufe/blob/f10eb66916e96af5ed9595366f9ca91fba2b9944/gufe/mapping/ligandatommapping.py#L310
##        # involves a change in net charge
##        charge_difference = get_alchemical_charge_difference(mapping)
##        if abs(charge_difference) > 1e-3:
##            # Error out - this shouldn't be happening
##            raise RuntimeError("charge change found!")
#        complex_rfe_settings = protocol
#        ligand_rfe_settings = protocol
#    
#        # create the transformations for the bound and solvent legs
#        for leg in ["solvent", "complex"]:
#            system_a_dict = {
#                "ligand": new_edge.componentA, "solvent": solvent
#            }
#            system_b_dict = {
#                "ligand": new_edge.componentB, "solvent": solvent
#            }
#            if leg == "complex":
#                system_a_dict["protein"] = protein
#                system_b_dict["protein"] = protein
#    
#                if cofactors is not None:
#                    for i, cofactor in enumerate(cofactors):
#                        cofactor_name = f"cofactor_{i}"
#                        system_a_dict[cofactor_name] = cofactor
#                        system_b_dict[cofactor_name] = cofactor
#    
#            system_a = openfe.ChemicalSystem(system_a_dict)
#            system_b = openfe.ChemicalSystem(system_b_dict)
#    
#            name = f"{leg}_{new_edge.componentA.name}_{new_edge.componentB.name}"
#            if leg == 'complex':
#                rbfe_protocol = RelativeHybridTopologyProtocol(settings=complex_rfe_settings)
#            else:
#                rbfe_protocol = RelativeHybridTopologyProtocol(settings=ligand_rfe_settings)
#
#            transformation = openfe.Transformation(
#                stateA=system_a,
#                stateB=system_b,
#                mapping=new_edge,
#                protocol=rbfe_protocol,
#                name=name
#            )
#            transformations.append(transformation)
#
#    return transformations


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