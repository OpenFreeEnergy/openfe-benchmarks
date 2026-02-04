import warnings

from rdkit import Chem
from loguru import logger

import openfe
from openfe.protocols.openmm_rfe.equil_rfe_methods import RelativeHybridTopologyProtocol
warnings.filterwarnings(
    "ignore", 
    message="Partial charges have been provided, these will preferentially be used instead of generating new partial charges"
)


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