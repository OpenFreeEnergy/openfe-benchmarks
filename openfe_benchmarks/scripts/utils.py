import logging
from rdkit import Chem

from openff.toolkit.utils.toolkit_registry import ToolkitRegistry
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper
import openfe

logger = logging.getLogger(__name__)

toolkit_registry = ToolkitRegistry([RDKitToolkitWrapper()])


def process_sdf(
    filename: str, return_dict: bool = False, key_type: str = "name"
) -> list[openfe.SmallMoleculeComponent] | dict[str, openfe.SmallMoleculeComponent]:
    """
    Process an SDF file and return a list or dictionary of SmallMoleculeComponent objects.

    Parameters
    ----------
    filename : str
        Path to the SDF file to process.
    return_dict : bool, optional
        If True, return a dictionary mapping keys to SmallMoleculeComponent objects.
        Defaults to False.
    key_type : str, optional
        The type of key to use when `return_dict` is True. One of:
        - ``"name"`` (default): use the molecule name as the key.
        - ``"smiles"``: use the canonical SMILES string as the key.
        - ``"inchikey"``: use the InChIKey as the key.

    Returns
    -------
    list[SmallMoleculeComponent] or dict[str, SmallMoleculeComponent]
        A list of SmallMoleculeComponent objects if `return_dict` is False.
        A dictionary mapping keys (determined by `key_type`) to
        SmallMoleculeComponent objects if `return_dict` is True.

    Raises
    ------
    ValueError
        If the SDF file cannot be loaded or if an unsupported `key_type` is
        provided.
    """
    supplier = Chem.SDMolSupplier(str(filename), removeHs=False)
    if supplier is None:
        raise ValueError(
            f"Failed to load molecules from the provided SDF file: {filename}"
        )
    else:
        molecules = [openfe.SmallMoleculeComponent(mol) for mol in supplier]
        if return_dict:
            if key_type == "name":
                molecules = {mol.name: mol for mol in molecules}
            elif key_type == "smiles":
                molecules = {
                    mol.to_openff().to_smiles(
                        explicit_hydrogens=True, toolkit_registry=toolkit_registry
                    ): mol
                    for mol in molecules
                }
            elif key_type == "inchikey":
                molecules = {
                    mol.to_openff().to_inchikey(
                        fixed_hydrogens=True, toolkit_registry=toolkit_registry
                    ): mol
                    for mol in molecules
                }
            else:
                raise ValueError(
                    f"Unsupported key_type '{key_type}'. "
                    "Must be one of: 'name', 'smiles', 'inchikey'."
                )

    return molecules
