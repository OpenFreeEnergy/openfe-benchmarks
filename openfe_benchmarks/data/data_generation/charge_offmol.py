"""
Generate partial charges for a set of molecules using OpenFE's bulk charge assignment utility will also add software version metadata to each ligand
as an sdf property.
"""

import pathlib
import os

from openfe.protocols.openmm_utils.charge_generation import (
    assign_offmol_partial_charges,
)
import click
from rdkit import Chem
from gufe import SmallMoleculeComponent
import openfe
from openff import toolkit
from openff.utilities.provenance import get_ambertools_version
import json
import tqdm


CHARGE_METHODS = {
    "am1bcc_at": {
        "method": "am1bcc",
        "backend": "ambertools",
        "generate_n_conformers": None,
        "output_name": "antechamber_am1bcc",
    },
    "am1bcc_oe": {
        "method": "am1bcc",
        "backend": "openeye",
        "generate_n_conformers": None,
        "output_name": "openeye_am1bcc",
    },
    "am1bccelf10_oe": {
        "method": "am1bccelf10",
        "backend": "openeye",
        "generate_n_conformers": 500,
        "output_name": "openeye_am1bccelf10",
    },
    "nagl_off": {
        "method": "nagl",
        "backend": "rdkit",
        "generate_n_conformers": None,
        "output_name": None,
    },
}

@click.command()
@click.option(
    "--input-path",
    type=click.Path(exists=True, dir_okay=False, path_type=pathlib.Path),
    required=True,
    help="Path to the input SDF file containing the molecules to be charged.",
)
@click.option(
    "--output-dir",
    type=click.Path(
        dir_okay=True, file_okay=False, exists=True, path_type=pathlib.Path
    ),
    required=True,
    help="Path to the output folder the SDF file with charged molecules will be saved.",
)
@click.option(
    "--charge-method",
    type=click.Choice(["am1bcc_at", "am1bccelf10_oe", "nagl_off", "am1bcc_oe"]),
    default="am1bcc_at",
    help="The method to use for charge assignment.",
)
@click.option(
    "--nagl-model",
    type=str,
    default=None,
    help="Path to the NAGL model to use for charge assignment when using the 'nagl_off' method if None the latest model will be used.",
)
def main(
    input_path: pathlib.Path,
    output_dir: pathlib.Path,
    charge_method: str,
    nagl_model: None | str,
):
    """Generate partial charges for a set of molecules using OpenFE's bulk charge assignment utility.

    Parameters
    ----------
    input_path : pathlib.Path
        Path to the input SDF file containing the molecules.
    output_dir : pathlib.Path
        Directory where the output SDF file with charged molecules will be saved.
    charge_method : str
        The method to use for charge assignment. Options include:

        - 'am1bcc_at': AM1BCC applied with AmberTools on the input conformer
        - 'am1bcc_oe': AM1BCC applied with OpenEye Toolkit on the input conformer
        - 'am1bccelf10_oe': AM1BCC Elf10 applied with OpenEye Toolkit using 500 conformers
        - 'nagl_off': NAGL charges applied with OpenFF-Toolkit

    nagl_model : str
        Model *.pt file (i.e., "openff-gnn-am1bcc-1.0.0.pt"), optionally with path, for the NAGL model to use for
        charge assignment when using the ``'nagl'`` method. If None the latest model will be used. See
        [OpenFF NAGL](https://docs.openforcefield.org/projects/nagl-models) documentation for more detail.

    Notes
    -----
    - Antechamber will be used for the am1bcc_at charge assignment method, the charges are calculated at the input geometry.
    - OpenEye toolkit is required for am1bccelf10_oe charge assignment method and am1bcc_oe.
    - The output SDF file will include software version metadata as a property for each ligand and will be named <input_name>_<charge_method>.sdf
    - Water is always skipped

    """
    off_mols = toolkit.Molecule.from_file(
        input_path.as_posix(), allow_undefined_stereo=True
    )
    # Handle both single molecule and multiple molecules
    if not isinstance(off_mols, list):
        off_mols = [off_mols]
    mols = [SmallMoleculeComponent.from_openff(mol) for mol in off_mols]

    charge_settings = CHARGE_METHODS[charge_method]

    charged_ligands = []
    failed_molecules = []
    for mol in tqdm.tqdm(mols, desc="Generating charges", ncols=80):
        # water partial charges should come from a specific water model
        if mol.to_openff().to_inchikey(fixed_hydrogens=True) == "XLYOFNOQVPJJNP-UHFFFAOYNA-N": # water
            continue

        try:
            charged_molecule = assign_offmol_partial_charges(
                offmol=mol.to_openff(),
                overwrite=True,
                method=charge_settings["method"],
                toolkit_backend=charge_settings["backend"],
                generate_n_conformers=charge_settings["generate_n_conformers"],
                nagl_model=nagl_model,
            )
        except Exception as e:
            print(f"Failed to generate charges for molecule {mol.name} with error: {e}")
            failed_molecules.append(mol)
            continue

        charged_ligands.append(SmallMoleculeComponent.from_openff(charged_molecule))

    # for each ligand stamp the provenance info as sdf property and write to output sdf
    # generate the provenance info
    provenance = {
        "openfe_version": openfe.__version__,
        "openff_toolkit_version": toolkit.__version__,
        "rdkit_version": Chem.rdBase.rdkitVersion,
        "charge_method": charge_method,
    }

    backend = charge_settings["backend"]
    output_name = charge_settings["output_name"]

    if backend == "ambertools":
        provenance["ambertools_version"] = get_ambertools_version()

    elif backend == "openeye":
        from openeye import oeomega, oequacpac

        provenance["oeomega"] = str(oeomega.OEOmegaGetVersion())
        provenance["oequacpac"] = str(oequacpac.OE_OEQUACPAC_VERSION)
    elif backend == "rdkit" and charge_method == "nagl_off":
        from openff import nagl
        from openff.nagl_models import get_models_by_type

        if nagl_model is None:
            # get the latest production nagl model
            nagl_model = get_models_by_type(
                model_type="am1bcc", production_only=True
            )[-1].name
        else:
            nagl_model = os.path.split(nagl_model)
        provenance["nagl_version"] = str(nagl.__version__)
        provenance["nagl_model"] = nagl_model
        output_name = f"nagl_{nagl_model_name}"

    # construct the output path
    output_path = (
        output_dir / f"{input_path.stem}_{output_name}.sdf"
    )
    with Chem.SDWriter(str(output_path)) as writer:
        for ligand in charged_ligands:
            rdkit_mol = ligand.to_rdkit()
            # add software version metadata as sdf property
            rdkit_mol.SetProp("charge_provenance", json.dumps(provenance))
            writer.write(rdkit_mol)

    for failed_mol in failed_molecules:
        print(f"Failed molecules: {failed_mol.name}")

if __name__ == "__main__":
    main()
