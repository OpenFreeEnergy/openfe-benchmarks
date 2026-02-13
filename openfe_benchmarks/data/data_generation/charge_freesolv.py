"""
Generate partial charges for a set of molecules using OpenFE's bulk charge assignment utility will also add software version metadata to each ligand
as an sdf property.
"""
import pathlib
import os

from openfe.protocols.openmm_utils.charge_generation import assign_offmol_partial_charges
import click
from rdkit import Chem
from gufe import SmallMoleculeComponent
import openfe
from openff import toolkit
from openff.utilities.provenance import get_ambertools_version
import json
import tqdm


@click.command()
@click.option("--input-dir", type=click.Path(exists=True, dir_okay=True, path_type=pathlib.Path), required=True,
              help="Path to the input SDF files containing the molecules to be charged.")
@click.option("--output-dir", type=click.Path(dir_okay=True, file_okay=False, exists=True, path_type=pathlib.Path),
              required=True, help="Path to the output folder the SDF files with charged molecules will be saved.")
@click.option("--charge-method", type=click.Choice(["am1bcc_at", "am1bccelf10_oe", "nagl_off", "am1bcc_oe"]),
              default="am1bcc_at", help="The method to use for charge assignment.")
@click.option("--nagl-model", type=str, default=None,
              help="Path to the NAGL model to use for charge assignment when using the 'nagl_off' method if None the latest model will be used.")
def main(input_dir: pathlib.Path, output_dir: pathlib.Path, charge_method: str, nagl_model: None | str):
    """Generate partial charges for a set of molecules using OpenFE's charge assignment utility.

    Parameters
    ----------
    input_dir : pathlib.Path
        Path to the input SDF files containing the molecules.
    output_dir : pathlib.Path
        Directory where the output SDF file with charged molecules will be saved.
    charge_method : str
        The method to use for charge assignment. Options include:

        - 'am1bcc_at': AM1BCC applied with AmberTools on the input conformer
        - 'am1bcc_oe': AM1BCC applied with OpenEye Toolkit on the input conformer
        - 'am1bccelf10_oe': AM1BCC Elf10 applied with OpenEye Toolkit using 500 conformers
        - 'nagl_off': NAGL charges applied with OpenFF-Toolkit

    n_cores : int
        Number of CPU cores to use for parallel processing.
    nagl_model : str
        Model *.pt file (i.e., "openff-gnn-am1bcc-1.0.0.pt"), optionally with path, for the NAGL model to use for
        charge assignment when using the ``'nagl'`` method. If None the latest model will be used. See 
        [OpenFF NAGL](https://docs.openforcefield.org/projects/nagl-models) documentation for more detail.

    Notes
    -----
    - Molecules are loaded using the OpenFF-Toolkit to avoid issues with sterochemistry perception in other toolkits, the input SDF files should have a single molecule with a single conformer.
    - Antechamber will be used for the am1bcc_at charge assignment method, the charges are calculated at the input geometry.
    - OpenEye toolkit is required for am1bccelf10_oe charge assignment method and am1bcc_oe.
    - The output SDF file will include software version metadata as a property for each ligand and will be named <input_name>_<charge_method>.sdf

    """
    mols = []
    for input_path in input_dir.glob("*.sdf"):
        # should be a single molecule per file
        off_mol = toolkit.Molecule.from_file(input_path.as_posix(), allow_undefined_stereo=True)
        mols.append(SmallMoleculeComponent.from_openff(off_mol))

    # construct the toolkit backend
    method_to_backend = {
        "am1bcc_at": "ambertools",
        "am1bcc_oe": "openeye",
        "am1bccelf10_oe": "openeye",
        "nagl_off": "rdkit"
    }
    backend = method_to_backend[charge_method]

    # convert the charge method to the expected format for openff
    openff_charge_method = charge_method.split("_")[0]

    # we need to generate conformers for am1bccelf10_oe or use the input conformer for other methods which is the None case
    generate_n_conformers = None if charge_method != "am1bccelf10_oe" else 500
    charged_ligands = []
    failed_molecules = []
    for mol in tqdm.tqdm(mols, desc="Generating charges", ncols=80):
        try:
            charged_molecule = assign_offmol_partial_charges(
                offmol=mol.to_openff(),
                overwrite=True,
                method=openff_charge_method,
                toolkit_backend=backend,
                generate_n_conformers=generate_n_conformers,
                nagl_model=nagl_model
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
    prov_nagl_model = None
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
            prov_nagl_model = get_models_by_type(model_type="am1bcc", production_only=True)[-1].name
        else:
            prov_nagl_model = os.path.split(nagl_model)
        provenance["nagl_version"] = str(nagl.__version__)
        provenance["nagl_model"] = prov_nagl_model

    # construct the output path
    method_to_name = {
        "am1bcc_at": "antechamber_am1bcc",
        "am1bccelf10_oe": "openeye_am1bccelf10",
        "nagl_off": f"nagl_{prov_nagl_model}",
        "am1bcc_oe": "openeye_am1bcc"
    }

    output_path = output_dir / f"ligands_{method_to_name[charge_method]}.sdf"
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
