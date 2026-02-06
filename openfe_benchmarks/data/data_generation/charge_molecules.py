"""
Generate partial charges for a set of molecules using OpenFE's bulk charge assignment utility will also add software version metadata to each ligand
as an sdf property.
"""

import pathlib
import os

from openfe.protocols.openmm_utils.charge_generation import bulk_assign_partial_charges
import click
from rdkit import Chem
from gufe import SmallMoleculeComponent
import openfe
from openff import toolkit
from openff.utilities.provenance import get_ambertools_version
import json


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
@click.option(
    "--n-cores",
    type=int,
    default=1,
    help="Number of CPU cores to use for parallel processing.",
)
def main(
    input_path: pathlib.Path,
    output_dir: pathlib.Path,
    charge_method: str,
    nagl_model: None | str,
    n_cores: int,
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

    n_cores : int
        Number of CPU cores to use for parallel processing.
    nagl_model : str
        Model *.pt file (i.e., "openff-gnn-am1bcc-1.0.0.pt"), optionally with path, for the NAGL model to use for
        charge assignment when using the ``'nagl'`` method. If None the latest model will be used. See
        [OpenFF NAGL](https://docs.openforcefield.org/projects/nagl-models) documentation for more detail.

    Notes
    -----
    - Antechamber will be used for the am1bcc_at charge assignment method, the charges are calculated at the input geometry.
    - OpenEye toolkit is required for am1bccelf10_oe charge assignment method and am1bcc_oe.
    - The output SDF file will include software version metadata as a property for each ligand and will be named <input_name>_<charge_method>.sdf

    """
    with Chem.SDMolSupplier(input_path.as_posix(), removeHs=False) as supplier:
        mols = [
            SmallMoleculeComponent.from_rdkit(mol)
            for mol in supplier
            if mol is not None
        ]
        input_order = [mol.name for mol in mols]
        # construct the toolkit backend
        method_to_backend = {
            "am1bcc_at": "ambertools",
            "am1bcc_oe": "openeye",
            "am1bccelf10_oe": "openeye",
            "nagl_off": "rdkit",
        }
        backend = method_to_backend[charge_method]

        # convert the charge method to the expected format for openff
        openff_charge_method = charge_method.split("_")[0]

        # we need to generate conformers for am1bccelf10_oe or use the input conformer for other methods which is the None case
        generate_n_conformers = None if charge_method != "am1bccelf10_oe" else 500

        charged_ligands = bulk_assign_partial_charges(
            molecules=mols,
            overwrite=True,
            method=openff_charge_method,
            toolkit_backend=backend,
            processors=n_cores,
            generate_n_conformers=generate_n_conformers,
            nagl_model=nagl_model,
        )

        # the multiprocessing can shuffle the order so we need to restore the input order
        charged_ligands = sorted(
            charged_ligands, key=lambda x: input_order.index(x.name)
        )

        # for each ligand stamp the provenance info as sdf property and write to output sdf
        # generate the provenance info
        provenance = {
            "openfe_version": openfe.__version__,
            "openff_toolkit_version": toolkit.__version__,
            "rdkit_version": Chem.rdBase.rdkitVersion,
            "charge_method": charge_method,
        }
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

        # construct the output path
        method_to_name = {
            "am1bcc_at": "antechamber_am1bcc",
            "am1bccelf10_oe": "openeye_am1bccelf10",
            "nagl_off": f"nagl_{nagl_model}",
            "am1bcc_oe": "openeye_am1bcc",
        }

        output_path = (
            output_dir / f"{input_path.stem}_{method_to_name[charge_method]}.sdf"
        )
        with Chem.SDWriter(str(output_path)) as writer:
            for ligand in charged_ligands:
                rdkit_mol = ligand.to_rdkit()
                # add software version metadata as sdf property
                rdkit_mol.SetProp("charge_provenance", json.dumps(provenance))
                writer.write(rdkit_mol)


if __name__ == "__main__":
    main()
