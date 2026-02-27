"""
Generate partial charges for all molecules in a ``ligands.sdf`` file using
OpenFE's bulk charge assignment utility, then write a charged SDF file with
software version metadata stamped as an SDF property on each molecule.

Supported charge methods
------------------------
- ``am1bcc_at``      — AM1BCC via AmberTools antechamber (input conformer)
- ``am1bcc_oe``      — AM1BCC via OpenEye Toolkit (input conformer)
- ``am1bccelf10_oe`` — AM1BCC-ELF10 via OpenEye Toolkit (500 conformers)
- ``nagl_off``       — NAGL graph-neural-network charges via OpenFF-Toolkit
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


@click.command()
@click.option(
    "--dir-path",
    type=click.Path(exists=True, dir_okay=True, path_type=pathlib.Path),
    required=True,
    help="Path to ligands.sdf file containing the molecules to be charged.",
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
def main(dir_path: pathlib.Path, charge_method: str, nagl_model: None | str):
    """Generate partial charges for all molecules in ``<dir_path>/ligands.sdf``.

    Reads the input SDF, assigns partial charges with the chosen method and
    toolkit backend, then writes ``<dir_path>/ligands_<method_label>.sdf``
    where each record carries a ``charge_provenance`` SDF property containing
    a JSON object with software versions and the charge method used.

    Failed molecules are skipped and reported at the end.

    Parameters
    ----------
    dir_path : pathlib.Path
        Directory containing ``ligands.sdf``. The output SDF is also written
        here as ``ligands_<method_label>.sdf``.
    charge_method : str
        Charge assignment method. One of:

        - ``'am1bcc_at'``: AM1BCC via AmberTools antechamber, calculated
          at the input conformer geometry.
        - ``'am1bcc_oe'``: AM1BCC via OpenEye Toolkit, calculated at the
          input conformer geometry (requires OpenEye licence).
        - ``'am1bccelf10_oe'``: AM1BCC-ELF10 via OpenEye Toolkit using 500
          generated conformers (requires OpenEye licence).
        - ``'nagl_off'``: NAGL graph-neural-network AM1BCC charges via
          OpenFF-Toolkit + RDKit.

    nagl_model : str or None
        Path to or filename of a NAGL model ``*.pt`` file, e.g.
        ``"openff-gnn-am1bcc-1.0.0.pt"``. Only used when
        ``charge_method='nagl_off'``. If ``None``, the latest production model
        is selected automatically. See the
        `OpenFF NAGL Models <https://docs.openforcefield.org/projects/nagl-models>`_
        documentation for available models.

    Notes
    -----
    - Molecules are loaded via the OpenFF Toolkit to avoid stereochemistry
      perception issues in other toolkits. The input SDF should contain one
      conformer per molecule.
    - For ``am1bccelf10_oe``, 500 conformers are generated internally before
      charge assignment; all other methods use the input conformer.
    - The ``charge_provenance`` SDF property records ``openfe_version``,
      ``openff_toolkit_version``, ``rdkit_version``, ``charge_method``, and
      any toolkit-specific version fields (e.g. ``ambertools_version``,
      OpenEye module versions, or ``nagl_version`` + ``nagl_model``).
    """
    input_sdf = dir_path / "ligands.sdf"
    off_mols = toolkit.Molecule.from_file(
        input_sdf.as_posix(), allow_undefined_stereo=True
    )
    # Handle both single molecule and multiple molecules
    if not isinstance(off_mols, list):
        off_mols = [off_mols]
    mols = [SmallMoleculeComponent.from_openff(mol) for mol in off_mols]

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
            prov_nagl_model = get_models_by_type(
                model_type="am1bcc", production_only=True
            )[-1].name
        else:
            prov_nagl_model = os.path.split(nagl_model)
        provenance["nagl_version"] = str(nagl.__version__)
        provenance["nagl_model"] = prov_nagl_model

    # construct the output path
    method_to_name = {
        "am1bcc_at": "antechamber_am1bcc",
        "am1bccelf10_oe": "openeye_am1bccelf10",
        "nagl_off": f"nagl_{prov_nagl_model}",
        "am1bcc_oe": "openeye_am1bcc",
    }

    output_path = dir_path / f"ligands_{method_to_name[charge_method]}.sdf"
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
