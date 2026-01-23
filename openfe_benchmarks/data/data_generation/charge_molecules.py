"""
Generate partial charges for a set of molecules using OpenFE's bulk charge assignment utility will also add software version metadata to each ligand
as an sdf property.
"""
import pathlib

from openfe.protocols.openmm_utils.charge_generation import bulk_assign_partial_charges
import click
from rdkit import Chem
from gufe import SmallMoleculeComponent
import openfe
import openff.toolkit
from openff.utilities.provenance import get_ambertools_version
import json

@click.command()
@click.option("--input-path", type=click.Path(exists=True, dir_okay=False, path_type=pathlib.Path), required=True, help="Path to the input SDF file containing the molecules to be charged.")
@click.option("--output-dir", type=click.Path(dir_okay=True, file_okay=False, exists=True, path_type=pathlib.Path), required=True, help="Path to the output folder the SDF file with charged molecules will be saved.")
@click.option("--charge-method", type=click.Choice(["am1bcc", "am1bccelf10"]), default="am1bcc", help="The method to use for charge assignment.")
@click.option("--n-cores", type=int, default=1, help="Number of CPU cores to use for parallel processing.")
def main(input_path: pathlib.Path, output_dir: pathlib.Path, charge_method: str, n_cores: int):
    """Generate partial charges for a set of molecules using OpenFE's bulk charge assignment utility.

    Parameters
    ----------
    input_path : pathlib.Path
        Path to the input SDF file containing the molecules.
    output_dir : pathlib.Path
        Directory where the output SDF file with charged molecules will be saved.
    charge_method : str
        The method to use for charge assignment. Options are 'am1bcc', 'am1bccelf10', 'gasteiger', and 'am1-mulliken'.
    n_cores : int
        Number of CPU cores to use for parallel processing.

    Notes
    -----
    - Antechamber will be used for the am1bcc charge assignment method, the charges are calculated at the input geometry.
    - openeye toolkits is required for am1bccelf10 charge assignment method.
    - The output SDF file will include software version metadata as a property for each ligand and will be named <input_name>_<charge_method>.sdf

    """
    with Chem.SDMolSupplier(input_path.as_posix(), removeHs=False) as supplier:
        mols = [SmallMoleculeComponent.from_rdkit(mol) for mol in supplier if mol is not None]

        # construct the toolkit backend
        method_to_backend = {
            "am1bcc": "ambertools",
            "am1bccelf10": "openeye",
        }
        # we need to generate conformers for am1bccelf10 or use the input conformer for other methods which is the None case
        generate_n_conformers = None if charge_method != "am1bccelf10" else 500
        charged_ligands = bulk_assign_partial_charges(
            molecules=mols,
            overwrite=True,
            method=charge_method,
            toolkit_backend=method_to_backend[charge_method],
            processors=n_cores,
            generate_n_conformers=generate_n_conformers,
            nagl_model=None
        )
        # for each ligand stamp the provenance info as sdf property and write to output sdf
        # generate the provenance info
        provenance = {
            "openfe_version": openfe.__version__,
            "openff_toolkit_version": openff.toolkit.__version__,
            "rdkit_version": Chem.rdBase.rdkitVersion,
            "charge_method": charge_method,
        }
        if method_to_backend[charge_method] == "ambertools":
            provenance["ambertools_version"] = get_ambertools_version()

        elif method_to_backend[charge_method] == "openeye":
            from openeye import oeomega, oequacpac
            provenance["oeomega"] = str(oeomega.OEOmegaGetVersion())
            provenance["oequacpac"] = str(oequacpac.OE_OEQUACPAC_VERSION)


        method_to_name = {
            "am1bcc": "antechamber_am1bcc",
            "am1bccelf10": "openeye_am1bccelf10",
        }
        output_path = output_dir / f"{input_path.stem}_{method_to_name[charge_method]}.sdf"
        with Chem.SDWriter(str(output_path)) as writer:
            for ligand in charged_ligands:
                rdkit_mol = ligand.to_rdkit()
                # add software version metadata as sdf property
                rdkit_mol.SetProp("charge_provenance", json.dumps(provenance))
                writer.write(rdkit_mol)

if __name__ == "__main__":
    main()