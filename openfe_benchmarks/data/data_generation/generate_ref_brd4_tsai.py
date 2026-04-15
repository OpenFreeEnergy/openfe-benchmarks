"""
For the BRD4 system from Tsai et al. extract the reference data for the ligands in
the system and save it to the output directory.

The reference data will be saved to JSON using gufe which allows us to include units. The structure is a dictionary
where each key is the ligand name and the value is a dict of the experimental DG and the associated uncertainty.
"""
import click
import pathlib
from gufe.tokenization import JSON_HANDLER
import pandas as pd
from openff.units import unit
import json
from openff.toolkit import Molecule


@click.command()
@click.option("--csv-name", type=str, required=True, help="The csv file name with exp. DG values.")
@click.option("--out-dir", type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=pathlib.Path), required=True, help="The output dir name")
@click.option("--input-sdf", type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=pathlib.Path), required=False, help="The input sdf file containing the ligands from this network, used to extract identifiers.")
def main(csv_name: str, out_dir: pathlib.Path, input_sdf: pathlib.Path):
    """
    Extract the reference data for the given system.

    Parameters
    ----------
    csv_name: str
        The csv file name with exp. DG values.
    out_dir : pathlib.Path
        The output dir name where the extracted reference data will be saved.
    input_sdf
        The input sdf file containing the ligands from this network, used to extract identifiers which will be stored in the reference data.
    """
    tsai_et_al_doi = "https://doi.org/10.1021/acs.jcim.5c02204"
    # load the ref data from a csv file
    ref_dg_data = pd.read_csv(csv_name)
    # load the ligands with openff toolkit to extract the identifiers
    molecules = Molecule.from_file(input_sdf.as_posix(), allow_undefined_stereo=True)
    molecule_by_name = {molecule.name.lower(): molecule for molecule in molecules}
    print(molecule_by_name)
    ref_dataset = {}
    for _, row in ref_dg_data.iterrows():
        ligand_name = row["Node"]
        if ligand_name not in molecule_by_name:
            # this ligand was not included in the benchmark dataset so skip it
            continue
        ligand_data = {
            "dg": row["Expt"] * unit.kilocalories_per_mole,
            "canonical_smiles": molecule_by_name[ligand_name].to_smiles(isomeric=True),
            # do we need to be tautomer specific?
            "inchikey": molecule_by_name[ligand_name].to_inchikey(fixed_hydrogens=True),
            "reference": tsai_et_al_doi,
        }
        ref_dataset[ligand_name] = ligand_data

    # save the extracted reference data to the output directory as json using gufe serialization
    out_file = out_dir / "experimental_binding_data.json"
    with open(out_file, "w") as f:
        json.dump(ref_dataset, f, indent=4, cls=JSON_HANDLER.encoder)
    print(f"Saved extracted reference data to {out_file}")

if __name__ == "__main__":
    main()