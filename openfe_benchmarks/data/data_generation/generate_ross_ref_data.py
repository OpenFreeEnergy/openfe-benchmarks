"""
For the given Schrodinger benchmark system extract the reference data for the ligands
 and save it to the output directory.

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
@click.option("--system-group", type=str, required=True, help="The system group of system names used to determine the edges, i.e., JACS/MERCK.")
@click.option("--system-name", type=str, required=True, help="The system name used to determine the edges, i.e., TYK2.")
@click.option("--csv-name", type=str, required=True, help="The name of the exp ligand csv file in the Schrodinger repo.")
@click.option("--out-dir", type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=pathlib.Path), required=True, help="The output dir name")
@click.option("--input-sdf", type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=pathlib.Path), required=False, help="The input sdf file containing the ligands from this network, used to extract identifiers.")
def main(system_group: str, system_name: str, csv_name: str, out_dir: pathlib.Path, input_sdf: pathlib.Path):
    """
    Extract the reference data for the given industry benchmark system.

    Parameters
    ----------
    system_group : str
        The industry system group ie JACS/MERCK used to find the reference data.
    system_name : str
        The industry system name ie TYK2 used to find the reference data.
    csv_name: str
        The name of the exp ligand csv file in the Schrodinger repo.
    out_dir : pathlib.Path
        The output dir name where the extracted reference data will be saved.
    input_sdf
        The input sdf file containing the ligands from this network, used to extract identifiers which will be stored in the reference data.
    """
    # tag each entry with the source, currently we link to the aggregated reference data but this should be changed if we refine it in future.
    ross_et_al_doi = "https://doi.org/10.1038/s42004-023-01019-9"
    # load the ref data from the industry benchmark results stored on github
    csv_file = f"https://raw.githubusercontent.com/schrodinger/public_binding_free_energy_benchmark/refs/heads/main/21_4_results/ligand_predictions/{system_group}/{csv_name}"
    ref_dg_data = pd.read_csv(csv_file)
    print(f"Extracted reference data for system group {system_group} and system name {system_name} with {len(ref_dg_data)} ligands.")
    # load the ligands with openff toolkit to extract the identifiers
    molecules = Molecule.from_file(input_sdf.as_posix(), allow_undefined_stereo=True)
    molecule_by_name = {molecule.name: molecule for molecule in molecules}
    ref_dataset = {}
    for _, row in ref_dg_data.iterrows():
        ligand_name = row["Ligand name"]
        if ligand_name not in molecule_by_name:
            # this ligand was not included in the benchmark dataset so skip it
            continue
        ligand_data = {
            "dg": row["Exp. dG (kcal/mol)"] * unit.kilocalories_per_mole,
            "canonical_smiles": molecule_by_name[ligand_name].to_smiles(isomeric=True),
            # do we need to be tautomer specific?
            "inchikey": molecule_by_name[ligand_name].to_inchikey(fixed_hydrogens=True),
            "reference": ross_et_al_doi,
        }
        value = row.get("Exp. dG error (kcal/mol)")
        if pd.notna(value):
            dg_uncertainty = row["Exp. dG error (kcal/mol)"] * unit.kilocalories_per_mole
            if dg_uncertainty.m != 0.0:
                ligand_data["uncertainty"] = dg_uncertainty
        ref_dataset[ligand_name] = ligand_data

    # save the extracted reference data to the output directory as json using gufe serialization
    out_file = out_dir / "experimental_binding_data.json"
    with open(out_file, "w") as f:
        json.dump(ref_dataset, f, indent=4, cls=JSON_HANDLER.encoder)
    print(f"Saved extracted reference data to {out_file}")

if __name__ == "__main__":
    main()