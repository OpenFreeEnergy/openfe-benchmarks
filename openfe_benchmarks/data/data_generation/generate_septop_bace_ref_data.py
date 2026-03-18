import click
import pathlib
from gufe.tokenization import JSON_HANDLER
import pandas as pd
from openff.units import unit
import json
from openff.toolkit import Molecule

@click.command()
@click.option("--out-dir", type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=pathlib.Path), required=True, help="The output dir name")
@click.option("--input-sdf", type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=pathlib.Path), required=False, help="The input sdf file containing the ligands from this network, used to extract identifiers.")
def main(out_dir: pathlib.Path, input_sdf: pathlib.Path):
    """
    Extract the reference data

    Parameters
    ----------
    out_dir : pathlib.Path
        The output dir name where the extracted reference data will be saved.
    input_sdf
        The input sdf file containing the ligands from this network, used to extract identifiers which will be stored in the reference data.
    """
    # tag each entry with the source, currently we link to the aggregated reference data but this should be changed if we refine it in future.
    baumann_et_al_doi = "https://doi.org/10.1021/acs.jctc.3c00282"
    # load the ref data from the industry benchmark results stored on github
    ref_dg_data = pd.read_csv("/Users/hannahbaumann/test_septop/bace/exp.csv")
    # load the ligands with openff toolkit to extract the identifiers
    molecules = Molecule.from_file(input_sdf.as_posix(), allow_undefined_stereo=True)
    molecule_by_name = {molecule.name: molecule for molecule in molecules}
    ref_dataset = {}
    for _, row in ref_dg_data.iterrows():
        print(row) 
        ligand_name = row["# Ligand"]
        if ligand_name not in molecule_by_name:
            # this ligand was not included in the benchmark dataset so skip it
            continue
        ligand_data = {
            "dg": row["exp_DG"] * unit.kilocalories_per_mole,
            "canonical_smiles": molecule_by_name[ligand_name].to_smiles(isomeric=True),
            # do we need to be tautomer specific?
            "inchikey": molecule_by_name[ligand_name].to_inchikey(fixed_hydrogens=True),
            "reference": baumann_et_al_doi,
        }
        dg_uncertainty = row["unc"] * unit.kilocalories_per_mole
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
