import io
import json
import pathlib
import zipfile

import click
import pandas as pd
import requests
from gufe.tokenization import JSON_HANDLER
from openff.toolkit import Molecule
from openff.units import unit

# Experimental data from the SepTop SI on Zenodo.
BAUMANN_ET_AL_DOI = "https://doi.org/10.1021/acs.jctc.3c00282"
ZENODO_SI_URL = "https://zenodo.org/records/8066404/files/SI.zip?download=1"
_BACE_CSV_MEMBER = "SI/BACE1/results/BACE_single_repeat.csv"
# the bJ_* series is named lig_* in our networks; map explicitly
_BJ_TO_LIG = {f"bJ_{n:02d}": f"lig_{n:02d}" for n in range(2, 8)}


def _fetch_bace_reference(url: str = ZENODO_SI_URL) -> pd.DataFrame:
    """
    Download the BACE1 experimental data from the Zenodo SI and clean it.

    The CSV has experimental and calculated values separated by
    ``# Experimental block`` / ``# Calculated block`` markers; we slice out the
    experimental rows, strip the ``# `` header prefix, and rename the bJ_*
    series into the lig_* convention used by our networks.
    """
    response = requests.get(url, timeout=300)
    response.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        raw_csv = archive.read(_BACE_CSV_MEMBER).decode()

    lines = raw_csv.splitlines()
    start = lines.index("# Experimental block") + 1
    end = lines.index("# Calculated block")
    block = lines[start:end]
    block[0] = block[0].lstrip("# ")  # "# Ligand,exp_DG,exp_dDG" -> "Ligand,..."
    exp_data = pd.read_csv(io.StringIO("\n".join(block)))
    exp_data["Ligand"] = exp_data["Ligand"].replace(_BJ_TO_LIG)
    return exp_data


def _load_reference_data(ref_data: pathlib.Path | None) -> pd.DataFrame:
    """Read a pre-cleaned reference CSV, or fetch and clean one from Zenodo."""
    if ref_data is not None:
        return pd.read_csv(ref_data)
    return _fetch_bace_reference()


@click.command()
@click.option(
    "--out-dir",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=pathlib.Path),
    required=True,
    help="The output dir where the extracted reference data will be saved.",
)
@click.option(
    "--input-sdf",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=pathlib.Path),
    required=True,
    help="The input sdf file containing the ligands from this network, used to extract identifiers.",
)
@click.option("--verbose", is_flag=True, help="Report network ligands with no matching experimental data.")
def main(out_dir: pathlib.Path, input_sdf: pathlib.Path, verbose: bool):
    """
    Extract experimental reference data.

    Parameters
    ----------
    out_dir
        Directory where the extracted reference data will be written.
    input_sdf
        SDF of the ligands in this network, used to extract identifiers (SMILES,
        InChIKey) that are stored alongside each experimental value.
    """
    # Pull the experimental data from the Zenodo SI
    ref_dg_data = _fetch_bace_reference()

    # load the ligands to extract their identifiers
    molecules = Molecule.from_file(input_sdf.as_posix(), allow_undefined_stereo=True)
    molecule_by_name = {molecule.name: molecule for molecule in molecules}

    ref_dataset = {}
    for _, row in ref_dg_data.iterrows():
        ligand_name = row["Ligand"]
        if ligand_name not in molecule_by_name:
            # this ligand was not included in the benchmark dataset so skip it
            continue
        molecule = molecule_by_name[ligand_name]
        ligand_data = {
            "dg": row["exp_DG"] * unit.kilocalories_per_mole,
            "canonical_smiles": molecule.to_smiles(isomeric=True),
            "inchikey": molecule.to_inchikey(fixed_hydrogens=True),
            "reference": BAUMANN_ET_AL_DOI,
        }
        dg_uncertainty = row["exp_dDG"] * unit.kilocalories_per_mole
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
