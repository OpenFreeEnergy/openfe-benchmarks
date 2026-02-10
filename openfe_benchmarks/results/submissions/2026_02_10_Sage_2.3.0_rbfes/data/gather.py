import pandas as pd
import pathlib
import click
from cinnabar import FEMap
from openff.units import unit


def get_exp_ddg_fep_plus(
    experimental_data: pd.DataFrame, ligand_a: str, ligand_b: str
) -> tuple[float, float]:
    """
    Given the experimental dataframe reported by fep plus and the ligand a & b identifiers get the reference DDG and
    error.

    Parameters
    ----------
    experimental_data: pd.DataFrame
        The experimental dataframe reported by FEP plus
    ligand_a, ligand_b: str
        The identifiers of the ligands involved in the edge, computed as A to B.
    Returns
    -------
        The experimental DDG and error estimate for this transformation.
    """
    # make it clear which result can not be found!
    try:
        ligand_a_dg = experimental_data[
            experimental_data["Ligand name"] == ligand_a
        ].iloc[0]["Exp. dG (kcal/mol)"]
        if "Exp. dG error (kcal/mol)" in experimental_data.columns:
            ligand_a_ddg = experimental_data[
                experimental_data["Ligand name"] == ligand_a
            ].iloc[0]["Exp. dG error (kcal/mol)"]
        else:
            ligand_a_ddg = 0.0
    except IndexError as e:
        print(ligand_a, " not found!")
        raise e
    try:
        ligand_b_dg = experimental_data[
            experimental_data["Ligand name"] == ligand_b
        ].iloc[0]["Exp. dG (kcal/mol)"]
        if "Exp. dG error (kcal/mol)" in experimental_data.columns:
            ligand_b_ddg = experimental_data[
                experimental_data["Ligand name"] == ligand_b
            ].iloc[0]["Exp. dG error (kcal/mol)"]
        else:
            ligand_b_ddg = 0.0
    except IndexError as e:
        print(ligand_b, " not found!")
        raise e
    ddg = ligand_b_dg - ligand_a_dg
    ddg_error = (ligand_a_ddg**2 + ligand_b_ddg**2) ** 0.5
    return ddg, ddg_error


def get_dataset_ddG(template: str) -> pd.DataFrame:
    glob_search = pathlib.Path(".").glob("[!.]*")
    systems = [p for p in glob_search if p.is_dir()]

    all_data = []
    for system in systems:
        exp_data = pd.read_csv(system / f"{system.name}_ref.csv")
        calc_data = pd.read_csv(system / template, sep="\t")

        for row in calc_data.to_dict(orient="records"):
            edge_data = {}
            edge_data["dataset"] = system.name
            edge_data["ligand_A"] = row["ligand_i"]
            edge_data["ligand_B"] = row["ligand_j"]
            exp_ddG, exp_err = get_exp_ddg_fep_plus(
                exp_data, row["ligand_i"], row["ligand_j"]
            )
            edge_data["exp_ddG"] = exp_ddG
            edge_data["exp_err"] = exp_err
            edge_data["calc_ddG"] = row["DDG(i->j) (kcal/mol)"]
            edge_data["calc_err"] = row["uncertainty (kcal/mol)"]
            all_data.append(edge_data)
    return pd.DataFrame(all_data)


def get_dataset_dG(edges_data) -> pd.DataFrame:
    glob_search = pathlib.Path(".").glob("[!.]*")
    systems = [p for p in glob_search if p.is_dir()]

    all_datasets = []
    for system in systems:
        experimental_data = pd.read_csv(system / f"{system.name}_ref.csv")
        edge_data = edges_data[edges_data["dataset"] == system.name]
        fe_map = FEMap()

        for row in edge_data.to_dict(orient="records"):
            fe_map.add_relative_calculation(
                labelA=row["ligand_A"],
                labelB=row["ligand_B"],
                value=row["calc_ddG"] * unit.kilocalorie_per_mole,
                uncertainty=row["calc_err"] * unit.kilocalorie_per_mole,
            )

        for row in experimental_data.to_dict(orient="records"):
            fe_map.add_experimental_measurement(
                label=row["Ligand name"],
                value=row["Exp. dG (kcal/mol)"] * unit.kilocalorie_per_mole,
                uncertainty=0.0 * unit.kilocalorie_per_mole,
            )

        fe_map.generate_absolute_values()

        # get the absolute calculated values
        abs_df = fe_map.get_absolute_dataframe()
        abs_calc = abs_df[abs_df["computational"]].copy(deep=True)
        exp_data = abs_df[not abs_df["computational"]].copy(deep=True)

        # get the mean exp value to shift the dG values
        mean_exp = exp_data["DG (kcal/mol)"].mean()
        abs_calc["DG (kcal/mol)"] += mean_exp

        exp_values, exp_error = [], []
        for _, row in abs_calc.iterrows():
            exp_row = exp_data[exp_data["label"] == row["label"]].iloc[0]
            exp_values.append(exp_row["DG (kcal/mol)"])
            exp_error.append(exp_row["uncertainty (kcal/mol)"])

        abs_calc["exp_dG"] = exp_values
        abs_calc["exp_err"] = exp_error
        abs_calc["dataset"] = system.name
        abs_calc.rename(
            columns={"DG (kcal/mol)": "calc_dG", "uncertainty (kcal/mol)": "calc_err"},
            inplace=True,
        )
        all_datasets.append(abs_calc.drop(columns=["computational", "source"]))

    return pd.concat(all_datasets)


@click.command()
@click.option(
    "-n",
    "--name",
    help="The label name",
    type=click.STRING,
    default=None,
)
def main(name: str):
    ddG = get_dataset_ddG(f"{name}.ddG")
    dG = get_dataset_dG(ddG)

    ddG.to_csv(f"{name}_ddG.csv")
    dG.to_csv(f"{name}_dG.csv")


if __name__ == "__main__":
    main()
