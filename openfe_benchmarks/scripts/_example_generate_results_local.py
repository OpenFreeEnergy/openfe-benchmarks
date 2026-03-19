import json
import click
from gufe.tokenization import JSON_HANDLER
from gufe import AlchemicalNetwork, ProteinComponent
from openff.units import unit
import pathlib
from cinnabar import FEMap
from collections import defaultdict
import numpy as np
import zstandard as zstd
from openfecli.commands.gather import _get_names, _get_type

def _get_simulation_key(result: dict) -> tuple[tuple[str, str], str]:
    lig_a_name, lig_b_name = _get_names(result)
    phase = _get_type(result)
    return (lig_a_name, lig_b_name), phase


@click.command()
@click.option(
    "--results_dir", 
    help="The directory containing the transformation results", 
    multiple=True, 
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=pathlib.Path)
)
@click.option(
    "--network", help="The path to the alchemical network JSON file",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=pathlib.Path)
)
@click.option(
    "--output_dir", help="The directory to write the archive to",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=pathlib.Path)
)
def main(results_dir, network, output_dir):
    """
    Gather all the results for the transformations in the network and write the DDG/DG to a json file with units and metadata which can be used 
    for down stream analysis.
    """
    # load the network to make sure all transformations are present in the results
    network = AlchemicalNetwork.from_json(network.as_posix())
    # make a list of transformation keys using (ligand_a.name, ligand_b.name, phase)
    transformations_to_run = set()
    for transformation in network.edges:
        ligand_a_name = transformation.mapping.componentA.name
        ligand_b_name = transformation.mapping.componentB.name
        # get the phase
        if transformation.stateA.contains(ProteinComponent):
            phase = "complex"
        else:
            phase = "solvent"
        transformations_to_run.add((ligand_a_name, ligand_b_name, phase))


    # make a key using the (lig_a.name, lig_b.name)
    raw_results = defaultdict(list)
    for result_dir in results_dir:
        # search for the results json files
        for result_file in result_dir.glob("*.json.zst"):
            with open(result_file, "rb") as f:
                dctx = zstd.ZstdDecompressor()
                with dctx.stream_reader(f) as reader:
                    result = json.load(reader, cls=JSON_HANDLER.decoder)
                    # make a key for this result
                    key, phase = _get_simulation_key(result)
                    raw_results[key].append((phase, result))

    # now loop over the raw results and extract the ddg/dg and metadata
    gathered_results = {"DG": [], "DDG": []}
    # first workout the system group and system name which should be stored on all egdes
    transformation = list(network.edges)[0]
    mapping_annotations = transformation.mapping.annotations
    system_group = mapping_annotations["system_group"]
    system_name = mapping_annotations["system_name"]

    # check that all simulations in the alchemical network have an associated result
    found_results = set()
    for key, results in raw_results.items():
        lig_a_name, lig_b_name = key
        for phase, result in results:
            key = (lig_a_name, lig_b_name, phase)
            if key not in transformations_to_run:
                raise ValueError(f"Found results for transformation {key} which is not in the alchemical network")
            found_results.add((lig_a_name, lig_b_name, phase))

    # now check that we have results for all transformations
    missing_transformations = transformations_to_run - found_results
    if missing_transformations:
        raise ValueError(f"Missing results for transformations: {missing_transformations}")

    # also build a femap so we can get DGs if possible
    fe_map = FEMap()
    for key, results in raw_results.items():
        lig_a_name, lig_b_name = key
        entry_data = {
            "ligand_a": lig_a_name,
            "ligand_b": lig_b_name,
            "system_group": system_group,
            "system_name": system_name,
            # assuming this is the RBFE protocol we should have the same number of repeats for each phase
            "repeats": len(results) // 2,
        }
        # group the results by phase
        complex_results = [result for phase, result in results if phase == "complex"]
        solvent_results = [result for phase, result in results if phase == "solvent"]
        assert len(complex_results) == len(solvent_results), f"Found different number of complex and solvent results for {key}"
        # get the estimated values for each repeat
        complex_data = [result["estimate"].m_as(unit.kilocalories_per_mole) for result in complex_results]
        complex_dg = np.mean(complex_data) * unit.kilocalories_per_mole
        complex_dg_uncertainty = np.std(complex_data) * unit.kilocalories_per_mole
        # get the solvent data
        solvent_data = [result["estimate"].m_as(unit.kilocalories_per_mole) for result in solvent_results]
        solvent_dg = np.mean(solvent_data) * unit.kilocalories_per_mole
        solvent_dg_uncertainty = np.std(solvent_data) * unit.kilocalories_per_mole

        # get the combinded ddg and uncertainty
        entry_data["DDG"] = complex_dg - solvent_dg
        entry_data["DDG_uncertainty"] = np.sqrt(complex_dg_uncertainty**2 + solvent_dg_uncertainty**2)

        # add the raw values for debugging?
        entry_data["DGs_complex"] = complex_data
        entry_data["DGs_solvent"] = solvent_data

        # calculate the smallest off diagonal element of the mabr overlap matrix and the replica mixing matrix for each result
        for phase_results, label in zip([complex_results, solvent_results], ["Complex", "Solvent"]):
            mbar_overlap_elements = []
            replica_mixing_elements = []
            for phase_result in phase_results:
                # get the key for this protocol result
                result_key = list(phase_result["protocol_result"]["data"].keys())[0]
                # calculate the smallest overlap for the mbar overlap matrix
                overlap_matrix = phase_result["protocol_result"]["data"][result_key][0]["outputs"]["unit_mbar_overlap"]["matrix"]
                mbar_overlap_elements.append(np.diagonal(overlap_matrix, offset=1).min())

                # calculate the smallest off diagonal element of the replica mixing matrix
                mixing_matrix = phase_result["protocol_result"]["data"][result_key][0]["outputs"]["replica_exchange_statistics"]["matrix"]
                replica_mixing_elements.append(np.diagonal(mixing_matrix, offset=1).min())

            entry_data[f"{label}_smallest_mbar_overlaps"] = mbar_overlap_elements
            entry_data[f"{label}_smallest_replica_mixing"] = replica_mixing_elements

        gathered_results["DDG"].append(entry_data)

        # also add the DDG data to the FEMAP
        fe_map.add_relative_calculation(
            labelA=lig_a_name,
            labelB=lig_b_name,
            value=entry_data["DDG"],
            uncertainty=entry_data["DDG_uncertainty"]
        )

    # check if the network is connected and we can calculate the DGs
    if fe_map.check_weakly_connected():
        # generate the absolute values for the map centered around zero these should be shifted when comparing with experiment
        fe_map.generate_absolute_values()
        abs_df = fe_map.get_absolute_dataframe()
        for _, row in abs_df.iterrows():
            entry_data = {
                "ligand": row["label"],
                "DG": row["DG (kcal/mol)"] * unit.kilocalories_per_mole,
                "DG_uncertainty": row["uncertainty (kcal/mol)"] * unit.kilocalories_per_mole,
                "system_group": system_group,
                "system_name": system_name,
                "source": row["source"],
            }
            gathered_results["DG"].append(entry_data)
    
    # write out the data to a json file
    output_file = output_dir / "computaional_results.json"
    with open(output_file, "w") as w:
        json.dump(gathered_results, w, cls=JSON_HANDLER.encoder, indent=4)

if __name__ == "__main__":
    main()



        
