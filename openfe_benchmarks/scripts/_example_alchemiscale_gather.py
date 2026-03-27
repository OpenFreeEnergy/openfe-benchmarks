import click
import pathlib
from gufe.archival import AlchemicalArchive
from gufe import AlchemicalNetwork
from alchemiscale import AlchemiscaleClient, ScopedKey, Scope
import logging
import bz2

logger = logging.getLogger(__name__)

def _configure_example_logging(level=logging.INFO):
    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(name)s %(levelname)s: %(message)s")
    )
    # Attach to this module's logger so output appears when running the example
    logger.addHandler(handler)
    logger.setLevel(level)
    # Optionally enable package-wide logs:
    logging.getLogger("openfe_benchmarks").setLevel(level)

@click.command()
@click.option("--network", type=click.Path(dir_okay=False, file_okay=True, exists=True, path_type=pathlib.Path), help="Path to the submitted alchemical network with the scope key saved as the name for which an archive should be made.", default=None)
@click.option("--org", type=str, help="Organization name used to build the scope to check can not be used with --inputs.", default=None)
@click.option("--campaign", type=str, help="Campaign name used to build the scope to check can not be used with --inputs.", default=None)
@click.option("--project", type=str, help="Project name used to build the scope to check can not be used with --inputs.", default=None)
@click.option("--allow-partial/--no-allow-partial", is_flag=True, default=False, help="Allow archive creation for an incomplete network.")
@click.option("--output", type=click.Path(dir_okay=True, path_type=pathlib.Path), help="Path to the output directory where the archive(s) should be made.", required=True)
def main(network, org, campaign, project, allow_partial, output):
    """
    Create an AlchemicalArchive for a given alchemical network which has been submitted to Alchemiscale.
    The network can be specified either by providing the path to the submitted network JSON file or by specifying the
    scope (org, campaign, project). Archives for all networks
    found within that scope will be created.

    Parameters
    ----------
    network : pathlib.Path optional
        The path to the submitted alchemicalnetwork JSON file which should have the network key set as the name
    org : str optional
        The organization name used to build the scope to check can not be used with --network.
    campaign : str optional
        The campaign name used to build the scope to check can not be used with --network.
    project : str optional
        The project name used to build the scope to check can not be used with --network.
    allow_partial : bool
        Whether to allow archive creation for an incomplete network, default is False which means only networks with all edges completed will be archived.
    output : pathlib.Path
        Path to the output directory where the archive(s) should be made.

    Notes
    -----
    - The archive will be saved to the output folder with the same name as the network key in compressed JSON format.

    Examples
    --------
    - Create an archive for a specific network specified by scope
    ```
    python _example_alchemiscale_gather.py --org openfe --campaign my_first_run --project tyk2 --output local_archives
    ```
    - Create an archive for each of the networks in a given campaign
    ```
    python _example_alchemiscale_gather.py --org openfe --campaign my_first_run --output local_archives
    ```
    """
    # we can not use scope with inputs
    if network is not None and (org is not None or campaign is not None or project is not None):
        raise ValueError(
            "Scope searching (org, campaign, project) can not be used when --network is provided these are mutually exclusive.")
    output.mkdir(exist_ok=True, parents=True)
    _configure_example_logging()
    client = AlchemiscaleClient(api_url="https://api.alchemiscale.org")

    networks_by_keys = {}
    # load the local network if provided
    if network is not None:
        logger.info(f"Loading local network")
        network = AlchemicalNetwork.from_json(network.as_posix())
        network_key = ScopedKey.from_str(network.name)
        networks_by_keys[network_key] = network
    else:
        # use the scope to query for network keys
        query_scope = Scope(org=org, campaign=campaign, project=project)
        logger.info(f"Querying Alchemiscale for submissions with scope: {query_scope}")
        network_keys = client.query_networks(scope=query_scope)

        logger.info(f"Number of networks found: {len(network_keys)}")
        for network_key in network_keys:
            logger.info(f"Loading network with key: {network_key}")
            network = client.get_network(network_key, visualize=False)
            networks_by_keys[network_key] = network

    logger.info(f"Creating archives")
    for network_key, network in networks_by_keys.items():
        logger.info(f"Creating archive for network: {network_key}")
        # first check the status
        status = client.get_network_status(network_key)
        status.pop("complete")
        if any(v == 0 for v in status.values()):
            if not allow_partial:
                raise ValueError(f"Network {network_key} is incomplete to archive run with --allow_partial.")
            else:
                logger.warning(f"Network {network_key} has incomplete edges but will be archived anyway as --allow-partial is set.")


        # we need the transformations and the results to build the archive
        transformation_results = []
        alchemiscale_network_results = client.get_network_results(network_key, return_protocoldagresults=True).items()
        for transformation_key, raw_results in alchemiscale_network_results:
            transformation = client.get_transformation(transformation_key, visualize=False)
            transformation_results.append((transformation, raw_results))

        archive = AlchemicalArchive(
            network=network,
            transformation_results=transformation_results,
        )
        # write the archive to a compressed format
        with bz2.open(output / f"{network_key}.json.bz2", "wb") as outfile:
            outfile.write(archive.to_json().encode())

        logger.info(f"Archive for network {network_key} created and saved to {output / f'{network_key}.json.bz2'}")

if __name__ == "__main__":
    main()
