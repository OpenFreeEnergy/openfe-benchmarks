import click
import os
import pathlib
from tqdm import tqdm
import openfe
from alchemiscale import AlchemiscaleClient, Scope


def get_network(filename: pathlib.Path) -> openfe.AlchemicalNetwork:
    """
    Read a serialized AlchemicalNetwork

    Parameters
    ----------
    filename : pathlib.Path
      A path to a file with the serialized AlchemicalNetwork

    Returns
    -------
    AlchemicalNetwork
      An AlchemicalNetwork with the desired Transformations.
    """
    return openfe.AlchemicalNetwork.from_json(filename)


@click.command
@click.option(
    '--org_scope',
    type=str,
    required=True,
    help='The organization scope name',
)
@click.option(
    '--scope_name_campaign',
    type=str,
    required=True,
    help='The campaign transformation scope name',
)
@click.option(
    '--base_scope_name_project',
    type=str,
    required=True,
    help='The base project transformation scope name',
)
def run(
    org_scope: str,
    scope_name_campaign: str,
    base_scope_name_project: str,
):
    """
    Submit and action an AlchemicalNetwork on Alchemiscale.

    Parameters
    ----------
    org_scope : str
      The organization Scope name.
    scope_name_campaign : str
      The campaign Scope name.
    base_scope_name_project : str
      The project Scope name.
    """
    # Get the alchemiscale bits
    user_id = os.environ['ALCHEMISCALE_ID']
    user_key = os.environ['ALCHEMISCALE_KEY']
    asc = AlchemiscaleClient(
        'https://api.alchemiscale.org',
        user_id,
        user_key
    )

    for charge_method in ['elf10']:
        forcefield = "openff-2.1.1"
        network_filename = f"AlchemicalNetworks/openff-2.1.1_{charge_method}_alchemical_network.json"
        print(network_filename)
        # Get the alchemical network
        alchemical_network = get_network(network_filename)

        # Set the scope for the transformation
        ffname = forcefield.replace('-', '_').replace('.', '_')
        scope_name_project = f"{ffname}_{charge_method}_{base_scope_name_project}"
        scope = Scope(org_scope, scope_name_campaign, scope_name_project)
        print(f"Scope is set to: {scope}")

        # Create a network and get a scope key
        an_sk = asc.create_network(alchemical_network, scope)

        # store the scoped key
        with open(f"{forcefield}_{charge_method}.key", 'w') as f:
            f.write(str(an_sk))

        # action out tasks
        repeats = 3
        print(f"Actioning {repeats} repeats for {len(alchemical_network.edges)} transformation")
        for transform in tqdm(alchemical_network.edges):
            transform_sk = asc.get_scoped_key(transform, scope)
            tasks = asc.create_tasks(transform_sk, count=repeats)
            asc.action_tasks(tasks, an_sk)
    
        # check what the network status looks like
        asc.get_network_status(an_sk)


if __name__ == "__main__":
    run()
