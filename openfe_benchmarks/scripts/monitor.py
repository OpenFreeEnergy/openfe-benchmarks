import click
import os
import pathlib
from alchemiscale import AlchemiscaleClient


@click.command
@click.option(
    '--restart',
    is_flag=True,
    help="Restart any failures",
)
def run(
    restart: bool,
):
    """
    Monitor a network of transformation tasks and
    optionally restart failed tasks.

    Parameters
    ----------
    restart : bool
      Whether or not to attempt to restart failed tasks.
    """
    # Get the alchemiscale bits
    user_id = os.environ['ALCHEMISCALE_ID']
    user_key = os.environ['ALCHEMISCALE_KEY']
    asc = AlchemiscaleClient(
        'https://api.alchemiscale.org',
        user_id,
        user_key
    )

    keys = [p for p in pathlib.Path('.').glob('*.key')]

    for key in keys:
        with open(key, 'r') as f:
            network_sk = f.read()

        asc.get_network_status(network_sk)

        if restart:
            err_tasks = asc.get_network_tasks(network_sk, status="error")
            print(f"Number of errored tasks found: {len(err_tasks)}")
            if len(err_tasks) > 0:
                print("Will attempt to restart tasks")
                asc.set_tasks_status(err_tasks, 'waiting')

                print("Printing new network status")
                asc.get_network_status(network_sk)
            else:
                print("No errored tasks were found, no further action")


if __name__ == "__main__":
    run()
