# Notes:

- The charges and networks are generated using the conda-lock environment specified in `data_generation/conda-lock_linux-64.yml`. Future generation of future partial charges should also be run on the same architecture.
- Charges are generated using the `charge_offmol.py` script in the `data_generation` folder, which accepts either a single SDF file or a directory of SDF files. Molecules that fail charging are skipped and reported. The bulk-charging `charge_molecules.py` script is also available.
- Industry benchmark networks are pulled for relevant systems using the `generate_industry_lomap_networks.py` script in the `data_generation` folder
