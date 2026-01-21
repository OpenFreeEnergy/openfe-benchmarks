# Open Free Energy Benchmark Systems

A set of benchmark systems to validate the OpenFE components.

## Free Energy Calculations

This repository is meant to store the inputs and results of free energy calculations in a systematic way so that others can, view, repeat, or build off of previous work.

There are three types of free energy calculations currently considered:
1. Relative Binding Free Energy (RBFE)
2. Absolute Binding Free Energy (ABFE)
3. Solvation Free Energy (SFE)

## Contents
### `submissions` Directory

Each benchmark record of any free energy calculation has three main components:
1. `inputs`: 
    - python modules / jupyter notebooks that create inputs
    - json files for each transformation
    - **Not** *.pdb, *.sdf, or *.graphml files, see [`data` directory](#data-directory) information.
2. `transformation_outputs`:
    - output json file for each input transformation json
    - output directory of analysis for each input transformation
3. `gathered_outputs`:
    - python modules / jupyter notebooks that generate collective results (likely using cinnabar)
    - openfe ddg, dg outputs and plots
    - log files
4. README.md
    - Contains a description of the network
    - Force field used (including partial charge scheme)
5. env.yaml
    - Fully specified environment used to create the network transformations.

If this description sounds foreign to you, consider starting with the [Open Free Energy Tutorials](https://docs.openfree.energy/en/latest/tutorials/index.html).

### `scripts` Directory

Example generate python modules / jupyter notebooks can be found in `???` and adapted to the particular needs of a calculation before including in the `inputs` and `gathered_outputs` of a directory.

### `data` Directory

Notice that the definition of the benchmark `inputs` directory above does not include *.pdb or *.sdf files. In order to reduce duplication information, all *.pdb, *.sdf, and *graphml files are cataloged in the `data` directory. For more information about the data organization and available structures, see [data/README.md](data/README.md).
