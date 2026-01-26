# Remediated Input Structures

This directory holds remediated benchmark inputs ready for use with the OpenFE toolkit.

## Structure

Data is organized in the following structure:

```bash
  -- <benchmark_set>
      -- <benchmark_system>
          -- protein.pdb # PDB of protein + cocrystalized waters & ions
          -- ligands.sdf / cofactors.sdf # SDFs of ligands and cofactors with no charges
          -- ligands_antechamber_am1bcc.sdf # SDF containing ligands to be transformed with AM1-BCC charges calculated with antechamber at the input conformation 
          -- cofactors_antechamber_am1bcc.sdf # (Optional) SDF containing any system cofactors with AM1-BCC charges calculated with antechamber at the input conformation
          -- industry_benchmarks_network # A lomap style network with Kartograf atom mappings extracted from the industry benchmarking results
```

## Notes:

- Charges are generated using the `charge_molecules.py` script in the `data_generation` folder
- Industry benchmark networks are generated using the `generate_industry_lomap_networks.py` script in the `data_generation` folder
- The charges and networks are generated using the conda-lock environment specified in `data_generation/conda-lock.yml`
