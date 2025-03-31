# Remediated Input Structures

This directory holds remediated benchmark inputs ready for use with the OpenFE toolkit.

## Structure

Data is organized in the following structure:

```bash
  -- <benchmark_set>
      -- <benchmark_system>
          -- protein.pdb # PDB of protein + cocrystalized waters & ions
          -- ligands_antechamber_am1bcc.sdf # SDF containing ligands to be transformed with AM1-BCC charges calculated with antechamber at the input conformation 
          -- cofactors_antechamber_am1bcc.sdf # (Optional) SDF containing any system cofactors with AM1-BCC charges calculated with antechamber at the input conformation
          -- lomap_network.graphml # A lomap network with Kartograf atom mappings
```

## Notes:

- Antechamber charges generated with ambertools 23.6 cuda_None_nompi_py312hc98840c_10
