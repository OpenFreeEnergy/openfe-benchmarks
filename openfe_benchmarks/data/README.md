# OpenFE Benchmarks Data

This directory contains benchmark system data for OpenFE (Open Free Energy) calculations, including protein structures, ligands, networks, and remediated industry benchmark systems.

## Industry Benchmark Systems

### Structure

The `industry_benchmark_systems/` directory contains remediated benchmark inputs organized as:

```
industry_benchmark_systems/
├── <benchmark_set>/           # e.g., jacs_set, fragments, janssen_bace, mcs_docking_set
│   └── <system_name>/         # e.g., p38, tyk2, mcl1
│       ├── protein.pdb        # Protein structure with cocrystallized waters & ions
│       ├── ligands_<charge_type>.sdf       # Ligands with specified partial charges
│       ├── cofactors_<charge_type>.sdf     # (Optional) System cofactors
│       └── <network_name>.graphml          # Network mappings (e.g., lomap_network.graphml)
```

**Supported Partial Charge Types:**
- `antechamber_am1bcc` - AM1-BCC charges calculated with AmberTools 23.6
- `openeye_elf10` - ELF10 charges from OpenEye toolkit

**Notes**:
- Antechamber charges generated with ambertools 23.6 cuda_None_nompi_py312hc98840c_10

### Using Industry Benchmark Systems

The module provides a Python API for accessing benchmark systems:

```python
from openfe_benchmarks.data.industry_benchmark_systems import (
    get_benchmark_system,
    list_benchmark_sets,
    list_systems,
    BenchmarkSystem,
    PARTIAL_CHARGE_TYPES,
)

# Discover available benchmark sets
benchmark_sets = list_benchmark_sets()
print(f"Available sets: {benchmark_sets}")

# List systems in a benchmark set
systems = list_systems('jacs_set')
print(f"Systems in jacs_set: {systems}")

# Load a specific benchmark system
system = get_benchmark_system('jacs_set', 'p38')

# Access system components
print(f"Protein: {system.protein}")
print(f"Ligands: {system.ligands}")  # Dict mapping charge type to file path
print(f"Cofactors: {system.cofactors}")  # Dict mapping charge type to file path
print(f"Network mappings: {system.mappings}")  # List of graphml file paths
```

For a complete tutorial, see `examples/using_industry_benchmark_systems.ipynb`.

## Inspecting Ligand Structures

View ligands in a Jupyter notebook with RDKit:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Load ligands from SDF file
ligand_rdmols = [m for m in Chem.SDMolSupplier('ligands.sdf', removeHs=False)]

# Generate 2D coordinates for visualization
[AllChem.Compute2DCoords(ligand) for ligand in ligand_rdmols]

# Display as grid
Chem.Draw.MolsToGridImage(ligand_rdmols, molsPerRow=7)
```

### Loading with OpenFE

```python
import openfe
from openfe import SmallMoleculeComponent, ProteinComponent

# Load protein
protein = ProteinComponent.from_pdb_file('protein.pdb')

# Load ligands
from rdkit import Chem
ligands = [SmallMoleculeComponent(mol) 
           for mol in Chem.SDMolSupplier('ligands.sdf', removeHs=False)]
```

## Networks

Networks define the transformations between ligands for relative free energy calculations.

### LigandNetwork

A `LigandNetwork` is a directed graph connecting ligands via atom mappings. Networks are defined in [gufe](https://github.com/openfreeenergy/gufe).

#### Creating Networks

```python
import openfe
from openfe import LigandNetwork, LigandAtomMapping

# Create mappings between ligands
mappings = [
    LigandAtomMapping(ligand1, ligand2, atom_mapping_dict),
    LigandAtomMapping(ligand2, ligand3, atom_mapping_dict),
    # ...
]

# Create network from mappings
network = LigandNetwork(mappings)
```

#### Serialization - GraphML Format

**Save to GraphML:**
```python
# Save network to file
with open("network.graphml", "w") as f:
    f.write(network.to_graphml())
```

**Load from GraphML:**
```python
# Load network from file
with open("network.graphml", "r") as f:
    graphml_str = f.read()

network = LigandNetwork.from_graphml(graphml_str)
```

GraphML is the primary serialization format for `LigandNetwork`. It preserves:
- Ligand structures (as serialized SmallMoleculeComponent objects)
- Atom mappings between ligands
- Edge annotations (scores, metadata, etc.)

#### Creating AlchemicalNetwork from LigandNetwork

Convert a `LigandNetwork` to a full `AlchemicalNetwork` for RBFE calculations:

```python
from openfe import SolventComponent, ProteinComponent

# Define system components
solvent = SolventComponent()  # Default water with ions
protein = ProteinComponent.from_pdb_file('protein.pdb')

# Create AlchemicalNetwork for RBFE
alchemical_network = network.to_rbfe_alchemical_network(
    solvent=solvent,
    protein=protein,
    protocol=protocol,  # Your chosen Protocol
    autoname=True,
    autoname_prefix="rbfe"
)
```

The `AlchemicalNetwork` contains all `Transformation` objects (edges) and `ChemicalSystem` objects (nodes) needed for execution.

### AlchemicalNetwork

An `AlchemicalNetwork` represents a complete simulation campaign with:
- Nodes: `ChemicalSystem` objects (protein + ligand + solvent combinations)
- Edges: `Transformation` objects (the calculations to perform)

**Note:** `AlchemicalNetwork` does not currently support GraphML serialization (only `LigandNetwork` does). Use JSON serialization for `AlchemicalNetwork`:

```python
# Serialize to dict (then JSON)
network_dict = alchemical_network.to_dict()

# Deserialize from dict
alchemical_network = AlchemicalNetwork.from_dict(network_dict)
```

## Additional Resources

- **OpenFE Documentation**: https://docs.openfree.energy
- **GUFE Documentation**: https://gufe.openfree.energy
- **Example Notebooks**: See `examples/` directory for:
  - `using_industry_benchmark_systems.ipynb` - How to use the benchmark systems API
  - `building_networks.ipynb` - Creating ligand networks
  - `ExploringEdges.ipynb` - Working with atom mappings

## Contributing

When adding new benchmark systems:

1. Follow the directory structure in `industry_benchmark_systems/`
2. Name files according to conventions:
   - `protein.pdb` for protein structures
   - `ligands_<charge_type>.sdf` for ligands
   - `cofactors_<charge_type>.sdf` for cofactors (if present)
3. Include network mapping files (`.graphml` format preferred)
4. Add new charge types to `PARTIAL_CHARGE_TYPES` in `industry_benchmark_systems/__init__.py`
5. Document preparation details in `PREPARATION_DETAILS.md` within the system directory