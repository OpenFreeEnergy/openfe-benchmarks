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
│       ├── ligands.sdf        # Ligands without charges
│       ├── ligands_<charge_type>.sdf       # Ligands with specified partial charges
│       ├── cofactors.sdf      # (Optional) System cofactors without charges
│       ├── cofactors_<charge_type>.sdf     # (Optional) System cofactors with charges
│       ├── PREPARATION_DETAILS.md          # System preparation documentation
│       └── <network_name>.json             # Network mappings (e.g., industry_benchmarks_network.json)
```

**Supported Partial Charge Types:**
- `antechamber_am1bcc` - AM1-BCC charges calculated with AmberTools 23.6
- `nagl_openff-gnn-am1bcc-1.0.0.pt` - NAGL AM1-BCC charges
- `openeye_am1bcc` - OpenEye AM1-BCC charges
- `openeye_am1bccelf10` - OpenEye AM1-BCC ELF10 charges

**Notes**:
- Antechamber charges generated with ambertools 23.6 cuda_None_nompi_py312hc98840c_10

### Using Industry Benchmark Systems

The module provides a Python API for accessing benchmark systems:

```python
from openfe_benchmarks.data import (
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
systems = list_systems('industry_benchmark_systems.jacs_set')
print(f"Systems in jacs_set: {systems}")

# Load a specific benchmark system
system = get_benchmark_system('industry_benchmark_systems.jacs_set', 'p38')

# Access system components
print(f"Protein: {system.protein}")
print(f"Ligands: {system.ligands}")  # Dict mapping charge type to file path
print(f"Cofactors: {system.cofactors}")  # Dict mapping charge type to file path (may be empty)
print(f"Networks: {system.networks}")  # List of network JSON file paths
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
from openfe import SmallMoleculeComponent, ProteinComponent, LigandNetwork
from rdkit import Chem

# Load protein
protein = ProteinComponent.from_pdb_file('protein.pdb')

# Load ligands
ligands = [SmallMoleculeComponent(mol) 
           for mol in Chem.SDMolSupplier('ligands.sdf', removeHs=False)
           if mol is not None]

# Load network from JSON file
network = LigandNetwork.from_json(file='industry_benchmarks_network.json')
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

#### Serialization - JSON Format

Network files in the benchmark systems use JSON serialization format.

**Load from JSON file:**
```python
# Load network from JSON file (preferred method)
network = LigandNetwork.from_json(file='network.json')

# Or load from JSON string
with open("network.json", "r") as f:
    json_str = f.read()
network = LigandNetwork.from_json(content=json_str)
```

**Save to JSON:**
```python
# Save network to JSON file
with open("network.json", "w") as f:
    f.write(network.to_json())
```

JSON format preserves:
- Ligand structures (as serialized SmallMoleculeComponent objects)
- Atom mappings between ligands
- Edge annotations (scores, metadata, etc.)

**Note:** GraphML format is also supported for LigandNetwork serialization using `.to_graphml()` and `.from_graphml()` methods.

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

**Note:** `AlchemicalNetwork` uses JSON serialization (not GraphML). Use JSON serialization for `AlchemicalNetwork`:

```python
# Serialize to JSON
network_json = alchemical_network.to_json()

# Save to file
with open("alchemical_network.json", "w") as f:
    f.write(network_json)

# Deserialize from JSON
with open("alchemical_network.json", "r") as f:
    alchemical_network = AlchemicalNetwork.from_json(f.read())
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
   - `ligands.sdf` for ligands without charges
   - `ligands_<charge_type>.sdf` for ligands with specific charge types
   - `cofactors.sdf` for cofactors without charges (if present)
   - `cofactors_<charge_type>.sdf` for cofactors with charges (if present)
   - `PREPARATION_DETAILS.md` (required) - system preparation documentation
   - `*network.json` for network files (JSON format)
3. Include network mapping files (`.json` format preferred)
4. Add new charge types to `PARTIAL_CHARGE_TYPES` in `data/__init__.py`
5. Document preparation details in `PREPARATION_DETAILS.md` within the system directory