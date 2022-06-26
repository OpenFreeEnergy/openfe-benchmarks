import os
import copy
from openfe_benchmarks import tyk2
from openfe.setup import ChemicalSystem
from openfe.setup.methods.openmm.equil_rbfe_methods import RelativeLigandTransform
from openff.units import unit

# Set the cuda device
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"]="0"

# Get systems and extract edge
tyk2_system = tyk2.get_system()

target_edges = ("ligand_23", "ligand_27")

for entry in tyk2_system.ligand_network.edges:
    if (entry.molA.name in target_edges and
        entry.molB.name in target_edges):
            edge = entry
            
# Set settings
settings = RelativeLigandTransform.get_default_settings()
settings.simulation_settings.equilibration_length = 1000 * unit.picosecond
settings.simulation_settings.production_length = 5000 * unit.picosecond
settings.system_settings.hydrogen_mass = 3.0
settings.integrator_settings.timestep = 4.0 * unit.femtosecond
settings.integrator_settings.n_steps = 250 * unit.timestep

# Create solvent state
stateA = ChemicalSystem({'ligand': edge.molA,
                         'solvent': tyk2_system.solvent_component})
stateB = ChemicalSystem({'ligand': edge.molB,
                         'solvent': tyk2_system.solvent_component})

simset = copy.deepcopy(settings)
simset.simulation_settings.output_filename = f"{edge.molA.name}_{edge.molB.name}_solvent.nc"
simset.simulation_settings.checkpoint_storage = f"{edge.molA.name}_{edge.molB.name}_solvent_checkpoint.nc"

solvent_transform = RelativeLigandTransform(
    stateA=stateA, stateB=stateB, ligandmapping=edge, settings=simset
)

# Create complex state
stateA = ChemicalSystem({'ligand': edge.molA,
                         'solvent': tyk2_system.solvent_component,
                         'protein': tyk2_system.protein_component})
stateB = ChemicalSystem({'ligand': edge.molB,
                         'solvent': tyk2_system.solvent_component,
                         'protein': tyk2_system.protein_component})

simset = copy.deepcopy(settings)

simset.simulation_settings.output_filename = f"{edge.molA.name}_{edge.molB.name}_complex.nc"
simset.simulation_settings.checkpoint_storage = f"{edge.molA.name}_{edge.molB.name}_complex_checkpoint.nc"

complex_transform = RelativeLigandTransform(
    stateA=stateA, stateB=stateB, ligandmapping=edge, settings=simset
)

# run complex
complex_transform.run(verbose=True)

# run solvent
solvent_transform.run(verbose=True)
