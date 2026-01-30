from abc import ABC, abstractmethod
import pathlib
import json

from rdkit import Chem

from openff.units import unit
import openfe
from openfe import (
    SmallMoleculeComponent, SolventComponent, ProteinComponent,
)
from openfe.protocols.openmm_rfe.equil_rfe_methods import RelativeHybridTopologyProtocol

from openfe_benchmarks.data import BenchmarkData
from openfe_benchmarks.scripts import utils as ofebu

class BenchmarkSystem(ABC):
    """
    Abstract base class defining the components of an alchemical network.

    Parameters
    ----------
    bmd_obj : BenchmarkData
        Valid BenchmarkData object
    partial_charge_scheme : str
        Key that corresponds to ``BenchmarkData.ligands``
    solvent : SolventComponent
        OpenFE SolventComponent object
    forcefield : str
        Force field name supported by ``openff.forcefields``, i.e., "openff-2.3.0"
    settings : RelativeHybridTopologyProtocol
        Object defining the settings in the free energy calculation

    Attributes
    ----------
    benchmark_data : BenchmarkData
        Object containing paths to necessary data
    benchmark_set : str
        Fully qualified name of the benchmark set this data system belongs to
        (e.g., 'industry_benchmark_systems.charge_annihilation_set')
    forcefield : str
        Force field name supported by ``openff.forcefields``, i.e., "openff-2.3.0"
    ligand_dict : List[SmallMoleculeComponent]
        List of SmallMoleculeComponent objects for each ligand in the benchmark
        system.
    name : str
        The name / identifier of the benchmark system
    network : LigandNetwork
        Network of SmallMoleculeComponent transformations.
    partial_charge_scheme : str
        Partial charge type used in alchemical network for ligand and cofactors.
    protein : ProteinComponent
        ProteinComponent defining the host molecule of the benchmark system
    settings : RelativeHybridTopologyProtocol
        Object defining the settings in the free energy calculation
    solvent : SolventComponent
        SolventComponent defining the solvent used for the benchmark system
    """
    def __init__(self, bmd_obj: BenchmarkData, partial_charge_scheme: str, solvent: SolventComponent, forcefield, settings=None):
        
        self.solvent = solvent
        self.benchmark_data = bmd_obj
        self.name = bmd_obj.name
        self.benchmark_set = bmd_obj.benchmark_set
        self.partial_charge_scheme = partial_charge_scheme
        self.forcefield = forcefield
        self.settings = self._settings(forcefield) if settings is None else settings
        
        self._process_inputs(bmd_obj.ligands, bmd_obj.protein, bmd_obj.cofactors)
        self.initial_network = openfe.LigandNetwork.from_json(file=str(bmd_obj.network))
        self.network = None

    @staticmethod
    @abstractmethod
    def _settings(forcefield) -> gufe.Protocol:
        """Return the alchemical topology settings for the benchmark system."""
        pass


    def _process_inputs(self, ligands, protein, cofactors):
        """Process ligands, proteins, and cofactors if required and provide warning if not. """
        self._process_ligands(ligands)
        self._process_protein(protein)
        self._process_cofactors(cofactors)


    def _process_ligands(self, ligands):
        """Process ligands if required and provide warning if not. """
        self.ligand_dict = ofebu.process_sdf(ligands[self.partial_charge_scheme], return_dict=True)


    def _process_protein(self, protein):
        """Process protein if required and provide warning if not. """
        self.protein = None if protein is None else ProteinComponent.from_pdb_file(str(protein), name=protein.name)


    def _process_cofactors(self, cofactors):
        """Process cofactors if required and provide warning if not. """
        self.cofactors = None if cofactors is None else ofebu.process_sdf(cofactors[self.partial_charge_scheme])

    
    def generate_alchemical_network(self):
        """
        Generate an alchemical network based on the initial network, settings, 
        solvent, ligand dictionary, protein, and cofactors.
        This method compiles the transformations required for the alchemical 
        network using the specified parameters and creates an 
        `openfe.AlchemicalNetwork` object.
        Returns
        -------
        None
            The generated alchemical network is stored in the `self.network` attribute.
        Notes
        -----
        - The `self.initial_network` should define the initial state of the network.
        - The `self.settings` specifies the settings to be used for the transformations.
        - The `self.solvent`, `self.ligand_dict`, `self.protein`, and `self.cofactors` 
          are used to define the environment and components of the network.
        """
        
        transformations = ofebu.compile_network_transformations(
            self.initial_network, self.settings, self.solvent, self.ligand_dict, self.protein, self.cofactors
        )
        self.network = openfe.AlchemicalNetwork(edges=transformations)
    
    def export_alchemical_network(self, filename="alchemical_network.json"):
        """
        Export the alchemical network to a JSON file.

        This method serializes the alchemical network into a JSON format and saves 
        it to the specified file.

        Parameters
        ----------
        filename : str, optional
            The name of the file to which the alchemical network will be exported. 
            Defaults to "alchemical_network.json".

        Notes
        -----
        The serialization process uses the `gufe.tokenization.JSON_HANDLER.encoder` 
        to encode the network data into JSON format.
        """
        """"""
        
        self.network.to_json(file=filename)
        
    def print_alchemical_network(self):
        """
        Pretty print the alchemical network.

        This method prints the alchemical network in a human-readable JSON format.
        """
        print(json.dumps(self.network.to_dict(), indent=4))

    
class RBFEBenchmarkSystem(BenchmarkSystem):
    """
    Class defining the components and alchemical network of a relative free
    energy benchmark system.

    Parameters
    ----------
    bmd_obj : BenchmarkData
        Valid BenchmarkData object
    partial_charge_scheme : str
        Key that corresponds to ``BenchmarkData.ligands``
    solvent : SolventComponent
        SolventComponent defining the solvent used for the benchmark system
    forcefield : str
        Force field name supported by ``openff.forcefields``, i.e., "openff-2.3.0"
    settings : RelativeHybridTopologyProtocol
        Object defining the settings in the free energy calculation

    Attributes
    ----------
    benchmark_data : BenchmarkData
        Object containing paths to necessary data
    benchmark_set : str
        Fully qualified name of the benchmark set this data system belongs to
        (e.g., 'industry_benchmark_systems.charge_annihilation_set')
    forcefield : str
        Force field name supported by ``openff.forcefields``, i.e., "openff-2.3.0"
    ligand_dict : List[SmallMoleculeComponent]
        List of SmallMoleculeComponent objects for each ligand in the benchmark
        system.
    name : str
        The name / identifier of the benchmark system
    network : LigandNetwork
        Network of SmallMoleculeComponent transformations.
    partial_charge_scheme : str
        Partial charge type used in alchemical network for ligand and cofactors.
    protein : ProteinComponent
        ProteinComponent defining the host molecule of the benchmark system
    settings : RelativeHybridTopologyProtocol
        Object defining the settings in the free energy calculation
    solvent : SolventComponent
        SolventComponent defining the solvent used for the benchmark system
    """

    @staticmethod
    def _settings(forcefield) -> RelativeHybridTopologyProtocol:
        """Utility method for getting RFEProtocol settings for non charge changing transformations.
        """
        # Are there additional settings we should specify here?
        settings = RelativeHybridTopologyProtocol.default_settings()
        settings.engine_settings.compute_platform = "CUDA"

        # Set ligand FF
        settings.forcefield_settings.small_molecule_forcefield = forcefield

        # Fast settings
        settings.simulation_settings.time_per_iteration = 2.5 * unit.picosecond
        settings.simulation_settings.real_time_analysis_interval = 1 * unit.nanosecond
        settings.output_settings.checkpoint_interval = 1 * unit.nanosecond
        settings.solvation_settings.box_shape = 'dodecahedron'
        settings.forcefield_settings.nonbonded_cutoff = 0.9 * unit.nanometer
        settings.solvation_settings.solvent_padding = 1.0 * unit.nanometer # for complex
        #settings.solvation_settings.solvent_padding = 1.5 * unit.nanometer # for ligand

        # Only run one repeat per input json file
        settings.protocol_repeats = 1
        
        ## For charge transformations
        #settings = RelativeHybridTopologyProtocol.default_settings()
        #settings.engine_settings.compute_platform = "CUDA"
        ## Should we use this new OpenFF version or the default?
        #settings.forcefield_settings.small_molecule_forcefield = forcefield
        #settings.alchemical_settings.explicit_charge_correction = True
        #settings.simulation_settings.production_length = 20 * unit.nanosecond
        #settings.simulation_settings.n_replicas = 22
        #settings.lambda_settings.lambda_windows = 22
        ## Only run one repeat per input json file
        #settings.protocol_repeats = 1
            
        return settings

    def _process_inputs(self, ligands, protein, cofactors):
        super()._process_inputs(ligands, protein, cofactors)
        if self.protein is None:
            raise ValueError(f"Protein data is missing or could not be processed for {self.benchmark_set} / {self.name}")