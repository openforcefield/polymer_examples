'''Charge Averaging functions - made into .py file to allow for debugging'''
import csv, json
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass, field

from rdkit import Chem
from openff.units import unit
from openff.interchange import Interchange

from openff.toolkit.topology import Topology
from openff.toolkit.topology.molecule import FrozenMolecule, Molecule, Atom
from openff.toolkit.utils import toolkit_registry
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper, AmberToolsToolkitWrapper
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.typing.engines.smirnoff import parameters as offtk_parameters

from openmm import LangevinMiddleIntegrator
from openmm.app import Simulation, PDBReporter, StateDataReporter
from openmm.unit import kelvin, picosecond, picoseconds, nanometer # need to do some unit conversion with both packages

# Charge calculation methods
def search_mol_files(filename : str, parent_path : Path=Path.cwd()/'compatible_pdbs', extensions : tuple[str, ...]=('pdb', 'json')) -> dict[str, Path]:
    '''Search file tree for a pdb and monomer file with matching names'''
    mol_files = {
        ext : path
            for path in parent_path.glob('**/*.*')
                for ext in extensions
                    if path.name == f'{filename}.{ext}'
    }

    for ext in extensions:
        if ext not in mol_files:
            raise FileNotFoundError(f'Could not find a(n) {ext} file \"{filename}.{ext}\"')
    else:
        return mol_files

def load_mol_and_topo(pdb_path : Path, json_path : Path, verbose : bool=False):
    '''Load Molecule and Topology from a pdb and a monomer json file, performing residue matching on monomer units
    Assumes that the pdb only contains has ONE simple homopolymer (will only load first molecule if multiple are present'''
    off_topology, _, error = Topology.from_pdb_and_monomer_info(str(pdb_path), json_path, strict=True, verbose=verbose)
    mol = next(off_topology.molecules) # get the first molecule (assumed to be the polymer of interest)

    return mol, off_topology

def poll_and_count_molecules(pdb_folder : Path, outname : str=None) -> dict[str, int]:
    '''Takes a path to a folder containing multiple .pdb files and produces
    a csv listing all found molecules and how many atoms each contains'''
    mol_sizes = {}
    mol_names = {path.stem for path in pdb_folder.iterdir()}
    for name in mol_names:
        try:
            mol_files = search_mol_files(name)
            mol, topology = load_mol_and_topo(mol_files['pdb'], mol_files['json'])  # will raise exception if files for molecule are not found
            mol_sizes[name] = len(mol.atoms)
        except FileNotFoundError:
            pass

    if outname is not None: # also write to file if a name for the output is specified
        outpath = pdb_folder/f'{outname}.csv'
        outpath.touch()

        with outpath.open('w') as mol_file:
            writer = csv.writer(mol_file, delimiter=',')
            writer.writerow(['Molecule Name', '# Atoms']) # add columns headers
            for mol_name, mol_size in mol_sizes.items():
                writer.writerow([mol_name, mol_size])

    return mol_sizes

def generate_molecule_charges(mol : Molecule, toolkit_method : str='openeye', partial_charge_method : str='am1bcc') -> Molecule:
    '''Takes a Molecule object and computes partial charges with AM1BCC using toolkit method of choice. Returns charged molecule'''
    toolkits = {
        'rdkit' : RDKitToolkitWrapper,
        'openeye' : OpenEyeToolkitWrapper,
        'ambertools' : AmberToolsToolkitWrapper
    }

    mol.assign_partial_charges( # finally, assign partial charges using those 10 conformers generated 
        partial_charge_method=partial_charge_method, 
        toolkit_registry=toolkits.get(toolkit_method)()
    )
    charged_mol = mol # rename for clarity
    # get some conformers to run elf10 charge method. By default, `mol.assign_partial_charges`
    # uses 500 conformers, but we can generate and use 10 here for demonstration
    # charged_mol.generate_conformers(
    #     n_conformers=10,
    #     rms_cutoff=0.25 * unit.angstrom,
    #     make_carboxylic_acids_cis=True,
    #     toolkit_registry=RDKitToolkitWrapper()
    # ) # very slow for large polymers! 

    print(f'final molecular charges: {charged_mol.partial_charges}')
    # note: the charged_mol has metadata about which monomers were assigned where as a result of the chemicaly info assignment.
    # This can be a way to break up the molecule into repeating sections to partition the library charges 
    for atom in charged_mol.atoms:
        assert(atom.metadata['already_matched'] == True)
        # print(atom.metadata['residue_name'])
    
    return charged_mol # code for exact how thely above function works can be found in openff/toolkit/utils/openeye_wrapper.py under the assign_partial_charges()

# charge averaging methods
AveragedChargeMap = defaultdict[str, dict[int, float]] # makes typehinting clearer

@dataclass
class Accumulator:
    '''Compact container for accumulating averages'''
    sum : float = 0.0
    count : int = 0

    @property
    def average(self) -> float:
        return self.sum / self.count

@dataclass
class AvgResidueCharges:
    '''Dataclass for more conveniently storing averaged charges for a residue group'''
    charges : dict[int, float]
    residue_name : str
    SMARTS : str

def find_repr_residues(cmol : Molecule) -> dict[str, int]:
    '''Determine names and smallest residue numbers of all unique residues in charged molecule
    Used as representatives for generating labelled SMARTS strings '''
    rep_res_nums = defaultdict(set) # numbers of representative groups for each unique residue, used to build SMARTS strings
    for atom in cmol.atoms: 
        rep_res_nums[atom.metadata['residue_name']].add(atom.metadata['residue_number']) # collect unique residue numbers

    for res_name, ids in rep_res_nums.items():
        rep_res_nums[res_name] = min(ids) # choose group with smallest id of each residue to denote representative group

    return rep_res_nums

def get_averaged_charges(cmol : Molecule, index_by : str='SMARTS') -> list[AvgResidueCharges]:
    '''Takes a charged molecule and averages charges for each repeating residue. 
    Returns a list of AvgResidueCharge objects each of which holds:
        - A dict of the averaged charges by atom 
        - The name of the residue associated with the charges
        - A SMARTS string of the residue's structure'''
    rdmol = cmol.to_rdkit() # create rdkit representation of Molecule to allow for SMARTS generation
    rep_res_nums = find_repr_residues(cmol) # determine ids of representatives of each unique residue

    atom_ids_for_SMARTS = defaultdict(list)
    res_charge_accums   = defaultdict(lambda : defaultdict(Accumulator))
    for atom in cmol.atoms: # accumulate counts and charge values across matching substructures
        res_name, substruct_id, atom_id = atom.metadata['residue_name'], atom.metadata['substructure_id'], atom.metadata['pdb_atom_id']
        if atom.metadata['residue_number'] == rep_res_nums[res_name]: # if atom is member of representative group for any residue...
            atom_ids_for_SMARTS[res_name].append(atom_id)             # ...collect pdb id...
            rdmol.GetAtomWithIdx(atom_id).SetAtomMapNum(substruct_id) # ...and set atom number for labelling in SMARTS string

        curr_accum = res_charge_accums[res_name][substruct_id] # accumulate charge info for averaging
        curr_accum.sum += atom.partial_charge.magnitude # eschew units (easier to handle, added back when writing to XML)
        curr_accum.count += 1

    avg_charges_by_residue = []
    for res_name, charge_map in res_charge_accums.items():
        SMARTS = Chem.rdmolfiles.MolFragmentToSmarts(rdmol, atomsToUse=atom_ids_for_SMARTS[res_name]) # determine SMARTS for the current residue's representative group
        charge_map = {substruct_id : accum.average for substruct_id, accum in charge_map.items()}
        charge_container = AvgResidueCharges(charges=charge_map, residue_name=res_name, SMARTS=SMARTS)
        avg_charges_by_residue.append(charge_container)

    return avg_charges_by_residue

def write_new_library_charges(avgs : list[AvgResidueCharges], offxml_src : Path, output_path : Path) -> tuple[ForceField, list[offtk_parameters.LibraryChargeHandler]]:
    '''Takes dict of residue-averaged charges to generate and append library charges to an .offxml file of choice, creating a new xml with the specified filename'''
    assert(output_path.suffix == '.offxml') # ensure output path is pointing to correct file type
    forcefield = ForceField(offxml_src)     # simpler to add library charges through forcefield API than to directly write to xml
    lc_handler = forcefield["LibraryCharges"]

    lib_chgs = [] #  all library charges generated from the averaged charges for each residue
    for averaged_res in avgs:
        lc_entry = {f'charge{cid}' : f'{charge} * elementary_charge' for cid, charge in averaged_res.charges.items()} # stringify charges into form usable for library charges
        lc_entry['smirks'] = averaged_res.SMARTS # add SMIRKS string to library charge entry to allow for correct labelling
        lc_params = offtk_parameters.LibraryChargeHandler.LibraryChargeType(allow_cosmetic_attributes=True, **lc_entry) # must enable cosmetic params for general kwarg passing
        
        lc_handler.add_parameter(parameter=lc_params)
        lib_chgs.append(lc_params)   # record library charges for reference
    forcefield.to_file(output_path) # write modified library charges to new xml (avoid overwrites in case of mistakes)
    
    return forcefield, lib_chgs

# OpenMM simulation methods
def create_sim_from_interchange(interchange : Interchange) -> Simulation:
    '''Sets up a Simulation object using topology and force field data as specified by an Interchange object
    Converts topologies and positions to OpenMM format from OpenFF formats (can support GROMACS format too in future)'''
    openmm_sys = interchange.to_openmm(combine_nonbonded_forces=True) 
    openmm_top = interchange.topology.to_openmm()
    openmm_pos = interchange.positions.m_as(unit.nanometer) * nanometer
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.0005*picoseconds)

    simulation = Simulation(openmm_top, openmm_sys, integrator)
    simulation.context.setPositions(openmm_pos)

    return simulation

def run_simulation(simulation : Simulation, output_folder : Path, output_name : str='md_sim', num_steps=1000, record_freq=10) -> None:
    '''Takes a Simulation object, performs energy minimization, and runs simulation for specified number of time steps
    Recording PBD frames and numerical data to file at the specified frequency'''
    folder_name = str(output_folder) # for some reason OpenMM simulations don;t like Path objects (only take strings)

    # for saving pdb frames and reporting state/energy data
    pdb_rep = PDBReporter(f'{folder_name}/{output_name}_frames.pdb', record_freq)  # save frames at the specified interval
    state_rep = StateDataReporter(f'{folder_name}/{output_name}_data.csv', record_freq, step=True, potentialEnergy=True, temperature=True)
    reporters = (pdb_rep, state_rep)

    # minimize and run simulation
    simulation.minimizeEnergy()
    simulation.saveCheckpoint(f'{folder_name}/{output_name}_checkpoint.chk') # save initial minimal state to simplify reloading process
    for rep in reporters:
        simulation.reporters.append(rep) # add any desired reporters to simulaiton for tracking
    simulation.step(num_steps)