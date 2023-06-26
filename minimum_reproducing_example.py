from openff.toolkit import Topology, Molecule
import sys
import os
from pathlib import Path
from monomer_generation.partition import partition
from rdkit import Chem
from openeye import oechem
import openmm
from openmm.app import PDBFile
from copy import deepcopy
from simtk import unit
import json
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.typing.engines.smirnoff import vdWHandler
from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, ToolkitRegistry, OpenEyeToolkitWrapper
from openff.amber_ff_ports.amber_ff_ports import get_forcefield_dirs_paths

GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)

def simulate_polymer(pdbfile, off_topology, forcefield):
    # mol should already have one conformer...

    pdbfile = PDBFile(pdbfile)
    omm_topology = pdbfile.topology

    forcefield.deregister_parameter_handler('ToolkitAM1BCC')
    forcefield.get_parameter_handler('ChargeIncrementModel', {"version":0.3, "partial_charge_method":"gasteiger"})

    system = forcefield.create_openmm_system(off_topology, allow_nonintegral_charges=True)

    time_step = 2*unit.femtoseconds  # simulation timestep
    temperature = 1000*unit.kelvin  # simulation temperature
    friction = 1/unit.picosecond  # collision rate
    integrator = openmm.LangevinIntegrator(temperature, friction, time_step)
    simulation = openmm.app.Simulation(omm_topology, system, integrator)
    positions = pdbfile.getPositions() 
    simulation.context.setPositions(positions)
    
    simulation.minimizeEnergy(maxIterations=1000)
    simulation.step(1)
    st = simulation.context.getState(getPositions=True, getEnergy=True)
    return st.getPotentialEnergy()

def successfully_loaded(top):
    match_info = [atom.metadata["match_info"] for atom in top.atoms]
    return all([bool(match) for match in match_info])

current_dir = Path(__file__).parent.resolve()
os.chdir(current_dir)

pdb_file = "compatible_pdbs/proteins/6cww.pdb"
json_file = "monomer_generation/json_files/6cww.json"

substructs = dict()
with open(json_file, "r") as file:
    substructs = json.load(file)
substructs = substructs["monomers"]

top = Topology.from_pdb(str(pdb_file), _custom_substructures=substructs)
assert successfully_loaded(top)
print("\t sucessfully loaded!")

if partition(top):
    print("\t sucessfully partitioned!")

general_offxml = 'openff-2.0.0.offxml'
amber_offxml = 'ff14sb_off_impropers_0.0.3.offxml'
water_model = 'tip3p_fb-1.1.0.offxml'

forcefield = ForceField(general_offxml, amber_offxml, water_model)
# new_parameter = vdWHandler.vdWType(
#             smirks='[*:1]',
#             epsilon=0.0157*unit.kilocalories_per_mole,
#             rmin_half=0.6000*unit.angstroms,
#         )
# forcefield.get_parameter_handler('vdW').parameters.insert(0, new_parameter)

energy = simulate_polymer(str(pdb_file), top, forcefield)
print("\t sucessfully parameterized!")
print(energy)