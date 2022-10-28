from pathlib import Path
import sys, os
from pygments import highlight
from rdkit import Chem
from random import randint
from copy import deepcopy
from openff.toolkit.topology.molecule import FrozenMolecule, Molecule
# from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
from rdkit import Chem
import openmm
from openmm.app import PDBFile, Element
from simtk import unit
from openff.toolkit.topology import Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
import numpy
import time
import parmed

def simulate_polymer(pdbfile, substructure_file, offxml_file, output):
    # mol should already have one conformer...

    # mol_conf, _ = Molecule.from_pdb_and_monomer_info(pdbfile, substructure_file)
    omm_pdbfile = PDBFile(pdbfile)
    omm_topology = omm_pdbfile.topology

    off_topology, _, error = Topology.from_pdb_and_monomer_info(pdbfile, substructure_file, strict=True)
    if error:
        print("unmatched atoms in topology: estimating as single bonds")
    forcefield = ForceField(offxml_file)
    # forcefield.deregister_parameter_handler('ToolkitAM1BCC')
    # get better charges 
    # forcefield.get_parameter_handler('ChargeIncrementModel', {"version":0.3, "partial_charge_method":"gasteiger"})
    start = time.time()
    system = forcefield.create_openmm_system(off_topology, allow_nonintegral_charges=True) #WARNING: I have no idea that this means 
    end = time.time()
    difference = end - start
    time_step = 2*unit.femtoseconds  # simulation timestep
    temperature = 1000*unit.kelvin  # simulation temperature
    friction = 1/unit.picosecond  # collision rate
    integrator = openmm.LangevinIntegrator(temperature, friction, time_step)
    simulation = openmm.app.Simulation(omm_topology, system, integrator)
    positions = omm_pdbfile.getPositions() 
    simulation.context.setPositions(positions)

    pdb_reporter = openmm.app.PDBReporter(f'{output}.pdb', 10)
    dcd_reporter = openmm.app.DCDReporter(f'{output}.dcd', 10)
    simulation.reporters.append(pdb_reporter)
    simulation.reporters.append(dcd_reporter)
    
    simulation.minimizeEnergy(maxIterations=10000)
    simulation.step(1000)
    st = simulation.context.getState(getPositions=True, getEnergy=True)
    print(st.getPotentialEnergy())
    # print(st.getPositions())
    unitless_positions = []
    for vec in st.getPositions():
        x = vec.x * 10     # please don't let me forget this
        y = vec.y * 10     # how to do this in a... better... way 
        z = vec.z * 10
        unitless_positions.append([x, y, z])
    unitless_positions = numpy.array(unitless_positions)
    return st, difference

if __name__ == "__main__":
    os.chdir("openff-workspace/polymer_examples")
    name = "PEO_PLGA"
    pdb_file = None
    json_file = None
    for file in Path(Path.cwd() / Path('compatible_pdbs')).glob("**/*.pdb"):
        if file.stem != name:
            continue
        print(file.name)
        pdb_file = file.absolute()
        json_file = file.parent / Path(f"{name}.json")
    if pdb_file == None or not pdb_file.exists():
        print(f"could not find pdb file: {pdb_file}")
        sys.exit(0)
    if json_file == None or not json_file.exists():
        print(f"could not find json file file {json_file}")
        sys.exit(0)

    offxml_file = 'openff_unconstrained_with_library_charges-2.0.0.offxml'
    st, diff = simulate_polymer(str(pdb_file), json_file, offxml_file, name)
    print(f"time to create openmm system: {diff}")

