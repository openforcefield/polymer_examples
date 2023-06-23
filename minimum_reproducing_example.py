from openff.toolkit import Topology, Molecule
from substructure_generator import SubstructureGenerator
import sys
import os
from monomer_smiles_input import ALL_SMILES_INPUT
from pathlib import Path
from partition import partition
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


pdb_file = PDBFiles.search(file_name)
