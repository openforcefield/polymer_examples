"""
This file is meant to be run with two different versions of the openff toolkit
on a list of pdb files that work/load/perform the same between both. Parameter assignment
data is stored in the "data" folder with the name of the toolkit version used to make
that file and later compared in "compare_data.py"
"""

from openff.toolkit import Topology, Molecule
import openff
import pkg_resources
tk_version = pkg_resources.get_distribution("openff-toolkit").version

import sys
import os
import numpy as np
from pathlib import Path
from rdkit import Chem
from openeye import oechem
import time
import openmm
from openmm.app import PDBFile
from copy import deepcopy
from simtk import unit
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.typing.engines.smirnoff import vdWHandler, AngleHandler
from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, ToolkitRegistry, OpenEyeToolkitWrapper

GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)

current_dir = Path(__file__).parent.resolve()
os.chdir(current_dir)

sys.path.append(os.path.abspath(__file__ + "/../..")) # TODO: fix this mess
from pdb_file_search import PDBFiles

files = ["8gt9",
         "8fy3",
         "8f0x",
         "8e8i",
         "8d1b",
         "8bhw",
         "8ciq",
         "7yb4",
         "1lyd"]
for file_name in files:
    pdb_file = PDBFiles.search(file_name)
    assert pdb_file != None
    top = Topology.from_pdb(str(pdb_file))

    general_offxml = 'openff-2.0.0.offxml'
    amber_offxml = 'ff14sb_off_impropers_0.0.3.offxml'
    water_model = 'tip3p_fb-1.1.0.offxml'
    forcefield = ForceField(general_offxml, amber_offxml) #, water_model)

    file_dir = current_dir / Path("data") / Path(file_name)
    file_dir.mkdir(parents=True, exist_ok=True)

    file_path = file_dir / Path(f"{tk_version}.txt")
    with open(file_path, "w") as file:
        file.write("")

    forcefield.deregister_parameter_handler('ToolkitAM1BCC')
    forcefield.get_parameter_handler('ChargeIncrementModel', {"version":0.3, "partial_charge_method":"gasteiger"})

    system = forcefield.create_interchange(top, allow_nonintegral_charges=True)

    for type, handler in system.collections.items():
        params = handler.get_system_parameters()
        with open(file_path, "a") as file:
            for i in range(params.shape[0]):
                row = ", ".join([str(e) for e in params[i, :]]) + "\n"
                file.write(row)