from openff.toolkit import Topology, Molecule
import openff
import pkg_resources
tk_version = pkg_resources.get_distribution("openff-toolkit").version

import sys
import os
from pathlib import Path
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

# If both versions of the toolkit load the same topology, the outputted 
# parameters should return the exact values in the same order as loaded by
# interchange. Below, we compare the file contents exactly as an easy check for this. 
for file_name in files:
    file_dir = current_dir / Path("data") / Path(file_name)
    file_contents = ""
    for params_file in file_dir.iterdir():
        with open(params_file, "r") as file:
            new_file_contents = file.read()
            if file_contents:
                print(f"comparing {params_file.name} to previous...", end="")
                assert file_contents == new_file_contents
                print("IDENTICAL")
            else:
                print(f"loading {params_file.name} for the first time...", end="")
                file_contents = new_file_contents