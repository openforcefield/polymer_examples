"""
Generates jsons files using the new format up-to-date as of 3/14/23
"""

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
    # print(st.getPotentialEnergy())
    # # print(st.getPositions())
    # unitless_positions = []
    # for vec in st.getPositions():
    #     x = vec.x * 10     # please don't let me forget this
    #     y = vec.y * 10     # how to do this in a... better... way 
    #     z = vec.z * 10
    #     unitless_positions.append([x, y, z])
    # unitless_positions = numpy.array(unitless_positions)
    return st.getPotentialEnergy()

sys.path.append(os.path.abspath(__file__ + "/../..")) # TODO: fix this mess
from pdb_file_search import PDBFiles

def successfully_loaded(top):
    match_info = [atom.metadata["match_info"] for atom in top.atoms]
    return all([bool(match) for match in match_info])

# Make a file to store new jsons (TODO: change this to any new file structure)
current_dir = Path(__file__).parent.resolve()
json_dir = current_dir / Path("json_files")
json_dir.mkdir(parents=False, exist_ok=True)

# set flag if the script should try to test_load the new json file
test_load = True

# create object for json creation and loading:
with open("polymer_energies.txt", "w") as file:
    file.write("\n")

for file_name, monomer_info in ALL_SMILES_INPUT.items():
    engine = SubstructureGenerator()
    for name, substructure_and_caps in monomer_info.items():
        smarts, caps = substructure_and_caps
        
        if caps:
            engine.add_monomer_as_smarts_fragment(smarts, name, caps)
        else:
            engine.add_monomer(name, smarts)
    json_file = json_dir / Path(file_name + ".json")
    engine.output_monomer_info_json(json_file)

    pdb_file = PDBFiles.search(file_name)
    if test_load:
        print(f"testing {file_name}")
        assert pdb_file != None
        substructs = engine.get_monomer_info_dict()["monomers"]
        top = Topology.from_pdb(str(pdb_file), _custom_substructures=substructs)
        assert successfully_loaded(top)

        new_top = Topology()
        for mol in top.molecules:
            smiles = mol.to_smiles()
            new_top.add_molecule(Molecule.from_smiles(smiles))
        top = deepcopy(new_top)

        general_offxml = '/home/coda3831/openff-workspace/openff-forcefields/openforcefields/offxml/openff-2.0.0.offxml'
        amber_offxml = '/home/coda3831/openff-workspace/openff-amber-ff-ports/openff/amber_ff_ports/offxml/ff14sb_off_impropers_0.0.3.offxml'
        forcefield = ForceField(general_offxml, amber_offxml)
        
        new_parameter = vdWHandler.vdWType(
            smirks='[*:1]',
            epsilon=0.0157*unit.kilocalories_per_mole,
            rmin_half=0.6000*unit.angstroms,
        )
        forcefield.get_parameter_handler('vdW').parameters.insert(0, new_parameter)

        energy = simulate_polymer(str(pdb_file), top, forcefield)
        print(energy)
        with open("polymer_energies.txt", "a") as file:
            file.write(f"{file_name}, {energy}\n")

        # # if partition(top):
        # #     print("\t sucessfully partitioned!")

        # # try loading without explicit or implicit Hs
        # ofs = oechem.oemolostream("temp.pdb")
        # flavor = oechem.OEOFlavor_PDB_BONDS | oechem.OEOFlavor_PDB_HETBONDS
        # ofs.SetFlavor(oechem.OEFormat_PDB, flavor)

        # for mol in top.molecules:
        #     for atom in mol.atoms:
        #         atom.stereochemistry = None
        #     oemol = mol.to_openeye()
        #     res_info = []
        #     for atom in oemol.GetAtoms():
        #         r = oechem.OEAtomGetResidue(atom)
        #         res_info.append(tuple([r.GetName(), r.GetResidueNumber(), atom.GetName(), r.GetChainID()]))
        #         if atom.GetAtomicNum() == 1: # remove Hs
        #             oemol.DeleteAtom(atom)
        #     for res_i, atom in zip(res_info, oemol.GetAtoms()):
        #         rname, res_num, aname, chainid = res_i
        #         r = oechem.OEAtomGetResidue(atom)
        #         r.SetName(rname)
        #         r.SetResidueNumber(res_num)
        #         atom.SetName(aname)
        #         r.SetChainID(chainid)

        #     oechem.OEWriteMolecule(ofs, oemol)

        # new_substructs = {}
        # for name, s in substructs.items():
        #     rdmol = Chem.MolFromSmarts(s)
            
        #     for atom in rdmol.GetAtoms():
        #         if atom.GetAtomicNum() > 0:
        #             a_num = atom.GetAtomicNum()
        #             # D_num = len([0 for _ in atom.GetBonds()])
        #             D_num = atom.GetDegree()
        #             F_num = atom.GetFormalCharge()
        #             H_num = len([a for a in atom.GetNeighbors() if a.GetAtomicNum() == 1])
        #             query_string = f"[#{a_num}D{D_num}H{H_num}{F_num:+}:{atom.GetAtomMapNum()}]"
        #             query = Chem.AtomFromSmarts(query_string)
        #             atom.SetQuery(query)
        #         elif atom.GetAtomicNum() == 0:
        #             query_string = f"[*:{atom.GetAtomMapNum()}]"
        #             query = Chem.AtomFromSmarts(query_string)
        #             atom.SetQuery(query)
        #         else:
        #             raise Exception
        #     rdmol = Chem.RemoveAllHs(rdmol)
        #     smarts_string = Chem.MolToSmarts(rdmol)
        #     smarts_string = smarts_string.replace('&', '')

        #     new_substructs[name] = smarts_string
        # top = Topology.from_pdb("temp.pdb", _custom_substructures=new_substructs)
        # assert successfully_loaded(top)
        # print("successfully loaded without explicit Hs")