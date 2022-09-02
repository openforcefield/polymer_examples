from pathlib import Path
import os
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import OpenEyeToolkitWrapper, RDKitToolkitWrapper, toolkit_registry
from rdkit import Chem
import openeye
from openeye import oechem
from collections import defaultdict
from openff.toolkit.utils.utils import temporary_cd   # this is really cool btw
from tempfile import TemporaryDirectory
import sys
from openmm.app import PDBFile

if os.getcwd() == '/home/coda3831/openff-workspace':
    os.chdir('polymer_examples')

input_directory = Path("")
output_directory = Path("compatible_pdbs")
Path(output_directory).mkdir(parents=True, exist_ok=True)
if not output_directory.exists() or not output_directory.is_dir() or output_directory == Path(""):
    print("bad output directory")
    sys.exit(0)
log = ""
failed_files = []
#['nucleic_acids/7sb8_dna.pdb', 'simple_polymers/paam_drei_no_wtr.pdb', 'simple_polymers/peg_c35r_no_wtr.pdb', 'simple_polymers/pnipam_drei_no_wtr.pdb', 'simple_polymers/polythiophene.pdb', 'simple_polymers/polyvinylchloride.pdb']
for file in input_directory.glob('**/*.pdb'):
    if str(output_directory) in str(file): # prevents recursive nightmare
        continue
    if file.name == "xlinked.pdb":
        continue
    if file.name != "peg_c35r_no_wtr_modified.pdb":
        continue
    # if str(file) not in ['simple_polymers/polythiophene.pdb', 'simple_polymers/polyvinylchloride.pdb']:
    #     continue
    msg = f"processing {file.name}:\n" 
    print(msg)
    log = log + msg
    mol = Molecule.from_file(str(file), "PDB", allow_undefined_stereo=True, toolkit_registry=OpenEyeToolkitWrapper())
    oemol = mol.to_openeye()
    openeye.OEPerceiveBondOrders(oemol)
    mol = Molecule.from_openeye(oemol, allow_undefined_stereo=True)
    original_n_atoms = mol.n_atoms
    # output to a pdb file and attempt reading by openmm
    with TemporaryDirectory() as tmpdir:
        with temporary_cd(tmpdir):
            mol.to_file("temp.pdb", 'PDB', toolkit_registry=OpenEyeToolkitWrapper())
            rdmol = Chem.MolFromPDBFile("temp.pdb", sanitize=False, removeHs=False)
            if rdmol == None or rdmol.GetNumAtoms() != oemol.NumAtoms():
                msg = f"\t****unable to process****\n" 
                print(msg)
                log = log + msg
                failed_files.append(str(file))
                continue
            Chem.rdmolfiles.MolToPDBFile(rdmol, "temp.pdb")
            pdb = PDBFile("temp.pdb")
            topology = pdb.topology
            new_n_atoms = topology.getNumAtoms()
            if original_n_atoms != new_n_atoms:
                rdmol = Chem.MolFromPDBFile("temp.pdb", sanitize=False, removeHs=False)
                msg = "\tmaking unique atoms\n"
                print(msg)
                log = log + msg
                
                ids = " 123456789abcdefghijklmnopqrstuvwxzyABCDEFGHIJKLMNOPQRSTUVWXYZ!@#$%^&_+()<?>-=:}{][`~/*"
                element_counts = defaultdict(int)
                for atom in rdmol.GetAtoms():
                    ri = atom.GetPDBResidueInfo()
                    name = ri.GetName()
                    element_counts[name] += 1
                    id = element_counts[name] - 1
                    id1 = int(id / 84) + 1
                    id2 = id % 84 + 1
                    if id < 83:
                        str_id = ids[id]
                    else:
                        str_id = ids[id1] + ids[id2]
                    new_name = name[0:2] + str_id + name[(3-1+len(str_id)):] 
                    ri.SetName(new_name)

    Path(output_directory / file.parent).mkdir(parents=True, exist_ok=True)
    output_path = str(output_directory) + "/" + str(file)
    Chem.MolToPDBFile(rdmol, output_path)
    # finally, check to see all programs can read the output
    mol = Molecule.from_file(output_path, "PDB", allow_undefined_stereo=True, toolkit_registry=OpenEyeToolkitWrapper())
    rdmol = Chem.MolFromPDBFile(output_path, sanitize=False, removeHs=False)
    pdb = PDBFile(output_path)
    topology = pdb.topology
    if not (mol.n_atoms == rdmol.GetNumAtoms() == topology.getNumAtoms()):
        msg = f"\tprocessing error. toolkits not equivalent after cleaning:\n" 
        print(msg)
        log = log + msg

    
print("___________RESULTS/LOG___________")
print(log)
print(f"failed files: {failed_files}")