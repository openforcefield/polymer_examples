from pathlib import Path
import os
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import OpenEyeToolkitWrapper, RDKitToolkitWrapper, toolkit_registry
from rdkit import Chem
import openeye
from openeye import oechem
from collections import defaultdict
from openff.toolkit.utils.utils import temporary_cd
from tempfile import TemporaryDirectory
import sys
from openmm.app import PDBFile
from pdbfixer import PDBFixer
import time

os.chdir("polymer_examples")

class Log:
    def __init__(self):
        self.info = ""
        self.failed_files = {}
    def append_info(self, new_info):
        self.info += new_info
    def append_ffile(self, file, reason):
        self.failed_files[file] = reason
    def printout(self):
        print("___________RESULTS/LOG___________")
        print(self.info)
        if self.failed_files:
            "The following files failed to either load, or convert to rdkit"
            for file, reason in self.failed_files.items():
                print(f"\t{file} -> {reason}")
        else:
            print("No files failed, yay!")

cwd = Path(__file__).parent.absolute()
os.chdir(cwd)

input_directory = Path("uncleaned_pdbs")
output_directory = Path("compatible_pdbs")

if not output_directory.exists() or not output_directory.is_dir() or output_directory == Path(""):
    print("bad output directory")
    sys.exit(0)

log = Log()
#['nucleic_acids/7sb8_dna.pdb', 'simple_polymers/paam_drei_no_wtr.pdb', 'simple_polymers/peg_c35r_no_wtr.pdb', 'simple_polymers/pnipam_drei_no_wtr.pdb', 'simple_polymers/polythiophene.pdb', 'simple_polymers/polyvinylchloride.pdb']
skipped_files = ["xlinked.pdb", "6cww.pdb"]

standard_workflow = [
                    #  input_directory / Path("crosslinked_polymers"),
                    #  input_directory / Path("peptoids"),
                    #  input_directory / Path("simple_polymers"),
                    #  input_directory / Path("sugars"),
                    #  input_directory / Path("nucleic_acids"),
                     ]

protein_workflow = [
                    input_directory / Path("proteins")
                    ]
# for directory in standard_workflow:
#     for file in directory.glob('**/*.pdb'):
#         if file.name in skipped_files:
#             continue

#         log.append_info(f"processing {file.name}: ")
#         try:
#             start = time.time()
#             try:
#                 mol = Molecule.from_file(str(file), "PDB", allow_undefined_stereo=True, toolkit_registry=OpenEyeToolkitWrapper())
#             except Exception as e:
#                 log.append_info("openeye failed to read, trying with rdkit. ")
#                 rdmol = Chem.MolFromPDBFile(str(file), removeHs=False, sanitize=False)
#                 mol = Molecule.from_rdkit(rdmol, allow_undefined_stereo=True)

#             oemol = mol.to_openeye()
#             openeye.OEPerceiveBondOrders(oemol)
#             mol = Molecule.from_openeye(oemol, allow_undefined_stereo=True)
#             original_n_atoms = mol.n_atoms

#             rdmol = None
#             with TemporaryDirectory() as tmpdir:
#                 with temporary_cd(tmpdir):
#                     mol.to_file("temp.pdb", 'PDB', toolkit_registry=OpenEyeToolkitWrapper())
#                     rdmol = Chem.MolFromPDBFile("temp.pdb", sanitize=False, removeHs=False)
#                     if rdmol == None:
#                         log.append_ffile(str(file), "rdkit failed to load pdb file processed with openeye")
#                         continue
#                     if rdmol.GetNumAtoms() != oemol.NumAtoms():
#                         log.append_ffile(str(file), "rdkit read a different number of atoms compared to openeye")
#                         continue
#             # make atom and residue names unique
#             element_counts = defaultdict(int)
#             res_number = 1
#             for atom in rdmol.GetAtoms():
#                 ri = atom.GetPDBResidueInfo()
#                 name = ri.GetName()
#                 if len(atom.GetSymbol()) == 1:
#                     new_name = atom.GetSymbol() + f"{element_counts[atom.GetSymbol()]:02d}"
#                 else:
#                     new_name = atom.GetSymbol() + f"{element_counts[atom.GetSymbol()]:01d}"
#                 ri.SetResidueNumber(res_number)
#                 ri.SetResidueName("UNK")
#                 ri.SetName(new_name)

#                 element_counts[atom.GetSymbol()] += 1
#                 if element_counts[atom.GetSymbol()] >= 100 / (10**(len(atom.GetSymbol())-1)): # if any atoms exceed 100 
#                     res_number += 1
#                     element_counts = defaultdict(int)
                


#             relative_file_path = Path(os.path.relpath(file, input_directory))
#             output_path = str(output_directory / relative_file_path)
#             Chem.MolToPDBFile(rdmol, output_path)
#             # finally, check to see all programs can read the output
#             openmm_mol = PDBFile(output_path)
#             rdmol = Chem.MolFromPDBFile(output_path, sanitize=False, removeHs=False)
#             if not (mol.n_atoms == openmm_mol.topology.getNumAtoms()):
#                 log.append_ffile(str(file), "toolkits did not read an equivalent number of atoms in the final pdb")
#             log.append_info(f"{time.time()-start:2f}s\n")
#         except Exception as e:
#             log.append_ffile(str(file), f"Critical error: {e}")

for directory in protein_workflow:
    for file in directory.glob('**/*.cif'):
        if file.name in skipped_files:
            continue
        print(f"processing {file.name}: ")
        log.append_info(f"processing {file.name}: ")
        try:
            start = time.time()
            # use openeye to read the first alternate location from cif file. Save chemical info in
            # sdf file and manually save residue info for any uses later
            ifs = oechem.oemolistream(str(file.absolute()))
            ofs = oechem.oemolostream("temp.sdf")

            flavor = oechem.OEIFlavor_Generic_Default | oechem.OEIFlavor_MMCIF_NoAltLoc 
            ifs.SetFlavor(oechem.OEFormat_MMCIF, flavor)

            residue_info = {}
            for mol in ifs.GetOEGraphMols():
                for atom in mol.GetAtoms():
                    r = oechem.OEAtomGetResidue(atom)
                    residue_info[atom.GetIdx()] = tuple([r.GetName(), r.GetResidueNumber(), atom.GetName(), r.GetChainID()])
                oechem.OEWriteMolecule(ofs, mol)

            # load into rdkit to add explicit hydrogens for nonstandard residuse. 
            # Could also be done with OEAddExplicitHydrogens but 
            # I'm not sure if openeye can add residue info like rdkit. More testing maybe needed 
            rdmol = None
            nonstandard_atoms = []
            with Chem.SDMolSupplier('temp.sdf') as suppl:
                for mol in suppl:
                    rdmol = mol 

            stand_aminos = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", 
                            "SER", "THR", "VAL", "TRP", "TYR", "PYL", "SEC", "DA", "DC", "DG", "DT", "DI", "A", "C", "G", "U", "I", "HOH"]
            for atom_idx, atom_info in residue_info.items():
                res_name, res_num, atom_name, chain_id = atom_info
                if res_name not in stand_aminos:
                    nonstandard_atoms.append(atom_idx)
                atom = rdmol.GetAtomWithIdx(atom_idx)
                ri = Chem.AtomPDBResidueInfo()
                ri.SetResidueNumber(res_num)
                ri.SetChainId(chain_id)
                ri.SetResidueName(res_name)
                ri.SetName(atom_name)
                atom.SetPDBResidueInfo(ri)
            rdmol = Chem.AddHs(rdmol, addCoords=True, addResidueInfo=True, onlyOnAtoms=nonstandard_atoms)
            Chem.MolToPDBFile(rdmol, "pre_fixer_file.pdb")

            # run pdb fixer on the resulting file to close terminal groups and fill loops
            fixer = PDBFixer(filename='pre_fixer_file.pdb')
            fixer.findMissingResidues()
            # fixer.findNonstandardResidues()      # we often do contain Nonstandard residues
            # fixer.replaceNonstandardResidues()   # ^^
            # fixer.removeHeterogens(True)
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(7.0)
            # fixer.addSolvent(fixer.topology.getUnitCellDimensions())

            relative_file_path = Path(os.path.relpath(file, input_directory))
            output_path = str((output_directory / relative_file_path).parent / Path(f"{file.stem}.pdb"))
            PDBFile.writeFile(fixer.topology, fixer.positions, open(output_path, 'w'))
            log.append_info(f"{time.time()-start:2f}s\n")
        except Exception as e:
            log.append_ffile(str(file), f"Critical error: {e}")

log.printout()