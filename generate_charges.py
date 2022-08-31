from pathlib import Path
import sys
from openff.toolkit.utils import toolkit_registry
from rdkit import Chem
from openff.toolkit.topology.molecule import FrozenMolecule, Molecule
from openff.units import unit
from openff.toolkit.utils.toolkits import RDKitToolkitWrapper, OpenEyeToolkitWrapper
from rdkit import Chem
from openff.toolkit.topology import Topology

def generate_charged_molecule(pdbfile, substructure_file):
    # mol should already have one conformer...

    off_topology, _, error = Topology.from_pdb_and_monomer_info(pdbfile, substructure_file, strict=True)
    # here, we assume that the topology only has ONE simple homopolymer. Later, all molecules 
    # can be extracted and charged
    mol = next(off_topology.molecules) # get the first molecule
    # get some conformers to run elf10 charge method. By default, `mol.assign_partial_charges`
    # uses 500 conformers, but we can generate and use 10 here for demonstration
    mol.generate_conformers(
                    n_conformers=10,
                    rms_cutoff=0.25 * unit.angstrom,
                    make_carboxylic_acids_cis=True,
                    toolkit_registry=OpenEyeToolkitWrapper()
                ) # very slow for large polymers! 
    # finally, assign partial charges using those 10 conformers generated 
    mol.assign_partial_charges(
                    partial_charge_method='am1bccelf10', 
                    use_conformers=mol.conformers,
                    toolkit_registry=OpenEyeToolkitWrapper()
                )
    # code for exact how the above function works can be found in openff/toolkit/utils/openeye_wrapper.py under the assign_partial_charges() function
    return mol

if __name__ == "__main__":
    name = "naturalrubber" # finds this file name in the "compatible_pbds" folder
    pdb_file = None
    json_file = None
    for file in Path(Path.cwd() / Path('polymer_examples/compatible_pdbs')).glob("**/*.pdb"):
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

    charged_mol = generate_charged_molecule(str(pdb_file), json_file)
    print(f"final molecular charges: {charged_mol.partial_charges}")

    # note: the charged_mol has metadata about which monomers were assigned where
    # as a result of the chemicaly info assignment. This can be a way to break
    # up the molecule into repeating sections to partition the library charges 
    for atom in charged_mol.atoms:
        assert(atom.metadata['already_matched'] == True)
        print(atom.metadata['residue_name'])