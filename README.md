# polymer_examples
A repository for testing the PDB loading capabilities for the [Open Force Field toolkit](https://github.com/openforcefield/openff-toolkit). To run the monomer template generation and testing scripts, you will need to have RDKit and OpenMM installed.
**Note**: The testing scripts included are stable for version 0.13.2 of the OpenFF toolkit. 
## Developer Installation with Conda
To install and test the toolkit, first download and install the Anaconda Distribution from [here](https://www.anaconda.com/download). Then, run:
```sh
conda create -n openff-toolkit -c conda-forge openff-toolkit
conda activate openff-toolkit
```
If you would prefer to use a developer install to suggest changes to the toolkit or fix bugs, follow the below steps from the [OpenFF Documentation](https://docs.openforcefield.org/projects/toolkit/en/latest/users/developing.html#setting-up-a-development-environment).
```sh
git clone https://github.com/openforcefield/openff-toolkit
cd openff-toolkit/
# Create a conda environment with dependencies from env/YAML file
conda env create -n openff-dev -f devtools/conda-envs/test_env.yaml
conda activate openff-dev
# Perform editable/dirty dev install
pip install -e .
```
You may also wish to manually download the [openff-interchange](https://github.com/openforcefield/openff-interchange) and [openff-forcefields](https://github.com/openforcefield/openff-forcefields) repos, which can be installed with a similar procedure (`pip install -e .` inside the activated openff-dev environment) for a closer look at the toolkit parameterization process.
## Using the toolkit to load PDBs
Begin with the `minimum_reproducing_example.py` file for a self-contained demonstration of the toolkit:
```sh
conda activate openff-toolkit
python minimum_reproducing_example.py
```
If the parameterization was successful, you should be met with a successful message and final system energy. This example, like the parameterization done in the rest of this repo, uses the sage 2.0.0 forcefield (`general_offxml`), ported amberff14 forcefield (`amber_offxml`), and the tip3b water model (`water_model`). The simulation is performed in OpenMM 1000 steps following energy minimization. Note the location of the pdb files which can be found compatible_pdbs/proteins/6cww.pdb and the location of the monomer templates found at monomer_generation/json_files/6cww.json. These additional monomers are loaded with the default amino acid templates which come pre-installed with the tooklit and can be found [here](https://github.com/openforcefield/openff-toolkit/blob/main/openff/toolkit/data/proteins/aa_residues_substructures_explicit_bond_orders_with_caps_explicit_connectivity.json) (otherwise, we would have lots of templates to specify!). 
## Generating Monomer Templates and Testing PDBs
Now with a single example out of the way, we can go ahead with the rest of the files in the repository!

1. First, run the pdb cleaner (this is unecessary as of 7/31/2023 since the post-cleaned files are uploaded. Running this step requires the pdbfixer library and OpenEye (which is a licenced toolkit)):
```sh
python pdb_cleaner.py
```
This step takes the various files in \uncleaned_pdbs and transforms them into well-formated files in \compatible_pdbs. 
For peptoids, simple_polymers, and sugars, the files are simply loaded with OpenEye or RDKit, and the residue names are modified to make them unique (OpenMM requres unique atom names on unique residue numbers, which is important because Topology.from_pdb() requires OpenMM!).
For proteins, the CIF files are loaded with openeye to remove alternate locations and populate chemical information. We then output and re-load into rdkit to add hydrogens to non-standard residues. Finally, we load again into the PDBFixer from OpenMM to add hydrogens to amino acids and close loops (MD simulations require closed loops in the OpenFF Toolkit). 
For DNA and RNA, OpenEye is used to load CIF files. A problem with chemistry reading required the use of a manual substructure search for the P-bridge (`P(~[OD1])(~[OD1])(~[OD2]-*)(~[OD2]-*)`) followed by manual insertion of the correct formal charge and bond order around those substructures. 

2. Run the monomer generation and testing. [Note: this will take around 4 hours on a basic system such as a peronsal computer. If you would prefer to test files individually, see the minimum reproducing example in the root directory] Here, we input monomer smarts from monomer_generation/monomer_smiles_input.py, create the formated monomer templates with code in monomer_generation/substructure_generator.py, and output to monomer_generation/json_files. After the monomer templates are stored, we perform a trial load with `Topology.from_pdb()` followed by an energy minimization. The times to perform the loading and parameterization steps are stored for later.
```sh
python monomer_generation/generate_monomers.py
```

3. Plot the final runtimes.
```sh
python monomer_generation/plot_loading_times.py
```
## Other Notes
The code used to create the combined monomer library that can load all files except for polyphenylene monomers (conflicts with conjugated aromatic systems) and 7wcc (Fe+2 ion conflicts with Fe+3 in earlier file) can be found at monomer_generation/create_and_validate_library.py

The code used to validate the new toolkit version with the old is in \validation_with_old_tk. Here, we loaded 9 pdb files that did not require custom substructures (so only default amino acids) and parameterized them using the 0.13.1 version and 0.13.2 (post _custom_substructure comit) version of the toolkit. Comparing the parameters of each molecule, which did not require anything but an exact file comparison, showed that the tooklit output was exactly the same for large biopolymers. 
