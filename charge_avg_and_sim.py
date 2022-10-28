from chavg_tools import *
from IPython.display import clear_output
import os

def main():
    mol_name = 'polyvinylchloride'

    offxml_src = Path('xml examples/openff_unconstrained_with_library_charges-2.0.0.offxml')
    output_folder = Path(f'averaged_polymers/{mol_name}')
    output_folder.mkdir(exist_ok=True)
    lc_path = output_folder/f'new {mol_name} charges.offxml' # path to output library charges to

    mol, topology = fetch_mol(mol_name)  # will raise exception if files for molecule are not found
    if lc_path.exists(): # check if library charges have already been generated for this molecule and load them...
        forcefield = ForceField(lc_path, allow_cosmetic_attributes=True)
    else: # ...otherwise if they haven't, generate them
        cmol = generate_molecule_charges(mol, toolkit_method='openeye', partial_charge_method='am1bcc') # perform AM1BCC
        clear_output() # for Jupyter notebooks only, can freely comment this out
        avgs = averaged_charges_by_SMARTS(cmol) # average charges over unique residues - placed after clear so we can see what averages are computed
        for smiles, subdict in avgs.items():
            print(f'{smiles}\n\t{subdict}\n')

        forcefield = write_new_library_charges(avgs, offxml_src, output_path=lc_path)

    interchange = Interchange.from_smirnoff(force_field=forcefield, topology=topology) # generate Interchange with new library charges prior to writing to file
    sim = create_sim_from_interchange(interchange)
    run_simulation(sim, output_folder=output_folder, output_name=mol_name, num_steps=10000, record_freq=10)

if __name__ == '__main__':
    main()