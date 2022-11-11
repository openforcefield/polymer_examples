from chavg_tools import *
from IPython.display import clear_output


def main():
    # Perform charge averaging on all target molecules which don't already have averaged LCs; 
    # Load forcefield for those which already do 
    sample_mols = ['PEO_PLGA']
    run_sims = True

    offxml_src = Path('xml examples/openff_unconstrained_with_library_charges-2.0.0.offxml')
    polymer_folder = Path('compatible_pdbs/simple_polymers')

    for mol_name in sample_mols: #mols_to_use:
        # DEFINING PATHS, CREATING FOLDERS, AND FETCHING FILES
        print(mol_name)
        pdb_path = Path(f'compatible_pdbs/simple_polymers/{mol_name}.pdb')
        charged_json = Path(f'charged_jsons/{mol_name}_with_charges.json')
        default_json = polymer_folder/f'{mol_name}.json'
        json_path = charged_json if charged_json.exists() else default_json
        
        output_folder = Path(f'averaged_polymers/{mol_name}')
        output_folder.mkdir(exist_ok=True)
        lc_path = output_folder/f'new {mol_name} charges.offxml' # path to output library charges to

        # LOAD MOLECULE AND TOPOLOGY, ATTEMPT TO APPLY LIBRARY CHARGES
        mol, topology = load_mol_and_topo(pdb_path, json_path)  # will raise exception if files for molecule are not found
        
        if lc_path.exists(): # check if library charges have already been generated for this molecule
            forcefield = ForceField(lc_path, allow_cosmetic_attributes=True)
        else:
            cmol = generate_molecule_charges(mol, toolkit_method='openeye', partial_charge_method='am1bcc') # perform AM1BCC
            clear_output() # for Jupyter notebooks only, can freely comment this out
            avgs = get_averaged_charges(cmol) # average charges over unique residues - placed after clear so we can see what averages are computed
            for averaged_res in avgs:
                print(averaged_res, '\n')

            forcefield, lib_chgs = write_new_library_charges(avgs, offxml_src, output_path=lc_path)
            
            # CREATE JSON WITH AVERAGED CHARGES IF ONE DOES NOT ALREADY EXIST
            if not charged_json.exists():
                with default_json.open('r') as old_json:
                    json_dat = json.load(old_json)

                charge_entry = {avgd_res.residue_name : avgd_res.charges for avgd_res in avgs}
                json_dat['charges'] = charge_entry

                charged_json.touch()
                with charged_json.open('w') as new_json:
                    json.dump(json_dat, new_json, indent=4)

        # RUN OpenMM SIMULATION FOR TARGET MOLECULE
        if run_sims:
            forcefield = ForceField(lc_path, allow_cosmetic_attributes=True)
            interchange = Interchange.from_smirnoff(force_field=forcefield, topology=topology)#, charge_from_molecules=[cmol]) # generate Interchange with new library charges prior to writing to file
            sim = create_sim_from_interchange(interchange)
            run_simulation(sim, output_folder=output_folder, output_name=mol_name, num_steps=1000, record_freq=10)

if __name__ == '__main__':
    main()