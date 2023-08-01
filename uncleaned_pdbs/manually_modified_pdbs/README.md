The following files were either manually created or manually modified due to bad formatting. 

The paam_modified, peg_modified, and pnipam_modified files were manually generated and provided by the Coray Colina group. 

The rest had incorrect bonds or atoms. For example, often times when protein files are generated from X-Ray diffraction (which most of these files were) the original CIF file in combination with the OpenEye toolkit which reads them interprets extra ligands as just being "tacked" on to the protein without an actual bond that connects the ligand/small molecule to the protein. These files/proteins represent those cases were bonds were manually moved and created so that the pdb CONECT record was "correct".
