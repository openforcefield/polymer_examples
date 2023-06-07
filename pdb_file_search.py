"""
Importing this folder allows for basic data searching functionality for pdb files only
"""
from pathlib import Path

class PDBFiles:

    @classmethod
    def search(self, pdb_file_name):
        self.cwdir = Path(__file__).parent.resolve()
        file_path = None
        for file in (self.cwdir / Path("compatible_pdbs")).glob("**/*.pdb"):
            if file.is_file() and file.stem == pdb_file_name:
                file_path = file

        return file_path