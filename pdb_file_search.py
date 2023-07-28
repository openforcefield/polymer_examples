"""
Importing this folder allows for basic data searching functionality for pdb files only
"""
from pathlib import Path

class PDBFiles:

    @classmethod
    def search(self, pdb_file_name):
        self.cwdir = Path(__file__).parent.resolve()
        file_paths = []
        for file in (self.cwdir / Path("compatible_pdbs")).glob("**/*.pdb"):
            if "-s" in file.stem:
                polymer_name = file.stem.split("-s")[0]
            else:
                polymer_name = file.stem
            if file.is_file() and pdb_file_name == polymer_name:
                file_paths.append(file)

        return file_paths