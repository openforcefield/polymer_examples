from rdkit import Chem, Geometry
from rdkit.Chem import AllChem
from pathlib import Path
import os
import numpy as np
import math
from copy import deepcopy
import random
from collections import defaultdict

random.seed(6626)

current_dir = Path(__file__).parent.resolve()
os.chdir(current_dir)
os.chdir("..")

template_dir = Path("uncleaned_pdbs/homopolymer_templates")
for file in template_dir.glob("**/*.pdb"):
    # if file.stem != "peg_modified":
    #     continue
    print(file.stem)
    rdmol = Chem.MolFromPDBFile(str(file), removeHs=False, proximityBonding=False)
    if rdmol == None:
        continue
    pts = rdmol.GetConformer().GetPositions()
    n_pts, _ = pts.shape

    com = np.average(pts, axis=0) # center of mass

    # # transpose points to lie at origin 
    # pts = pts - com
    
    P = np.eye(n_pts) * (1-(1/n_pts)) + ((-(np.eye(n_pts)-1))*(-1/n_pts))
    cov = np.matmul(np.matmul(np.transpose(pts),P),pts)
    eval, evec = np.linalg.eig(cov)

    # get max eval
    princ_id = np.argmax(eval)
    princ_vec = evec[:,princ_id] # direction of poitns
    minor_ids = [i for i in range(3) if i != princ_id]

    evec2 = np.transpose(evec[:, minor_ids[0]])
    eval2 = eval[minor_ids[0]]

    evec3 = np.transpose(evec[:, minor_ids[1]])
    eval3 = eval[minor_ids[1]]

    # define a cylinder to stack molecules together
    r_max = 0
    for i in range(n_pts):
        AP = com - pts[i,:]
        r = np.linalg.norm(np.cross(AP, princ_vec)) / np.linalg.norm(princ_vec)
        if r > r_max:
            r_max = r

    # create new conformers based on the two minor eigenvalues
    
    n_new_polymers = random.randint(3,9)**2 # generate a random square between 9 and 100 to use as the size of our final polymers
    n_x = int(math.sqrt(n_new_polymers))
    n_y = int(math.sqrt(n_new_polymers))
    X, Y = np.meshgrid(np.linspace(-n_x/2.0, n_x/2.0, n_x), np.linspace(-n_y/2.0, n_y/2.0, n_y))

    new_conformers = np.array([])
    
    erdmol = Chem.Mol()
    for i in range(0, n_new_polymers): # keep the first conformer 
        x = X.flatten()[i] * 3 * r_max
        y = Y.flatten()[i] * 3 * r_max
        new_pts = pts + x*evec2 + y*evec3

        # final translation to help everthing fit in normal box sizes
        new_pts = new_pts - com

        new_rdmol = deepcopy(rdmol)
        new_rdmol.RemoveAllConformers()

        conf = Chem.Conformer(n_pts)
        for j in range(n_pts):
            conf.SetAtomPosition(j, Geometry.Point3D(new_pts[j][0], new_pts[j][1], new_pts[j][2]))
        new_rdmol.AddConformer(conf, assignId=True)
        erdmol = Chem.CombineMols(erdmol, new_rdmol)
        # new_pts = np.reshape(new_pts, (n_pts, 3, 1))
        # if new_conformers.size == 0:
        #     new_conformers = new_pts
        # else:
        #     new_conformers = np.concatenate([new_conformers, new_pts], axis=2)

    # make atom and residue names unique since are are pasting these files straight into "compatible_pdbs/"
    element_counts = defaultdict(int)
    res_number = 1
    for atom in erdmol.GetAtoms():
        ri = atom.GetPDBResidueInfo()
        name = ri.GetName()
        if len(atom.GetSymbol()) == 1:
            new_name = " " + atom.GetSymbol() + f"{element_counts[atom.GetSymbol()]:02d}"
        else:
            new_name = " " + atom.GetSymbol() + f"{element_counts[atom.GetSymbol()]:01d}"
        ri.SetResidueNumber(res_number)
        ri.SetChainId("1")
        ri.SetResidueName("UNK")
        ri.SetName(new_name)

        element_counts[atom.GetSymbol()] += 1
        if element_counts[atom.GetSymbol()] >= 100 / (10**(len(atom.GetSymbol())-1)): # if any atoms exceed 100 
            res_number += 1
            element_counts = defaultdict(int)
    Chem.MolToPDBFile(erdmol, f"compatible_pdbs/simple_polymers/{file.stem}-s{n_new_polymers}.pdb")