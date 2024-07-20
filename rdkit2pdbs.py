#!/usr/bin/env python3
'''
(c)Copyright 2015, Maciej Wojcikowski Revision caf5d84a.
Wojcikowski, M., Zielenkiewicz, P. & Siedlecki, P. Open Drug Discovery Toolkit (ODDT): a new open-source player in the drug discovery field. J Cheminform 7, 26 (2015). https://doi.org/10.1186/s13321-015-0078-2
Modified by XU Ximing
xuximing@ouc.edu.cn
2022.12.10
'''
from __future__ import absolute_import, print_function
from math import isnan, isinf
from itertools import combinations

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry.rdGeometry import  Point3D

from rdkit2gmx import *

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='rdkit2pdbqt')
    parser.add_argument("-l", help="ligand.sdf")
    parser.add_argument('-o', type=str, help="output")

    args = parser.parse_args()
    #parsing the command

    sd = Chem.SDMolSupplier(args.l, removeHs=False)
    with open(args.o, "w") as f:
        for mol in sd:
            gro_block, itp_block, top_block, mol2_atoms, dummy_pairs  = MolToGMXBlock(mol, None)
            f.write(gro_block)

