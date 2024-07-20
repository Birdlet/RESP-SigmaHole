#https://github.com/rdkit/rdkit/discussions/6710

from rdkit import Chem
from rdkit.Chem import rdMolTransforms as Trans

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='rdkit2pdbqt')
    parser.add_argument("-l", help="ligand.mol2")
    # parser.add_argument('-s', type=str, help="SMARTS to define substructure", required=False)
    parser.add_argument('-b', type=str, help="dihedral difined 4,5,6,7")
    parser.add_argument('-o', type=str, help="output")
    
    args = parser.parse_args()
    #parsing the command

    mol = Chem.MolFromMol2File(args.l, removeHs=False)
    conf = mol.GetConformers()[0]


    bond_idx = [int(b)-1 for b in args.b.split(",")]
    if args.o.endswith("pdb"):
        sdf = ""
        for deg in range(0, 360, 10):
            Trans.SetDihedralDeg(conf, bond_idx[0], bond_idx[1], bond_idx[2], bond_idx[3], deg)
            sdf += "CRYST1  200.000  200.000  200.000  90.00  90.00  90.00 P 1           1" + Chem.MolToPDBBlock(mol).replace("END","ENDMDL")
        
        with open(args.o, "w") as f:
            f.write(sdf)
    else:
        sdf = ""
        for deg in range(0, 360, 10):
            Trans.SetDihedralDeg(conf, bond_idx[0], bond_idx[1], bond_idx[2], bond_idx[3], deg)
            sdf += Chem.MolToMolBlock(mol)+"$$$$\n"
        
        with open(args.o, "w") as f:
            f.write(sdf)
    
