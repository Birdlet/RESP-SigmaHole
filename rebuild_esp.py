from rdkit import Chem
import numpy as np

psi_bohr2angstroms = 0.52917720859 # Bohr to Angstroms conversion factor


if __name__ == "__main__":
    import os
    import argparse

    parser = argparse.ArgumentParser(description='resp')
    parser.add_argument("-l", help="ligand.mol2")
    parser.add_argument('-g', default='grid.dat', type=str, help="grid.dat")
    # parser.add_argument('-r', default='ANTECHAMBER.esp', type=str, help="reference ESP")
    parser.add_argument('-o', default='ffESP.mol2', type=str, help="default ouput ESP.mol2 name")


    args = parser.parse_args()

    atoms = []
    mol=Chem.MolFromMol2File(args.l, removeHs=False)
    conf =  mol.GetConformer()
    for i, at in enumerate(mol.GetAtoms()):
        aname = at.GetProp('_TriposAtomName')
        atype = at.GetProp('_TriposAtomType')
        achg = at.GetProp('_TriposPartialCharge')
        p = conf.GetAtomPosition(i)
        atoms.append( (np.array((p.x, p.y, p.z)), float(achg)))

    points = np.loadtxt(args.g)

    ofname = os.path.split(args.o)[-1]
    mol2s = "@<TRIPOS>MOLECULE\n%s\n %d 0 1\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n" % (ofname, points.shape[0])
    for point in points:
        potential = 0
        for atom, q in atoms:
            d = (point - atom)/psi_bohr2angstroms
            d2 = np.sum(d**2)
            potential += q/d2

        # "  %5d H          -4.1596   -0.0665   -0.0000 N.2   1    UNL        -0.5079"
        mol2_line = "   %4d %-8s  %8.4f%10.4f%10.4f %-5s %-4s %-8s  %8.4f\n" % (i+1, "H", 
                        point[0], point[1], point[2], "H", "1", "ESP", potential)
        mol2s += mol2_line
    
    
    with open(args.o, "w") as f:
        f.write(mol2s)
