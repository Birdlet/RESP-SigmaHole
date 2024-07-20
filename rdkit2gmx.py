#!/usr/bin/env python3
'''
(c)Copyright 202024, Ding Kang
2024.02.02
'''
from __future__ import absolute_import, print_function
from math import isnan, isinf
from itertools import combinations

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry.rdGeometry import Point3D

import sys
from utils import _sybyl_atom_type
nm2A = 10


# [ virtual_sites2 ]
# ; Site  from        funct  d
# dummyID  p1ID  p0ID  2  distance   
virtual_sites2 = "%4d  %4d  %4d     2      %.3f"


def MolFromITPFile(itp):
    top_blocks = [
                 "; $NAME_LP.top created by rdkit2gmx.py\n" + \
                  "[ defaults ]\n" + \
                  "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n" + \
                  "1               2               yes             0.5     0.8333333333\n",
    ]

    atom_lines = []

    flag = False
    name = ""
    top_block = ""

    with open(itp, 'r') as f:
        lines = list(f.readlines())
        for i, line in enumerate(lines):
            if line.strip().startswith("[ atomtypes ]"):
                top_block += line
                top_block += "    DLP    DLP     0.0000      0.050     A    0.00000000000e+00    0.000000e+00 \n"
                continue
            if line.strip().startswith("[ atoms ]"):
                flag = True
                top_block += line
                continue
            if line.strip().startswith("[ bonds ]"):
                flag = False
                top_block += "$ATOM\n\n$VSITE\n\n"
                top_block += line
                continue

            if line.strip().startswith("[ moleculetype ]"):
                j = 1
                while not name or name.startswith(";"):
                    name = lines[i+j].strip().split()[0]
                    j += 1

            if flag:
                if not line.strip(): continue
                if line.strip().startswith(";"): continue
                atom_lines.append(line.split(";")[0].strip())
            else:
                top_block += line

    atoms = []
    for atom in atom_lines:
        cols = atom.strip().split()
        atoms.append(cols)

    top_blocks.append(top_block)
    top_blocks.append( "\n\n[ system ]\n" + \
                 " %s\n\n" % name + \
                "[ molecules ]\n" + \
                "; Compound        nmols\n" + \
                " %s              1\n" % name)

    return atoms, top_blocks


def GMXAtomLines(mol):
    """Create a list with GMX atom lines for each atom in molecule.
    """
    gro_line = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
    gro_line = "    1%-5s%5s%5d%24s"
    
    atom_lines = [line.replace('HETATM', 'ATOM  ')
                  for line in Chem.MolToPDBBlock(mol).split('\n')
                  if line.startswith('HETATM') or line.startswith('ATOM')]
    
    gro_lines = []
    for idx, atom in enumerate(mol.GetAtoms()):
        coord = float(atom_lines[idx][30:38])/nm2A, float(atom_lines[idx][38:46])/nm2A, float(atom_lines[idx][46:54])/nm2A
        coords = "%8.4f%8.4f%8.4f" % coord
        resname = atom_lines[idx][17:20].strip()
        atname = atom_lines[idx][12:16].strip()
        atid = atom_lines[idx][6:11].strip()

        gro_lines.append( gro_line % (resname, atname, int(atid), coords) )
    return gro_lines


def Mol2AtomLines(mol, dummy_atoms, charges, resname=""):
    atoms = []
    mol2_line = "   %4d %-8s  %8.4f%10.4f%10.4f %-5s %-4s %-8s  %8.4f"
    assert mol.GetNumAtoms() == len(charges)
    conf =  mol.GetConformer()
    for idx, at in enumerate(mol.GetAtoms()):
        try:
            aname = at.GetProp("_TriposAtomName")
            atype = at.GetProp("_TriposAtomType")
        except:
            aname, atype = _sybyl_atom_type(at)
        p = conf.GetAtomPosition(idx)
        atoms.append( mol2_line % (idx+1, aname, p.x, p.y, p.z, atype, "1", resname, float(charges[idx])) )

    return atoms + dummy_atoms

def MolToGMXBlock(mol, itp=None):
    """Write RDKit Molecule to a GMX block
    Parameters
    ----------
        mol: rdkit.Chem.rdchem.Mol
            Molecule with a protein ligand complex
    Returns
    -------
        block: str
            String wit GMX encoded molecule
    """
    gro_line = "    1%-5s%5s%5d%8.4f%8.4f%8.4f"
    #      6 C6         44.3900   34.4900   23.6500 C.ar  100  AAA100      0.0000
    mol2_line = "   %4d %-8s  %8.4f%10.4f%10.4f %-5s %-4s %-8s  %8.4f"

    # make a copy of molecule
    mol = Chem.Mol(mol)

    atom_lines = GMXAtomLines(mol)
    assert len(atom_lines) == mol.GetNumAtoms()
    resname = atom_lines[0][5:10].strip()

    dummy_lines = []
    mol2_dummy = []
    gmx_virtual_sites = []
    dummy_idx = mol.GetNumAtoms()+1

    if itp:
        itp_atoms, top_blocks = MolFromITPFile(itp)
        assert len(atom_lines) == len(itp_atoms)

    d_scale = 0.8
    e_scale = 1.0 # extent distance
    dummy_pairs = []
    for idx in range(len(atom_lines)):
        # add dummy / virtual site
        _atom = mol.GetAtomWithIdx(idx)
        if _atom.GetAtomicNum() in (17,35,53):
            _nb = _atom.GetNeighbors()
            if len(_nb) == 1:
                _nb = _nb[0]
                p1 = mol.GetConformer().GetAtomPosition(_nb.GetIdx())
                p0 = mol.GetConformer().GetAtomPosition(_atom.GetIdx())
                # p2 = Point3D(0.6,0.6,0.6)*(p0-p1)+p0

                if _atom.GetAtomicNum() == 35:
                    d_scale = 1.89 # 1.37
                elif _atom.GetAtomicNum() == 53:
                    d_scale = 2.20 # 1.73
                else: # _atom.GetAtomicNum() == 17:
                    d_scale = 1.64
                d_scale *= e_scale
                v = p0-p1
                v.Normalize()
                p2 = v*d_scale + p0
                d = (p2-p0).Length()/nm2A # 10 for 1 nm = 10 A
                d = (d_scale+(p2-p0).Length())/nm2A
                dummy_line = gro_line % (resname, "DLP"+str(dummy_idx-mol.GetNumAtoms()), dummy_idx, p2.x/nm2A, p2.y/nm2A, p2.z/nm2A)
                mol2_dummy.append( mol2_line % (dummy_idx, "DLP"+str(dummy_idx-mol.GetNumAtoms()), p2.x, p2.y, p2.z, "Du", "1", resname, 0.05) )
                gmx_virtual_sites.append( virtual_sites2 % (dummy_idx, _nb.GetIdx()+1, _atom.GetIdx()+1, d) )
                
                if itp:
                    itp_atoms[_nb.GetIdx()][6] = "%.6f" % (float(itp_atoms[_nb.GetIdx()][6]) - 0.05)
                    itp_atoms.append( [str(dummy_idx), "DLP", itp_atoms[0][2], itp_atoms[0][3], "DLP"+str(dummy_idx-mol.GetNumAtoms()),
                                    str(dummy_idx), "0.050000", "0.0000"])

                dummy_pairs.append((dummy_idx,))
                dummy_idx += 1
                dummy_lines.append(dummy_line)
        elif _atom.GetAtomicNum() == 16:
            if _atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2 and len(_atom.GetNeighbors())==2:
                apl1 = _atom.GetNeighbors()[0]
                apr1 = _atom.GetNeighbors()[1]
                d_scale *= e_scale
                #if apl1.GetAtomicNum() == 6 and apr1.GetAtomicNum() == 6:
                if _atom.GetIsAromatic():
                    p0 = mol.GetConformer().GetAtomPosition(_atom.GetIdx())
                    pl1 = mol.GetConformer().GetAtomPosition(apl1.GetIdx())
                    pr1 = mol.GetConformer().GetAtomPosition(apr1.GetIdx())
                    pl2 = (p0-pl1)*d_scale + p0
                    pr2 = (p0-pr1)*d_scale + p0
                    dl = (p0-pl1).Length()/nm2A * (1+d_scale)
                    dr = (p0-pr1).Length()/nm2A * (1+d_scale)
                    
                    dummy_line = gro_line % (resname, "DLP"+str(dummy_idx-mol.GetNumAtoms()), dummy_idx, pl2.x/nm2A, pl2.y/nm2A, pl2.z/nm2A)
                    gmx_virtual_sites.append( virtual_sites2 % (dummy_idx, _atom.GetNeighbors()[0].GetIdx()+1, _atom.GetIdx()+1, dl) )
                    
                    _nb = _atom.GetNeighbors()[0]
                    mol2_dummy.append( mol2_line % (dummy_idx, "DLP"+str(dummy_idx-mol.GetNumAtoms()), pl2.x, pl2.y, pl2.z, "Du", "1", resname, 0.05) )
                    if itp:
                        itp_atoms[_nb.GetIdx()][6] = "%.6f" % (float(itp_atoms[_nb.GetIdx()][6]) - 0.05)
                        itp_atoms.append( [str(dummy_idx), "DLP", itp_atoms[0][2], itp_atoms[0][3], "DLP"+str(dummy_idx-mol.GetNumAtoms()),
                                    str(dummy_idx), "0.050000", "0.0000"])

                    dummy_idx += 1

                    dummy_lines.append(dummy_line)
                    dummy_line = gro_line % (resname, "DLP"+str(dummy_idx-mol.GetNumAtoms()), dummy_idx, pr2.x/nm2A, pr2.y/nm2A, pr2.z/nm2A)
                    gmx_virtual_sites.append( virtual_sites2 % (dummy_idx, _atom.GetNeighbors()[1].GetIdx()+1, _atom.GetIdx()+1, dr) )
                    
                    _nb = _atom.GetNeighbors()[1]                    
                    mol2_dummy.append( mol2_line % (dummy_idx, "DLP"+str(dummy_idx-mol.GetNumAtoms()), pr2.x, pr2.y, pr2.z, "Du", "1", resname, 0.05) )
                    if itp:
                        itp_atoms[_nb.GetIdx()][6] = "%.6f" % (float(itp_atoms[_nb.GetIdx()][6]) - 0.05)
                        itp_atoms.append( [str(dummy_idx), "DLP", itp_atoms[0][2], itp_atoms[0][3], "DLP"+str(dummy_idx-mol.GetNumAtoms()),
                                    str(dummy_idx), "0.050000", "0.0000"])

                    dummy_pairs.append((dummy_idx-1, dummy_idx))
                    dummy_idx += 1
                    dummy_lines.append(dummy_line)


    # print("\n".join(gmx_virtual_sites))
    if itp:
        itp_atom_format = "%6s%11s%7s%9s%7s%7s%13s%11s"
        if 1+mol.GetNumAtoms() == dummy_idx:
            charges = [ at[6] for at in itp_atoms]
        else:
            charges = [ at[6] for at in itp_atoms][:1+mol.GetNumAtoms() - dummy_idx]
        if 1+mol.GetNumAtoms() < dummy_idx:
            vsite_block = "[ virtual_sites2 ]\n; Site  from        funct  d\n" + "\n".join(gmx_virtual_sites)
        else:
            vsite_block = ""
        itp_atoms = [ itp_atom_format % tuple(at) for at in itp_atoms]
        atom_block = "\n".join(itp_atoms)
        mol2_atoms = Mol2AtomLines(mol, mol2_dummy, charges, resname)
        top_blocks[1] = top_blocks[1].replace("$ATOM", atom_block).replace("$VSITE", vsite_block)
    else:
        itp_atoms = []
        charges = [0.0 for _ in range(mol.GetNumAtoms())]
        mol2_atoms = ""
        top_blocks = ["", ""]
        atom_block = ""


    
    gro_block = '\n'.join([resname, "%d"%(dummy_idx-1)] + atom_lines + dummy_lines + ["  10.00000  10.00000  10.00000\n",])

    # print("-"*80)
    # print("[ atoms ]")
    # print(atom_block)
    # print("-"*80)
    # print("[ virtual_sites2 ]")
    # print(vsite_block)
    # print("-"*80)
    # print(gro_block)
    # print("-"*80)

    return gro_block, top_blocks[1], "".join(top_blocks), mol2_atoms, dummy_pairs


if __name__ == "__main__":
    import os
    import argparse
    import json
    parser = argparse.ArgumentParser(description='rdkit2pdbqt')
    parser.add_argument("-l", help="ligand.mol2")
    parser.add_argument('-t', type=str, help="ligand.itp")
    parser.add_argument('-o', type=str, help="default ouput base name")
    parser.add_argument('--mol2_only', action="store_true")
    # parser.add_argument('-ot', type=str, help="Topology output")
    # parser.add_argument('-oi', type=str, help="ITP output")
    # parser.add_argument('-os', type=str, help="Coord output")

    args = parser.parse_args()
    #parsing the command

    mol=Chem.MolFromMol2File(args.l, removeHs=False)
    if not mol:
        os.system("obabel -imol2 %s -omol -O %s" % (args.l, args.l[:-1]))
        mol=Chem.MolFromMolFile(args.l[:-1], removeHs=False)
    if not mol:
        os.system("obabel -imol2 %s -oxyz -O %s" % (args.l, args.l[:-4]+"xyz"))
        os.system("obabel -ixyz %s -omol -O %s" % (args.l[:-4]+"xyz", args.l[:-1]))
        mol=Chem.MolFromMolFile(args.l[:-1], removeHs=False)
    itp = args.t
    
    gro_block, itp_block, top_block, mol2_atoms, dummy_pairs = MolToGMXBlock(mol, itp)

    with open(args.l) as f:
        newline = ""
        lines = list(f.readlines())
        i = 0 
        while i < len(lines):
            if lines[i].startswith("@<TRIPOS>MOLECULE"):
                newline += lines[i]
                newline += lines[i+1]
                natoms = lines[i+2]
                newline += " " + " ".join([str(len(mol2_atoms)) ,] + natoms.strip().split()[1:]) + "\n"
                i += 3
                continue
            if lines[i].startswith("@<TRIPOS>ATOM"):
                newline += lines[i]
                newline += "\n".join(mol2_atoms) + "\n"
                while not lines[i].startswith("@<TRIPOS>BOND"):
                    i += 1
                newline += lines[i]
                i += 1
                continue
            newline += lines[i]
            i+=1
    
    omol2 = args.o + "_LP.mol2"
    otop = args.o + "_LP.top"
    oitp = args.o + "_LP.itp"
    ogro = args.o + "_LP.gro"
    odummy = args.o + "_LP.pairs"

    if args.mol2_only:
        with open(omol2, "w") as f:
            f.write(newline)
    else:
        with open(otop, "w") as f:
            f.write(top_block)
        with open(oitp, "w") as f:
            f.write(itp_block)
        with open(ogro, "w") as f:
            f.write(gro_block)
        with open(odummy, "w") as f:
            json.dump(dummy_pairs, f)
