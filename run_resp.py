#/!bin/bash

import os

psi_bohr2angstroms = 0.52917720859 # Bohr to Angstroms conversion factor
nm2A = 10

gas = """
cd $_pwd
mkdir gas;cd gas
antechamber -i ../gas.gesp -fi gesp -o gas.mol2 -fo mol2 -at sybyl -c resp -nc -1
parmchk2 -i Lig.mol2 -f mol2 -o Lig.frcmod

tleap -f - <<_EOF
source leaprc.ff14SB
source leaprc.gaff2
loadamberparams Lig.frcmod
lig=loadmol2 Lig.mol2
check lig
saveamberparm lig Lig.prmtop Lig.inpcrd
quit
_EOF

acpype -p Lig.prmtop -x Lig.inpcrd

"""


def MolFromITPFile(itp):
    top_blocks = [
                 "; $NAME_LP.top created by rdkit2gmx.py\n" + \
                  "[ defaults ]\n" + \
                  "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n" + \
                  "1               2               yes             0.5     0.8333333333\n",
    ]

    atom_lines = []

    flag = False
    flag2 = False
    name = ""
    top_block = ""

    with open(itp, 'r') as f:
        lines = list(f.readlines())
        for i, line in enumerate(lines):
            if line.strip().startswith("[ atomtypes ]"):
                top_block += line
                continue
            if line.strip().startswith("[ atoms ]"):
                flag = True
                top_block += line
                continue
            if line.strip().startswith("[ bonds ]"):
                flag = False
                flag2 = False
                top_block += line
                continue
            if line.strip().startswith("[ virtual_sites2 ]"):
                flag2 = True
                top_block += "$ATOM\n\n"
                top_block += line
                continue
            if flag2:
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

def MolFromGroFile(gro):
    with open(gro) as f:
        lines = list(f.readlines())
        name = lines[0]
        natom = int(lines[1].strip())
        atoms = lines[2:2+natom]
    Mol = {
           "name": name,
           "natom": natom,
           "atoms": atoms
           }
    return Mol


def INDummy(infile, modinfile, ndummy, stage=1, dummy_pairs=None, equ=True):
    with open(infile) as f:
        lines = list(f.readlines())

    natom_line = 0
    for i,l in enumerate(lines):
        if l.startswith(" &end"):
            natom_line = i + 3
    cs = lines[natom_line].strip().split()
    lines[natom_line] = "%5s%5d\n" % (cs[0], int(cs[1])+ndummy)
    lines = lines[:natom_line + int(cs[1])+1]

    if stage == 1:
        if dummy_pairs:
            _ndummy = 0
            for d in dummy_pairs:
                if len(d) > 1:
                    lines.append("    0    0\n")
                    _ndummy += 1
                    for _ in d[1:]:
                        _ndummy += 1
                        lines.append("    0%5d\n" % d[0])
                else:
                    _ndummy += 1
                    lines.append("    0    0\n")
            assert _ndummy == ndummy
        else:
            for _ in range(ndummy):
                lines.append("    0    0\n")
    else:
        if dummy_pairs and (not equ):
            _ndummy = 0
            for d in dummy_pairs:
                if len(d) > 1:
                    lines.append("    0    0\n")
                    _ndummy += 1
                    for _ in d[1:]:
                        _ndummy += 1
                        lines.append("    0%5d\n" % 0)
                else:
                    _ndummy += 1
                    lines.append("    0  -99\n")
            assert _ndummy == ndummy
        else:
            for _ in range(ndummy):
                lines.append("    0  -99\n")

    lines.append("\n") 

    with open(modinfile, "w") as f:
        f.write("".join(lines))

def ESPDummy(espfile, modespfile, dummy_coords):
    with open(espfile) as f:
        lines = list(f.readlines())
    cs = lines[0][0:5], lines[0][5:10], lines[0][10:15]
    natoms = int(cs[0]) + len(dummy_coords)

    lines[0] = "%5d%5s%5s\n" % (natoms, cs[1], cs[2])
    dummys = []
    for coord in dummy_coords:
        _line = "%16s%16.7e%16.7e%16.7e\n" % ("", coord[0]/psi_bohr2angstroms*nm2A, 
                                              coord[1]/psi_bohr2angstroms*nm2A, 
                                              coord[2]/psi_bohr2angstroms*nm2A) 
        dummys.append( _line.replace("e", "E") )
    lines = lines[0:natoms+1-len(dummy_coords)] + dummys + lines[natoms+1-len(dummy_coords):]
    
    with open(modespfile, "w") as f:
        f.write("".join(lines))

def WriteQIN(qin, natoms, ndummy):
    with open(qin, "w") as f:
        for i in range(1, natoms+1):
            if i > natoms-ndummy:
                f.write("%10.6f" % 1.0)
            else:
                f.write("%10.6f" % 0.0)
            if i % 8 == 0:
                f.write("\n")


def MolFromMol2File(mol2):
    with open(mol2) as f:
        lines = list(f.readlines())
    at = 0
    atoms = []
    nlines = []
    for  i in range(len(lines)):
        if lines[i].lstrip().startswith("@<TRIPOS>ATOM"):
            at = 1
            nlines.append(lines[i])
            nlines.append("$ATOM")
            continue
        elif lines[i].lstrip().startswith("@<TRIPOS>BOND"):
            at = 0
        if at:
            atoms.append(lines[i].strip().split())
        else:
            nlines.append(lines[i])
    return atoms, "".join(nlines)


def ITPDummyReCharges(gro, itp, mol2, dummy_pairs=None):
    
    atoms, top_blocks = MolFromITPFile(itp)
    mol2_atoms, mol2_blocks = MolFromMol2File(mol2)
    ndummy = 0
    for atom in atoms:
        if atom[1] == "DLP":
            ndummy += 1
    Mol = MolFromGroFile(gro)
    dummy_coords = []
    #     "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
    for atom in Mol["atoms"]:
        if atom[10:15].strip()[:3] == "DLP":
            dummy_coords.append((float(atom[20:28]), float(atom[28:36]), float(atom[36:44])))

    # atc = "antechamber -i ../gas.gesp -fi gesp -o gas.mol2 -fo mol2 -at sybyl -c resp -nc %d" % nc
    # # > ANTECHAMBER_RESP1.IN ANTECHAMBER_RESP2.IN
    # os.system(atc)

    if not os.path.exists("MOD1"):
        os.mkdir("MOD1")
    if not os.path.exists("MOD2"):
        os.mkdir("MOD2")

    # RESP Stage 1
    WriteQIN("MOD2/qin", len(atoms), ndummy)
    ESPDummy("ANTECHAMBER.ESP", "ANTECHAMBER_MOD.ESP", dummy_coords)

    INDummy("ANTECHAMBER_RESP1.IN", "MOD1/ANTECHAMBER_RESP1_MOD.IN", ndummy, stage=1, dummy_pairs=dummy_pairs)
    resp1 = "cd MOD1; rm output; resp -i ANTECHAMBER_RESP1_MOD.IN -o output -p punch -q qin -t qout -e ../ANTECHAMBER_MOD.ESP -s esout"
    os.system(resp1)
    # RESP Stage 2
    INDummy("ANTECHAMBER_RESP2.IN", "MOD2/ANTECHAMBER_RESP2_MOD.IN", ndummy, stage=2, dummy_pairs=dummy_pairs)
    resp2 = "cd MOD2; rm output; resp -i ANTECHAMBER_RESP2_MOD.IN -o output -p punch -q ../MOD1/qout -t qout -e ../ANTECHAMBER_MOD.ESP -s esout"
    os.system(resp2)    
    
    with open("MOD2/qout") as f:
        charges = []
        for l in f.readlines():
            if l: charges.extend( [float(x) for x in l.strip().split()] )

    for i in range(len(atoms)):
        atoms[i][6] = "%.6f" % charges[i]

    itp_atom_format = "%6s%11s%7s%9s%7s%7s%13s%11s"
    mol2_line = "   %4s %-8s  %8s%10s%10s %-5s %-4s %-8s  %8.4f\n"

    itp_atoms = [ itp_atom_format % tuple(at) for at in atoms]
    _mol2_atoms = []
    for at, chg in zip(mol2_atoms, charges):
        _mol2_atoms.append( mol2_line % tuple(at[:-1]+[chg,]) )

    atom_block = "\n".join(itp_atoms)

    top_blocks[1] = top_blocks[1].replace("$ATOM", atom_block)
    itp_bocks = top_blocks[1]
    top_blocks = "".join(top_blocks)
    mol2_blocks = mol2_blocks.replace("$ATOM", "".join(_mol2_atoms))

    return itp_bocks, top_blocks, mol2_blocks


def Psi4RespCharge(lig, nc=0):
    gas = "antechamber -i ../gas.gesp -fi gesp -o gas.mol2 -fo mol2 -at sybyl -c resp -nc %(nc)d"
    acpype = "acpype -i gas.mol2 -n %(nc)d -c user"
    os.system(sol % nc)
    os.system(acpype % nc)


def Psi4Resp2Charge(lig, d=0.5, nc=0):
    sol = "antechamber -i ../sol.gesp -fi gesp -o sol.mol2 -fo mol2 -at sybyl -c resp -nc %(nc)d"
    gas = "antechamber -i ../gas.gesp -fi gesp -o gas.mol2 -fo mol2 -at sybyl -c resp -nc %(nc)d"

    mol2_line = "   %4d %-8s  %8.4f%10.4f%10.4f %-5s %-4s %-8s  %8.4f"

    
    for idx, at in enumerate(mol.GetAtoms()):
        aname = at.GetProp("_TriposAtomName")
        atype = at.GetProp("_TriposAtomType")
    
    with open("gas/qout") as f:
        gas_chgs = []
        for l in f.readlines():
            if l: gas_chgs.extend( [float(x) for x in l.strip().split()] )
    with open("sol/qout") as f:
        sol_chgs = []
        for l in f.readlines():
            if l: sol_chgs.extend( [float(x) for x in l.strip().split()] )

    assert len(gas_chgs) == len(sol_chgs)

    chgs = [a*d+b*(1-d) for a,b in zip(gas_chgs, sol_chgs)]

    atoms = ""
    from rdkit import Chem
    mol=Chem.MolFromMol2Block("sol.mol2", removeHs=False)
    conf =  mol.GetConformer()
    for i, at in enumerate(mol.GetAtoms()):
        aname = at.GetProp('_TriposAtomName')
        atype = at.GetProp('_TriposAtomType')
        achg = at.GetProp('_TriposPartialCharge')
        p = conf.GetAtomPosition(idx)
        atoms.append( mol2_line % (i+1, aname, p.x, p.y, p.z, atype, "1", "LIG", chgs[i]) )

    with open("sol.mol2") as f:
        lines = list(f.readlines())
        c = lines[2].strip().split()[0]

    acpype = "acpype -i %(lig)s -n %(nc)d -c user"
    os.system(sol % nc)
    os.system(acpype % nc)
    

 
if __name__ == "__main__":
    import argparse
    import json
    parser = argparse.ArgumentParser(description='resp')
    parser.add_argument("-l", help="ligand.gro")
    parser.add_argument('-t', type=str, help="ligand.itp")
    parser.add_argument('-m', type=str, help="ligand.mol2")
    parser.add_argument('-dp', type=str, help="dummy.pairs")
    parser.add_argument('-n', type=int, help="netcharge")
    parser.add_argument('-o', type=str, help="default ouput base name")


    args = parser.parse_args()
    #parsing the command

    with open(args.dp) as f:
        dummy_pairs = json.load(f)

    itp_block, top_block, mol2_blocks = ITPDummyReCharges(args.l, args.t, args.m, dummy_pairs=dummy_pairs)

    otop = args.o + "_LP_RESP.top"
    oitp = args.o + "_LP_RESP.itp"
    omol2 = args.o + "_LP_RESP.mol2"

    with open(otop, "w") as f:
        f.write(top_block)
    with open(oitp, "w") as f:
        f.write(itp_block)
    with open(omol2, "w") as f:
        f.write(mol2_blocks)
