import sys
import os
infile =sys.argv[1]
outfile = sys.argv[2]
psi_bohr2angstroms = 0.52917720859 # Bohr to Angstroms conversion factor
ofname = os.path.split(outfile)[-1]


with open(infile) as f:
    lines = list(f.readlines())
    natoms = int(lines[0].strip().split()[0])
    nesps = int(lines[0].strip().split()[1])
    mol2s = "@<TRIPOS>MOLECULE\n%s\n %d 0 1\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n" % (ofname, nesps)
    for i, l in enumerate(lines[1+natoms:1+natoms+nesps]):
        cs = l.split()
        # "  %5d H          -4.1596   -0.0665   -0.0000 N.2   1    UNL        -0.5079"
        mol2_line = "   %4d %-8s  %8.4f%10.4f%10.4f %-5s %-4s %-8s  %8.4f\n" % (i+1, "H", 
                        psi_bohr2angstroms*float(cs[1]), psi_bohr2angstroms*float(cs[2]), psi_bohr2angstroms*float(cs[3]),
                        "H", "1", "ESP", float(cs[0]))
        mol2s += mol2_line

with open(outfile, "w") as f:
    f.write(mol2s)