#! Electrostatic potential and electric field evaluated on a grid around water.

bohr_to_angstrom = 0.52917721092

def ReadXYZ(fname):
    coords = []
    symbols = []
    with open(fname) as f:
        lines = list(f.readlines())
        natoms = int(lines[0].strip())
        atom_lines = lines[2:2+natoms]
        for at in atom_lines:
            cs = at.strip().split()
            symbols.append(cs[0])
            coords.append( (float(cs[1])/bohr_to_angstrom, float(cs[2])/bohr_to_angstrom, float(cs[3])/bohr_to_angstrom) )
    return coords, symbols

coords, symbols = ReadXYZ("opt.xyz")

xyzs = []
esps = []
try:
    # psi4 coordinates in Angstrom
    with open("grid_esp.dat") as f:
        for l in f.readlines():
            if l:
                esps.append(float(l.strip()))
    with open("grid.dat") as f:
        for l in f.readlines():
            xyz = l.strip().split()
            if xyz:
                xyz = float(xyz[0])/bohr_to_angstrom, float(xyz[1])/bohr_to_angstrom, float(xyz[2])/bohr_to_angstrom
                xyzs.append(xyz)
except:
    # xtb all in AU
    with open("xtb_esp.dat") as f:
        for l in f.readlines():
            if l:
                cs = l.strip().split()
                esps.append(float(cs[-1]))
                xyzs.append((float(cs[0]), float(cs[1]), float(cs[2])))

assert len(esps) == len(xyzs)
with open("Lig.gesp", "w") as f:
    f.write(" ESP FILE - ATOMIC UNIT\n")
    f.write(" CHARGE =   0 - MULTIPLICITY =   0\n")
    f.write(" ATOMIC COORDINATES AND ESP CHARGES. #ATOMS = %8d\n" % len(symbols))

    for i in range(len(symbols)):
        f.write("  %-7s%15.8e %15.8e %15.8e %15.8e\n" % (symbols[i], coords[i][0], coords[i][1], coords[i][2], 0.0))

    f.write(" DIPOLE MOMENT:\n\n")
    f.write(" TRACELESS QUADRUPOLE MOMENT:\n\n\n")
    f.write(" ESP VALUES AND GRID POINT COORDINATES. #POINTS = %7d\n" %len(esps))
    for e, c in zip(esps, xyzs):
        f.write(" %15.8e %15.8e %15.8e %15.8e\n" % (e, c[0], c[1], c[2]))
