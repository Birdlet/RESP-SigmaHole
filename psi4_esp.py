import sys
import os
import time
from math import exp
import numpy as np

from vdw_surface import ReadXYZ, vdw_surface



psi_hartree2kcalmol = 627.5095 # Hartree to kcal mol$^{-1}$ conversion factor
psi_hartree2kJmol = 2625.500   # Hartree to kilojoule mol$^{-1}$ conversion factor
psi_bohr2angstroms = 0.52917720859 # Bohr to Angstroms conversion factor
temprature = 298.15  # Kalven
minus_1_div_RT = -0.0019872 # mol/(kcal*K)

method = "GFN2"
forcefield = "MMFF94"
sample = False
n_sample = 50
n_dft = 4

# https://psi4.github.io/psi4docs/master/basissets_byelement.html
# https://psicode.org/psi4manual/4.0b5/dft_byfunctional.html

frozen = """
set optking {
  frozen_dihedral = ("
%s  ")
}
"""

inp_prev ="""
# Any line starting with the # character is a comment line
#! Sample HF/cc-pVDZ H2O computation

memory 2 GB


molecule CONF {
 %(charge)s 1
"""


high_dft_basis = {
        "OBASIS": "def2-tzvp",
        "BASIS": "def2-tzvp",
        "DFT": "m06-2x-d3zero2b",
        }
low_dft_basis = {
        "OBASIS": "6-31g*\n    assign I  lanl2dz",
        "BASIS": "6-311g*\n    assign I  def2-tzvp",
        "DFT": "b3lyp-d3",
        }

original_basis = {
        "OBASIS": "6-31g*\n    assign I  lanl2dz",
        "BASIS": "6-311g*\n    assign I  def2-tzvp",
        "DFT": "HF",
        }


inp_post_opt ="""

set {
    freeze_core True
}

basis {
    assign %(OBASIS)s
}

optimize('%(DFT)s')
CONF.save_xyz_file("opt.xyz", 1)

"""

inp_post_gas ="""

set {
    freeze_core True
}

basis {
    assign %(BASIS)s
}

# set cubeprop_tasks ['esp']

E, wfn = prop('%(DFT)s', properties=["GRID_ESP", "GRID_FIELD"], return_wfn=True)
oeprop(wfn, 'DIPOLE', 'QUADRUPOLE', title='DIPOLE+QUADRUPOLE')
# cubeprop(wfn)
"""

inp_post_sol ="""

set {
    freeze_core True
    pcm true
    pcm_scf_type total
}

basis {
    assign %(BASIS)s
}

pcm = {
   Units = Angstrom
   Medium {
   SolverType = IEFPCM
   Solvent = Water
   }

   Cavity {
   RadiiSet = UFF
   Type = GePol
   Scaling = False
   Area = 0.3
   Mode = Implicit
   }
}
# set cubeprop_tasks ['esp']

E, wfn = prop('%(DFT)s', properties=["GRID_ESP", "GRID_FIELD"], return_wfn=True)
oeprop(wfn, 'DIPOLE', 'QUADRUPOLE', title='DIPOLE+QUADRUPOLE')

# cubeprop(wfn)
"""


xtb_xyz = """%s
$constrain
  force constant=100.0
%s$end"""

# dihedral: 3, 4, 1, 7, auto

def Boltzman(energy):
    qs = []
    for e in energy:
        qs.append( exp( e*minus_1_div_RT*temprature ) )
    return [100*q/sum(qs) for q in qs]

def write_gesp(xyz="opt.xyz", grid="grid.dat", esp="grid_esp.dat", gesp="Lig.gesp.psi4"):
    coords, symbols = ReadXYZ(xyz)

    xyzs = []
    with open(grid) as f:
        for l in f.readlines():
            xyz = l.strip().split()
            if xyz:
                xyz = float(xyz[0])/psi_bohr2angstroms, float(xyz[1])/psi_bohr2angstroms, float(xyz[2])/psi_bohr2angstroms
                xyzs.append(xyz)
    esps = []
    with open(esp) as f:
        for l in f.readlines():
            if l:
                esps.append(float(l.strip()))

    assert len(esps) == len(xyzs)
    with open(gesp, "w") as f:
        global charge
        f.write(" ESP FILE - ATOMIC UNIT\n")
        f.write(" CHARGE = %3s - MULTIPLICITY =   1\n" % charge)
        f.write(" ATOMIC COORDINATES AND ESP CHARGES. #ATOMS = %8d\n" % len(symbols))

        for i in range(len(symbols)):
            f.write("  %-7s%15.8e %15.8e %15.8e %15.8e\n" % (symbols[i], coords[i][0]/psi_bohr2angstroms,
                coords[i][1]/psi_bohr2angstroms, coords[i][2]/psi_bohr2angstroms, 0.0))

        f.write(" DIPOLE MOMENT:\n\n")
        f.write(" TRACELESS QUADRUPOLE MOMENT:\n\n\n")
        f.write(" ESP VALUES AND GRID POINT COORDINATES. #POINTS = %7d\n" %len(esps))
        for e, c in zip(esps, xyzs):
            f.write(" %15.8e %15.8e %15.8e %15.8e\n" % (e, c[0], c[1], c[2]))

if __name__ == "__main__":
    higher_level = False
    original_level = False

    if len(sys.argv) > 2:
        init_conf = sys.argv[1]
        charge = sys.argv[2]
    else:
        print("Usage: python opiskit.py example.xyz charge[int]")
        exit(1)

    start = time.time()

    print("-"*60)

    frozen_dihedrals = ""

    if higher_level:
        inp_post_opt = inp_post_opt % high_dft_basis
        inp_post_sol = inp_post_sol % high_dft_basis
        inp_post_gas = inp_post_gas % high_dft_basis
    elif original_level:
        inp_post_opt = inp_post_opt % original_basis
        inp_post_sol = inp_post_sol % original_basis
        inp_post_gas = inp_post_gas % original_basis
    else:
        inp_post_opt = inp_post_opt % low_dft_basis
        inp_post_sol = inp_post_sol % low_dft_basis
        inp_post_gas = inp_post_gas % low_dft_basis
    inp_prev = inp_prev % locals()

    SOLV = False

    if True:
        with open(init_conf) as f:
            lines = f.readlines()
            natoms = int(lines[0].strip())
            atoms = "".join(lines[2:2+natoms])

        print( "Optimize Conformation")
        with open("psi4_opt.in", "w") as f:
            f.write(inp_prev)
            f.write(atoms)
            f.write("}\n")
            f.write(inp_post_opt)
        print("psi4 psi4_opt.in opt.out -n 40")
        os.system("psi4 psi4_opt.in opt.out -n 40")

        with open("opt.out") as f:
            lines = f.read().split("\n")
            for j, l in enumerate(lines):
                if l.startswith( "    Total Energy =" ):
                    ei = j
                elif l.startswith("    Geometry"):
                    gi = j
                elif l.startswith("PsiException: Could not converge geometry optimization in 50 iterations"):
                    gi = gi-2

        try:
            qmenergy = float(lines[ei].split()[-1]) * psi_hartree2kcalmol
        except:
            qmenergy = 0
        #optxyz = "%d\n energy: %f hartree\n" % (natoms-2, float(lines[ei].split()[-1]))
        #for l in lines[gi+4:gi+4+natoms-2]:
        #   xs = l.strip().split()
        #optxyz += "%-6s%24.12f%24.12f%24.12f\n" % (xs[0], float(xs[1]), float(xs[2]), float(xs[3]))
        #with open("opt.xyz", "w") as f:
        #    f.write(optxyz)
        if not os.path.exists("opt.xyz"):
            print("Optimization Failed! USE gfb-xtb2 to optimize!")
            os.system("xtb %s -o -P 40" % init_conf)
            os.system("mv xtbopt.xyz opt.xyz")

        os.system('sed -i "s/CL/Cl/g" opt.xyz')
        os.system('sed -i "s/BR/Br/g" opt.xyz')
        os.system("obabel -ixyz opt.xyz -omol2 -O opt.mol2 2>/dev/null")
        with open("opt.xyz") as f:
            lines = f.readlines()
            natoms = int(lines[0].strip())
            atoms = "".join(lines[2:2+natoms])

        options = {
            "VDW_POINT_DENSITY": 5.0,
            "VDW_SCALE_FACTORS": [1.4, 1.6, 1.8, 2.0],
            # "VDW_SCALE_FACTORS": [1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0],
            "VDW_RADII": {},
            }
        points = []
        coordinates, symbols = ReadXYZ("opt.xyz")
        for scale_factor in options['VDW_SCALE_FACTORS']:
            shell, radii = vdw_surface(coordinates, symbols, scale_factor,
                                    options['VDW_POINT_DENSITY'], options['VDW_RADII'])
            points.append(shell)
        points = np.concatenate(points)
        np.savetxt('grid.dat', points, fmt='%15.10f')


        print("Gas Phase SPE")
        with open("psi4_gas.in", "w") as f:
            f.write(inp_prev)
            f.write(atoms)
            f.write("}\n")
            f.write(inp_post_gas)
        print("psi4 psi4_gas.in opt.out -n 40")
        os.system("psi4 psi4_gas.in gas.out -n 40")

        os.system("mv grid_esp.dat grid_esp_gas.dat")
        os.system("mv grid_field.dat grid_field_gas.dat")
        write_gesp(xyz="opt.xyz", grid="grid.dat", esp="grid_esp_gas.dat", gesp="gas.gesp")

        with open("gas.out") as f:
            lines = f.read().split("\n")
            for j, l in enumerate(lines):
                if l.startswith( "    Total Energy =" ):
                    ei = j
                elif l.startswith("    Geometry"):
                    gi = j
                elif l.startswith("PsiException: Could not converge geometry optimization in 50 iterations"):
                    gi = gi-2
        spe_gas_energy = float(lines[ei].split()[-1]) * psi_hartree2kcalmol

        spe_sol_energy = 0.0
        print("Solvated Phase SPE")
        with open("psi4_sol.in", "w") as f:
            f.write(inp_prev)
            f.write(atoms)
            f.write("}\n")
            f.write(inp_post_sol)
        if SOLV:
            print("psi4 psi4_sol.in opt.out -n 40")
            os.system("psi4 psi4_sol.in sol.out -n 40")

            os.system("mv grid_esp.dat grid_esp_sol.dat")
            os.system("mv grid_field.dat grid_field_sol.dat")
            write_gesp(xyz="opt.xyz", grid="grid.dat", esp="grid_esp_sol.dat", gesp="sol.gesp")

            with open("sol.out") as f:
                lines = f.read().split("\n")
                for j, l in enumerate(lines):
                    if l.startswith( "    Total Energy =" ):
                        ei = j
                    elif l.startswith("    Geometry"):
                        gi = j
                    elif l.startswith("PsiException: Could not converge geometry optimization in 50 iterations"):
                        gi = gi-2
            spe_sol_energy = float(lines[ei].split()[-1]) * psi_hartree2kcalmol

        print("Done!")
        print("Gas Phase Energy: %8.4f  Solv Phase Energy: %8.4f" % (spe_gas_energy, spe_sol_energy)) 

    print("-"*60)
    print("Time Consuming: %d s" % (time.time()-start))

