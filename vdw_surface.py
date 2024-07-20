from __future__ import division, absolute_import, print_function

import numpy as np


"""
A sript to generate van der Waals surface of molecules.
"""

# Van der Waals radii (in angstrom) are taken from GAMESS.
# GAMESS https://raw.githubusercontent.com/streaver91/gamess/master/source/prplib.src
vdw_r = {'H': 1.20, 'HE': 1.20,
         'LI': 1.37, 'BE': 1.45, 'B': 1.45, 'C': 1.50,
         'N': 1.50, 'O': 1.40, 'F': 1.35, 'NE': 1.30,
         'NA': 1.57, 'MG': 1.36, 'AL': 1.24, 'SI': 1.17,
         'P': 1.80, 'S': 1.75, 'CL': 1.70,
         'BR': 1.92, 'I': 2.11} # Bondi radii; BR vdw rdii = 1.80 in GAMESS

def surface(n):
    """Computes approximately n points on unit sphere. Code adapted from GAMESS.

    Parameters
    ----------
    n : int
        approximate number of requested surface points

    Returns
    -------
    ndarray
        numpy array of xyz coordinates of surface points
    """

    u = []
    eps = 1e-10
    nequat = int(np.sqrt(np.pi*n))
    nvert = int(nequat/2)
    nu = 0
    for i in range(nvert+1):
        fi = np.pi*i/nvert
        z = np.cos(fi)
        xy = np.sin(fi)
        nhor = int(nequat*xy+eps)
        if nhor < 1:
            nhor = 1
        for j in range(nhor):
            fj = 2*np.pi*j/nhor
            x = np.cos(fj)*xy
            y = np.sin(fj)*xy
            if nu >= n:
                return np.array(u)
            nu += 1
            u.append([x, y, z])
    return np.array(u)

def vdw_surface(coordinates, elements, scale_factor, density, input_radii):
    """Computes points outside the van der Waals surface of molecules.

    Parameters
    ----------
    coordinates : ndarray
        cartesian coordinates of the nuclei, in units of angstrom
    elements : list
        The symbols (e.g. C, H) for the atoms
    scale_factor : float
        The points on the molecular surface are set at a distance of
        scale_factor * vdw_radius away from each of the atoms.
    density : float
        The (approximate) number of points to generate per square angstrom
        of surface area. 1.0 is the default recommended by Kollman & Singh.
    input_radii : dict
        dictionary of user's defined VDW radii

    Returns
    -------
    radii : dict
        A dictionary of scaled VDW radii
    surface_points : ndarray
        array of the coordinates of the points on the surface

    """
    radii = {}
    surface_points = []
    # scale radii
    vdw_keys_up = [k.upper() for k in vdw_r.keys()]
    for i in elements:
        if i in radii.keys():
            continue
        if i in input_radii.keys():
            radii[i] = input_radii[i] * scale_factor
        elif i in vdw_r.keys():
            radii[i] = vdw_r[i] * scale_factor
        elif i.upper() in vdw_r.keys():
            radii[i] = vdw_r[i.upper()] * scale_factor
        else:
            raise KeyError('%s is not a supported element; ' %i
                         + 'use the "VDW_RADII" option to add '
                         + 'its van der Waals radius.')
    # loop over atomic coordinates
    for i in range(len(coordinates)):
        # calculate approximate number of ESP grid points
        n_points = int(density * 4.0 * np.pi* np.power(radii[elements[i]], 2))
        # generate an array of n_points in a unit sphere around the atom
        dots = surface(n_points)
        # scale the unit sphere by the VDW radius and translate
        dots = coordinates[i] + radii[elements[i]] * dots
        for j in range(len(dots)):
            save = True
            for k in range(len(coordinates)):
                if i == k:
                    continue
                # exclude points within the scaled VDW radius of other atoms
                d = np.linalg.norm(dots[j] - coordinates[k])
                if d < radii[elements[k]]:
                    save = False
                    break
            if save:
                surface_points.append(dots[j])

    return np.array(surface_points), radii


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
            coords.append( (float(cs[1]), float(cs[2]), float(cs[3])) )
    coords = np.array(coords)
    return coords, symbols

if __name__ == "__main__":
    import os
    import sys

    if len(sys.argv) > 2:
        init_conf = sys.argv[1]
        bohr = True
    elif len(sys.argv) > 1:
        init_conf = sys.argv[1]
        bohr = False
    else:
        print("Usage: python vdw_surface.py example.xyz")
        exit(1)

    bohr_to_angstrom = 0.52917721092

    # Get the points at which we're going to calculate the ESP
    options = {
            "VDW_POINT_DENSITY": 5.0,
            "VDW_SCALE_FACTORS": [1.4, 1.6, 1.8, 2.0],
            "VDW_RADII": {},
            }
    points = []
    coordinates, symbols = ReadXYZ(sys.argv[1])
    for scale_factor in options['VDW_SCALE_FACTORS']:
        shell, radii = vdw_surface(coordinates, symbols, scale_factor,
                                    options['VDW_POINT_DENSITY'], options['VDW_RADII'])
        points.append(shell)
    points = np.concatenate(points)
    if bohr:
        # for XTB, use AU
        points /= bohr_to_angstrom
    np.savetxt('grid.dat', points, fmt='%15.10f')
