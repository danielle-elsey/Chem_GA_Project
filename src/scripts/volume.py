#!/usr/bin/env python

# Simple estimation of VdW molecular volume (in A^3)
# by Geoffrey Hutchison <geoffh@pitt.edu>
#
# Usage:
# volume.py *.g09
#
# Thanks to Wikipedia and StackExchange
# https://en.wikipedia.org/wiki/Spherical_cap
# https://math.stackexchange.com/a/297764

from __future__ import print_function

import sys
import os.path
import math
import glob

import openbabel as ob
import pybel

# vdw radii
vdw_radii = [
  # From Alvarez doi: 10.1039/C3DT50599E
  # Dalton Trans., 2013,42, 8617-8636
  # Dummy, 1st row
  0.69, 1.2, 1.43,
  # 2nd row (Li..Ne)
  2.12, 1.98, 1.91, 1.77, 1.66, 1.50, 1.46, 1.58,
  # 3rd row (Na .. Ar)
  2.50, 2.51, 2.25, 2.19, 1.90, 1.89, 1.82, 1.83,
  # 4th row (K, Ca)
  2.73, 2.62,
  # 1st row TM (Sc.. Zn)
  2.58, 2.46, 2.42, 2.45, 2.45, 2.44, 2.40, 2.40, 2.38, 2.39,
  # 4th row p-block (Ga .. Kr)
  2.32, 2.29, 1.88, 1.82, 1.86, 2.25,
  # 5th row Rb, Sr
  3.21, 2.84,
  # 2nd row TM (Y .. Cd)
  2.75, 2.52, 2.56, 2.45, 2.44, 2.46, 2.44, 2.15, 2.53, 2.49,
  # 5th row p-block (Sn .. Xe)
  2.43, 2.42, 2.47, 1.99, 2.04, 2.06,
  # 6th row Cs, Ba
  3.48, 3.03,
  # Lanthanides (La..Gd)
  2.98, 2.88, 2.92, 2.95, 2.90, 2.87, 2.83,
  # Lanthanides (Tb..Yb)
  2.79, 2.87, 2.81, 2.83, 2.79, 2.80,
  # 3rd row TM (Lu..Hg)
  2.74, 2.63, 2.53, 2.57, 2.49, 2.48, 2.41, 2.29, 2.32, 2.45,
  # 6th row p-block (Tl.. Bi)
  # 2.5 is a default here
  2.47, 2.60, 2.54, 2.5, 2.5, 2.5,
  # 7th row
  # 2.5 is a default here
  2.5, 2.5,
  # Actinides (default is now 3.0)
  2.8, 2.93, 2.88, 2.71, 2.82, 2.81, 2.83, 3.05, 3.38, 3.05, 3., 3., 3., 3.,
  # Trans-actinides
  3., 3., 3., 3., 3., 3., 3., 3., 3., 3.,
  # 7th row p-block
  3., 3., 3., 3., 3., 3.,
]



# Main Loop -- go through multiple files, e.g. volume.py *.g09
for filename in sys.argv[1:]:
#for filename in glob.iglob("*.mol"):
    extension = os.path.splitext(filename)[1][1:]
    mol = next(pybel.readfile(extension, filename))

    volume = 0.0
    # Add the volumes of all atoms
    for atom in mol:
        r = vdw_radii[atom.atomicnum]
        volume += 4.0*math.pi*r**3/3.0

    # now subtract overlapping region from bonds
    for bond in ob.OBMolBondIter(mol.OBMol):
        d = bond.GetLength()
        r1 = vdw_radii[bond.GetBeginAtom().GetAtomicNum()]
        r2 = vdw_radii[bond.GetEndAtom().GetAtomicNum()]
        r1r2 = r1 + r2
        volume -= math.pi/(12.0*d)*(r1r2 - d)**2*(d**2 + 2.0*d*(r1r2)-3.0*(r1-r2)**2)


    print(filename, "Volume: ", volume)
