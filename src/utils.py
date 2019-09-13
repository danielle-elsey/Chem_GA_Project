'''
Utility script for general functions such as making
smiles strings and files
'''
import pybel
import random
from copy import deepcopy
from itertools import product

ob = pybel.ob


def make3D(mol):
    '''
    Makes the mol object from SMILES 3D

    Parameters
    ---------
    mol: object
        pybel molecule object
    '''
    # make mol object 3D and add hydrogens
    pybel._builder.Build(mol.OBMol)
    mol.addh()

    ff = pybel._forcefields["mmff94"]
    success = ff.Setup(mol.OBMol)
    if not success:
        ff = pybel._forcefields["uff"]
        success = ff.Setup(mol.OBMol)
        if not success:
            sys.exit("Cannot set up forcefield")

    ff.ConjugateGradients(100, 1.0e-3)
    ff.WeightedRotorSearch(100, 25)
    ff.ConjugateGradients(250, 1.0e-4)

    ff.GetCoordinates(mol.OBMol)


def make_file_name(polymer, poly_size):
    '''
    Makes file name for a given polymer

    Parameters
    ---------
    polymer: list (specific format)
        [(#,#,#,#), A, B]

    Returns
    -------
    file_name: str
        polymer file name (w/o extension) showing monomer indicies and full sequence
        e.g. 100_200_101010 for a certain hexamer
    '''

    # capture monomer indexes as strings for file naming
    mono1 = str(polymer[1])
    mono2 = str(polymer[2])

    # make string of actual length sequence (e.g. hexamer when poly_size = 6)
    seq = list(polymer[0] * ((poly_size // 4) + 1))[:poly_size]
    true_length_seq = ''.join(map(str, seq))

    # make file name string
    file_name = '%s_%s_%s' % (mono1, mono2, true_length_seq)

    return file_name


def make_polymer_str(polymer, smiles_list, poly_size):
    '''
    Constructs polymer string from monomers, adds standard end groups [amino and nitro]

    Parameters
    ---------
    polymer: list (specific format)
        [(#,#,#,#), A, B]
    smiles_list: list
        list of all possible monomer SMILES
    poly_size: int
        number of monomers per polymer

    Returns
    -------
    poly_string: str
        polymer SMILES string
    '''
    poly_string = ''

    # add amino end group
    poly_string = poly_string + "N"

    # cycle over monomer sequence until total number of monomers in polymer is reached
    cycle = 0
    for i in range(poly_size):
        seq_cycle_location = i - cycle * (len(polymer[0]))
        seq_monomer_index = polymer[0][seq_cycle_location]
        if seq_cycle_location == (len(polymer[0]) - 1):
            cycle += 1

        # find monomer index from given sequence location in polymer
        monomer_index = polymer[seq_monomer_index + 1]
        poly_string = poly_string + smiles_list[monomer_index]

    # add nitro end group
    poly_string = poly_string + "N(=O)=O"

    return poly_string


def find_sequences(num_mono_species):
    '''
    Finds all possible sequences

    Parameters
    ---------
    num_mono_species: int
        number of monomer species in each polymer (e.g. copolymer = 2)

    Returns
    -------
    numer_seqs: list
        all possible sequences as a numerical list
    '''
    # find length of each sequence
    seq_len = num_mono_species**2

    # find all possible sequences as numerical lists [start index 0]
    # (use cartesian product of each type of monomer over sequence length)
    numer_seqs = list(product(range(num_mono_species), repeat=seq_len))

    return numer_seqs


def sort_mono_indicies_list(mono_list):
    '''
    Makes list of all monomer indicies ordered by frequency (highest first)

    Parameters
    ---------
    mono_list: list (specific format)
        [[monomer index, frequency]]

    Returns
    -------
    mono_indicies: list
        list of monomer indicies ordered by frequency (highest first)
    '''

    ordered_mono_list = deepcopy(mono_list)

    # sort based on frequencies, put in descending order (highest first)
    ordered_mono_list = sorted(
        ordered_mono_list, key=lambda monomer: monomer[1], reverse=True)

    # make list of monomer indicies only (in order of descending frequency)
    mono_indicies = []
    for x in range(len(ordered_mono_list)):
        mono_indicies.append(ordered_mono_list[x][0])

    return mono_indicies


def weighted_random_pick(mono_list):
    '''
    Selects random monomer based on weighted (frequency) system

    Parameters
    ---------
    mono_list: list (specific format)
        [[monomer index, frequency]]

    Returns
    -------
    mono_idx: int
        index of selected randomly chosen monomer
    '''

    ordered_mono_list = deepcopy(mono_list)

    # sort based on frequencies, put in ASCENDING order (LOWEST first)
    ordered_mono_list = sorted(
        ordered_mono_list, key=lambda monomer: monomer[1])

    # calculate sum of all weights [frequencies]
    sum = 0
    for idx in ordered_mono_list:
        freq = idx[1]
        sum += freq

    # generate random number from 0 to sum, lower endpoint included
    rand_num = random.randint(0, sum - 1)

    # loop over list of sorted weights, subtract weights from random number until random number is less than next weight
    for idx in ordered_mono_list:
        mono_idx = idx[0]
        freq = idx[1]
        if rand_num < freq:
            return mono_idx
        else:
            rand_num -= freq


def find_volume(filename):
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

    extension = filename.split('.')[1]
    print('extension', extension)
    mol = next(pybel.readfile(extension, filename))
    print('filename', filename)

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

    return volume


def make_vol_list(population, poly_size):
    vol_list = []
    for poly in population:
        temp_file_name = utils.make_file_name(poly, poly_size)
        print(temp_file_name)
        # try:
        volume = find_volume('%s.xyz' % temp_file_name)
        vol_list.append(volume)
        # catch:
        #     print('Error in finding volumes')

    return vol_list
