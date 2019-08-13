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
