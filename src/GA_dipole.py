'''
Notes:
-each polymer is defined in the format [(numerical sequence), monomer1 index, ... monomerX index]
-population is list of these polymers
'''

import string
import random
import math
import pybel#print SMILES string of max polymer
import numpy as np
from statistics import mean
from itertools import product
from copy import deepcopy

ob = pybel.ob


def make3D(mol):
    '''
    Makes the mol object from SMILES 3D

    Parameters
    ---------
    mol: object
        pybel molecule object
    '''
    pybel._builder.Build(mol.OBMol)
    mol.addh()


def find_sequences(num_type_mono):
    '''
    Finds all possible sequences

    Parameters
    ---------
    num_type_mono: int
        number of types of monomers in each polymer

    Returns
    -------
    numer_seqs: list
        all possible sequences as a numerical list
    '''
    # find length of each sequence
    seq_len = num_type_mono**2

    # find all possible sequences as numerical lists [start index 0]
    # (use cartesian product of each type of monomer over sequence length)
    numer_seqs = list(product(range(num_type_mono), repeat=seq_len))

    return numer_seqs


def find_poly_mw(population, poly_size, mw_list):
    '''
    Calculates molecular weight of polymers

    Parameters
    ---------
    population: list
        list of polymers in population
    poly_size: int
        number of monomers per polymer
    mw_list: list
        list of molecular weights

    Returns
    -------
    poly_mw_list: list
        list of polymer molecular weights
    '''
    poly_mw_list = []
    for polymer in population:
        temp_mw = 0
        # counter for cycling over monomer sequence
        cycle = 0
        # loop over total monomers in each polymer
        for i in range(poly_size):
            # cycle over monomer sequence until total number of monomers is reached
            seq_cycle_location = i - cycle * (len(polymer[0]))
            seq_monomer_index = polymer[0][seq_cycle_location]
            if seq_cycle_location == (len(polymer[0]) - 1):
                cycle += 1
            # find index in polymer of monomer index (from smiles_list) from given sequence value
            monomer_index = polymer[seq_monomer_index + 1]
            temp_mw += mw_list[monomer_index]
        poly_mw_list.append(temp_mw)
    return poly_mw_list


def find_poly_dipole(population, poly_size, smiles_list):
    '''
    Calculates dipole of polymers

    Parameters
    ---------
    population: list
        list of polymers in population
    poly_size: int
        number of monomers per polymer
    smiles_list: list
        list of monomer SMILES

    Returns
    -------
    poly_dipole_list: list
        list of polymer dipoles
    '''#print SMILES string of max polymer
    poly_dipole_list = []
    for polymer in population:
        poly_smiles = construct_polymer_string(polymer, smiles_list, poly_size)
        mol = pybel.readstring("smi", poly_smiles)
        make3D(mol)
        cm = ob.OBChargeModel_FindType('mmff94')
        cm.ComputeCharges(mol.OBMol)
        diptensor = cm.GetDipoleMoment(mol.OBMol)
        dipole = math.sqrt(diptensor.GetX()**2 +
                           diptensor.GetY()**2 + diptensor.GetZ()**2)
        poly_dipole_list.append(dipole)
    return poly_dipole_list


def make_pop_readable(population, sequence_list, smiles_list):
    '''
    Returns population with polymers in human readable format:
    [[alpha sequence], 'monomer SMILES string 1', .., 'monomer SMILES string X']

    Parameters
    ---------
    population: list
        list of polymers in population
    sequence_list: list
        list of sequences
    smiles_list: list
        list of monomer SMILES
    Returns
    -------
    readable_poly: list
        readable list of polymer information
    '''
    # create dict mapping numbers [start index 0] to uppercase English alphabet
    alpha_dict = dict(zip(range(26), string.ascii_uppercase))

    readable_poly = []
    for polymer in population:
        temp_poly = []

        # find numerical sequence from index of sequence_list
        temp_numer_seq = polymer[0]
        # translate numerical sequence into alphabetic sequence
        temp_alpha_seq = []
        for number in temp_numer_seq:
            temp_alpha_seq.append(alpha_dict[number])
        temp_poly.append(temp_alpha_seq)

        for monomer_index in polymer[1:]:
            # find monomer SMILES string from index of smiles_list
            temp_mono_smiles = smiles_list[monomer_index]
            temp_poly.append(temp_mono_smiles)

        readable_poly.append(temp_poly)

    return readable_poly


def fitness_fn(opt_property, population, poly_property_list):
    '''
    Fitness function ranking current population members
    Rank 0 = best

    Parameters
    ---------
    opt_property: str
        property being optimized
    population: list
        list of polymers in population
    poly_prop_list: list
        list of properties of the polymers

    Returns
    -------
    ranked_indicies: list
        ranked list of polymers with first entry being best
    '''
    # rank based on highest property value = best
    if opt_property == "mw" or "dip":
        # gives indicies in ascending order
        ranked_indicies = list(np.argsort(poly_property_list))
        # print(ranked_indicies)
        # reverse list so highest property value = 0th
        # ranked_indicies[::-1]
        ranked_indicies.reverse()
        # print('\n\n\n', ranked_indicies)
    else:
        print("Error: opt_property not recognized. trace:fitness_fn")
    return ranked_indicies


def parent_select(opt_property, population, poly_property_list):
    '''
    Return top half of the population

    Parameters
    ---------
    opt_property: str
        property being optimized
    population: list
        list of polymers in population
    poly_prop_list: list
        list of properties of the polymers

    Returns
    -------
    parent_list: list
        top half of population
    '''
    parent_count = int(len(population) / 2)
    parent_list = []

    fitness_list = fitness_fn(opt_property, population, poly_property_list)

    for x in range(parent_count):
        parent_list.append(population[fitness_list[x]])

    return parent_list


def mutate(polymer, sequence_list, smiles_list):
    '''
    Mutates polymer

    Parameters
    ---------
    polymer: [(#,#,#,#), A, B]
        ???
    sequence_list: list
        list of sequences
    smiles_list: list
        list of monomer SMILES

    Returns
    -------
    polymer: ???
        polymer string ???)
    '''

    # determine whether to mutate (mutation rate is 1/3*1/3 = 1/9 or ~11%)
    rand1 = random.randint(1,3)
    rand2 = random.randint(1,3)

    if rand1 == rand2:
        pass
    else:
        return polymer

    # choose point of mutation (sequence or specific monomer)
    point = random.randint(0, len(polymer) - 1)
    # replace sequence
    if point == 0:
        polymer[point] = sequence_list[random.randint(
            0, len(sequence_list) - 1)]
    # or replace specific monomer
    else:
        polymer[point] = random.randint(0, len(smiles_list) - 1)

    return polymer


def crossover_mutate(parent_list, pop_size, num_type_mono, sequence_list, smiles_list):
    '''
    Crossover of polymers

    Parameters
    ---------
    parent_list: list
        list of parents
    pop_size: int
        number of polymers in each generation
    num_type_mono: int
        number of types of monomer per polymer
    sequence_list: list
        list of sequences
    smiles_list: list
        list of monomer SMILES

    Returns
    -------
    new_pop: list
        crossover population list
    '''
    # calculate number of children to generate
    num_children = pop_size - len(parent_list)
    # initialize new population with parents
    new_pop = deepcopy(parent_list)

    while len(new_pop) < num_children:
        #for child in range(num_children):
        # randomly select two parents (as indexes from parent list) to cross
        parent_a = random.randint(0, len(parent_list) - 1)
        parent_b = random.randint(0, len(parent_list) - 1)

        # ensure parents are unique indiviudals
        if len(parent_list) > 1:
            while parent_b == parent_a:
                parent_b = random.randint(0, len(parent_list) - 1)

        # determine number of monomers taken from parent A
        num_mono_a = random.randint(1, num_type_mono)

        # randomly determine which parent's sequence will be used
        par_seq = random.randint(0, 1)

        # create hybrid child
        temp_child = []
        # give child appropriate parent's sequence
        if par_seq == 0:
            temp_child.append(parent_list[parent_a][0])
        else:
            temp_child.append(parent_list[parent_b][0])

        # give child first half monomers from A, second half from B
        for monomer in range(1, num_mono_a + 1):
            temp_child.append(parent_list[parent_a][monomer])
        if num_mono_a < num_type_mono:
            for monomer in range(num_mono_a + 1, num_type_mono + 1):
                temp_child.append(parent_list[parent_b][monomer])

        #new_pop.append(temp_child)
        #counter += 1

        # mutate
        temp_child = mutate(temp_child, sequence_list, smiles_list)

        #try to avoid duplicates in population, but prevent infinite loop if unique individual not found after so many attempts
        if temp_child in new_pop:
            pass
        else:
            new_pop.append(temp_child)

    return new_pop


def construct_polymer_string(polymer, smiles_list, poly_size):
    '''
    Construction of polymer from monomers, adds standard end groups [amino and nitro]

    Parameters
    ---------
    polymer: list (specific format)
        [(#,#,#,#), A, B]
    smiles_list: list
        list of monomer SMILES
    poly_size: int
        number of monomers per polymer

    Returns
    -------
    poly_string: str
        polymer SMILES string
    '''
    poly_string = ''

    #add amino end group
    poly_string = poly_string + "N"

    cycle = 0
    # loop over total monomers in each polymer
    for i in range(poly_size):
        # cycle over monomer sequence until total number of monomers is reached
        seq_cycle_location = i - cycle * (len(polymer[0]))
        seq_monomer_index = polymer[0][seq_cycle_location]
        if seq_cycle_location == (len(polymer[0]) - 1):
            cycle += 1
        # find index in polymer of monomer index (from smiles_list) from given sequence value
        monomer_index = polymer[seq_monomer_index + 1]
        poly_string = poly_string + smiles_list[monomer_index]

    #add nitro end group
    poly_string = poly_string + "N(=O)=O"

    return poly_string


def main():
    # number of polymers in population
    pop_size = 32
    # number of monomers per polymer
    poly_size = 6
    # number of types of monomers in each polymer
    num_type_mono = 2
    # minimum difference between two successive generation values to declare convergence
    convergence_std = 0.001

    min_std = 3434

    # property of interest (options: molecular weight 'mw', dipole moment 'dip')
    opt_property = "dip"

    # Read in monomers from input file
    #read_file = open('../input_files/1235MonomerList.txt', 'r')
    read_file = open('ihome/ghutchison/dch45/Chem_GA_Project/input_files/1235MonomerList.txt', 'r')


    # create list of monomer SMILES strings
    # assumes input file has one monomer per line
    smiles_list = []
    for line in read_file:
        # grab first token from each line
        temp = line.split()[0]
        smiles_list.append(temp)

    read_file.close()

    # create all possible numerical sequences for given number of monomer types
    sequence_list = find_sequences(num_type_mono)

    # create inital population as list of polymers
    population = []
    counter = 0
    while counter < pop_size:
        # for polymer in range(pop_size):
        temp_poly = []

        # select sequence type for polymer
        poly_seq = sequence_list[random.randint(0, len(sequence_list) - 1)]
        temp_poly.append(poly_seq)

        # select monomer types for polymer
        for num in range(num_type_mono):
            poly_monomer = random.randint(0, len(smiles_list) - 1)
            temp_poly.append(poly_monomer)

        # add polymer to population
        # check for duplication
        if temp_poly in population:
            pass
        else:
            population.append(temp_poly)
            counter += 1

    # find initial population properties
    if opt_property == "mw":
        # calculate MW for each monomer
        mw_list = []
        for monomer in smiles_list:
            temp_mono = pybel.readstring("smi", monomer).molwt
            mw_list.append(temp_mono)

        # Calculate polymer molecular weights
        poly_property_list = find_poly_mw(population, poly_size, mw_list)
    elif opt_property == "dip":
        # calculate dipole moments for each polymer
        poly_property_list = find_poly_dipole(population, poly_size, smiles_list)
    else:
        print("Error: opt_property not recognized. trace:main:initial pop properties")

    # set initial values for min, max, and avg polymer weights
    min_test = min(poly_property_list)
    max_test = max(poly_property_list)
    avg_test = (max_test - min_test) / 2.0

    # create new output read_file
    write_file = open("ga_polymer_output.txt", "w+")
    write_file.write("min max avg\n")

    #Loop
    last_min = 0
    #print(last_min, min_test)
    #while (min_test < min_std):
    #while(abs(last_min-min_test) > convergence_std):
    for x in range(20):

        #last_min = min_test

        #Selection - select best 50% as parents
        #select heaviest (best) 50% of polymers as parents
        population = parent_select(opt_property, population, poly_property_list)

        #Crossover - create children to repopulate bottom 50% of polymers in population
        population = crossover_mutate(population, pop_size, num_type_mono, sequence_list, smiles_list)

        #find starting place for children (second half of population)
        #child_start = int(len(population)/2)+1
        #Mutation
        #for polymer in population[child_start:]:
            #polymer = mutate(polymer, sequence_list, smiles_list)

        if opt_property == "mw":
            poly_property_list = find_poly_mw(population, poly_size, mw_list)
        elif opt_property == "dip":
            poly_property_list = find_poly_dipole(population, poly_size, smiles_list)
        else:
            print("Error: opt_property not recognized. trace:main:loop pop properties")

        #find minimum polymer weight
        min_test = min(poly_property_list)
        max_test = max(poly_property_list)
        avg_test = mean(poly_property_list)

        write_file.write("{} {} {}\n".format(min_test, max_test, avg_test))
        print(min_test, max_test, avg_test)

        #print SMILES string of max polymer
        #index = poly_property_list.index(max_test)
        #print(construct_polymer_string(population[index], smiles_list, poly_size))

        #print(min_test)


    write_file.close()

    '''
    #Print out SMILES strings meeting MW criteria
    polymer_smiles_list = []
    for polymer in population:
        polymer_smiles_list.append(construct_polymer_string(polymer, smiles_list, poly_size))
    print(polymer_smiles_list)
    '''

if __name__ == '__main__':
    main()
