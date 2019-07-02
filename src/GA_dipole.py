'''
Notes:
-each polymer is defined in the format [(numerical sequence), monomer1 index, ... monomerX index]
-population is list of these polymers
'''

import os
import subprocess
import string
import random
import math
import pybel
import numpy as np
import multiprocessing as mp
from scipy import stats
from statistics import mean
from itertools import product
from copy import deepcopy
import shlex
import pandas as pd

ob = pybel.ob

def analysis(data_file, perc):

    df = pd.read_csv(data_file, sep='\t', header=0)

    # Makes a pivot table with SMILES as the rows and columns of different row numbers
    pt = pd.pivot_table(df, index=["SMILES"], values=["Run"], columns=["Run"], fill_value=0)
    # Add across each row to get the total number of monomers and sort from largest to smallest
    pt_sum = pt.sum(1).sort_values(ascending=False)

    # Takes the values that are larger than the quartile being analyzed
    a = pt_sum[pt_sum > pt_sum.quantile(perc)]

    # Identifies the number of monomers within that top percentage which will be used in the analysis
    num_monomers = len(a.index)

    # Generates smaller pivot table which includes only the top monomers
    pt_subset = pt.loc[a.index.values]

    # Pairwise spearman correlation of all runs
    corr = pt_subset.corr(method='spearman')

    # Removes duplicate Spearman correlations and substituting with nan
    c = corr.copy()
    c.values[np.tril_indices_from(c)] = np.nan

    s = c.unstack().mean()

    # Generate ndarray of values of the Spearman correlations
    s1 = s.values

    # Removes nan values leaving just the desired  Spearman correlations
    sv = s1[~np.isnan(s1)]

    # Converts numpy array to list of values which can then be saved individually
    sv_values = list(sv)

    # Saves the values in a test file in the format percentage, Generations, Spearman Value which
    # can then be plotted. This file will just be appended to which allows many data sets to be combined into
    # one long list with all data which can then be plotted, analyzed, etc.

    for i in sv_values:
        #generations = int(sys.argv[3])
        generations = 2

        #monomers = int(sys.argv[4])
        monomers = 1234

        #with open(sys.argv[2], "a") as output_file:
        with open('spear_output.txt', 'a') as output_file:
            output_file.write("%i\t%i\t%s\t%f\n" % (monomers, generations, perc, i))

    return sv_values[0]

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

def run_geo_opt(i, polymer, smiles_list, poly_size):
    '''
    Runs geometry optimization calculation

    Parameters
    ---------
    i: int
        index of polymer
    polymer: [(#,#,#,#), A, B]
        ???
    poly_size: int
        number of monomers per polymer
    smiles_list: list
        list of monomer SMILES

    Returns
    -------
    N/A
    '''

    poly_smiles = construct_polymer_string(polymer, smiles_list, poly_size)

    # make polymer string into openbabel object
    mol = pybel.readstring('smi', poly_smiles)
    make3D(mol)

    # capture monomer indexes and numerical sequence as strings for file naming
    mono1 = str(polymer[1])
    mono2 = str(polymer[2])
    seq = polymer[0]
    # create string of actual length sequence (e.g. hexamer when poly_size = 6)
    triple_seq = seq + seq + seq
    true_length_seq = ''.join(str(triple_seq[i]) for i in range(poly_size))

    # write polymer .xyz file to containing folder
    mol.write('xyz', 'input/%s_%s_%s.xyz' % (mono1, mono2, true_length_seq), overwrite=True)

    xtb = subprocess.call('(/ihome/ghutchison/geoffh/xtb/xtb input/%s_%s_%s.xyz --opt >output/%s_%s_%s.out)' % (mono1, mono2, true_length_seq, mono1, mono2, true_length_seq), shell=True)
    save_file = subprocess.call('(cp xtbopt.xyz opt/%s_%s_%s_opt.xyz)' % (mono1, mono2, true_length_seq), shell=True)
    del_restart = subprocess.call('(rm -f *restart)', shell=True)

def find_poly_dipole(population, poly_size, smiles_list):
    '''
    Calculates dipole of polymers
    TODO: add functionality to use polarizability calculations

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
    '''

    poly_polar_list = []
    poly_dipole_list = []

    # TODO: get rid of enumerate after implementing new naming system
    for i, polymer in enumerate(population):
        run_geo_opt(i, polymer, smiles_list, poly_size)

        # capture monomer indexes and numerical sequence as strings for file naming
        mono1 = str(polymer[1])
        mono2 = str(polymer[2])
        seq = polymer[0]
        # create string of actual length sequence (e.g. hexamer when poly_size = 6)
        triple_seq = seq + seq + seq
        true_length_seq = ''.join(str(triple_seq[i]) for i in range(poly_size))

        read_output = open('output/%s_%s_%s.out' % (mono1, mono2, true_length_seq), 'r')

        # parse output file for static polarizability and dipole moment
        for line in read_output:
            # create list of tokens in line
            tokens = line.split()

            if line.startswith(" Mol. α(0)"):
                temp_polar = float(tokens[4])
                poly_polar_list.append(temp_polar)
            '''
            elif line.startswith("   full:"):
                try:
                    # dipole tensor - STILL A LIST OF STRINGS (not floats)
                    # TODO: add tensor functionality later
                    dipole_line = tokens
                    temp_dipole = float(tokens[4])
                    poly_dipole_list.append(temp_dipole)
                else:
                    poly_dipole_list.append(-10)
                # break inner for loop to avoid overwriting with other lines starting with "full"
                break
             '''

        read_output.close()

    # print("calculated dipoles:", poly_dipole_list)
    # print("calculated polarizabilities:", poly_polar_list)
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

def mutate(polymer, sequence_list, smiles_list, mono_list):
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

    # determine whether to mutate
    mut_rate = 1
    rand = random.randint(1,10)

    if rand <= (mut_rate*10):
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
        new_mono = random.randint(0, len(smiles_list) - 1)
        polymer[point] = new_mono
        # increase frequency count for monomer in mono_list
        mono_list[new_mono][1] += 1


    return polymer

def crossover_mutate(parent_list, pop_size, num_type_mono, sequence_list, smiles_list, mono_list):
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
    #num_children = pop_size - len(parent_list)
    # initialize new population with parents
    new_pop = deepcopy(parent_list)

    while len(new_pop) < pop_size:
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
        temp_child = mutate(temp_child, sequence_list, smiles_list, mono_list)

        #try to avoid duplicates in population, but prevent infinite loop if unique individual not found after so many attempts
        if temp_child in new_pop:
            pass
        else:
            new_pop.append(temp_child)

    return new_pop

def get_mono_ranks(mono_list):
    '''
    Make list of ranks of monomers by usage

    Parameters
    ---------
    mono_list: list (specific format)
        [[monomer index, frequency]]

    Returns
    -------
    mono_ranks: list
        list of monomer ranks (by usage) in order of monomer SMILES list indicies
    '''

    ordered_mono_list = deepcopy(mono_list)

    # sort based on frequencies, put in descending order (highest first)
    ordered_mono_list = sorted(ordered_mono_list, key = lambda monomer: monomer[1], reverse = True)

    # add ranks(sorted order number) to list of mono indicies and freqs(sorted order number)
    for x in range(len(ordered_mono_list)):
        ordered_mono_list[x].append(x)

    # reorder list by monomer index (lowest first)
    ordered_mono_list = sorted(ordered_mono_list, key = lambda monomer: monomer[0])

    #print("ordered_mono_list", ordered_mono_list)

    # extract list of ranks in order of monomer indicies
    mono_ranks = []
    for x in range(len(ordered_mono_list)):
        mono_ranks.append(ordered_mono_list[x][2])

    #print("mono ranks", mono_ranks)
    return mono_ranks

def get_spearman_ranks(mono_list, init_mono_indicies):
    '''
    Create nested list of numerical rank with monomer index (rank 0 = most freq momomer)
    Order by frequency (most freq monomer index first)
    Rank are initial ranks in new order

    Parameters
    ---------
    mono_list: list (specific format)
        [[monomer index, frequency]]

    Returns
    -------
    mono_spear_ranks: list (specific format)
        [[initial ranks], [final ranks]]
    '''

    init_mono_ranks = list(range(len(mono_list)))

    # make list of mono indicies in order of freq, highest first
    final_mono_indicies = make_mono_indicies_list(mono_list)

    # make list of monomer indicies and their ranks (most freq first with lowest rank)
    final_mono_ranks = []
    for x in range(len(final_mono_indicies)):
        mono_idx = final_mono_indicies[x]
        rank = init_mono_indicies.index(mono_idx)
        final_mono_ranks.append(rank)

    mono_spear_ranks = [init_mono_ranks, final_mono_ranks]

    '''
    print("initial ranks \n:", init_mono_ranks)
    print(init_mono_indicies)
    print("final ranks \n:", final_mono_ranks)
    print(final_mono_indicies)
    '''

    return mono_spear_ranks

def make_mono_indicies_list(mono_list):
    '''
    Make list of all monomer indicies ordered by frequency (highest first)

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
    ordered_mono_list = sorted(ordered_mono_list, key = lambda monomer: monomer[1], reverse = True)

    # make list of monomer indicies only (in order of descending frequency)
    mono_indicies = []
    for x in range(len(ordered_mono_list)):
        mono_indicies.append(ordered_mono_list[x][0])

    #print(mono_indicies)

    return mono_indicies

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
    # number of monomers in imported monomer list
    mono_list_size = 1234


    #true mw max for 1235MonomerList.txt w pop_size = 10 poly_size = 5
    min_std = 3434

    # property of interest (options: molecular weight 'mw', dipole moment 'dip')
    opt_property = "mw"

    # Read in monomers from input file
    read_file = open('../input_files/1235MonomerList.txt', 'r')
    # read_file = open('/ihome/ghutchison/dch45/Chem_GA_Project/input_files/1235MonomerList.txt', 'r')


    # create list of monomer SMILES strings
    # assumes input file has one monomer per line
    smiles_list = []
    for line in read_file:
        # grab first token from each line
        temp = line.split()[0]
        smiles_list.append(temp)

    read_file.close()

    #create monomer list [(mono index 1, frequency), (mono index 2, frequency),...]

    mono_list = []
    for x in range(len(smiles_list)):
        mono_list.append([x, 0])

    # create all possible numerical sequences for given number of monomer types
    sequence_list = find_sequences(num_type_mono)


    # create inital population as list of polymers
    population = []
    counter = 0
    while counter < pop_size:
        # for polymer in range(pop_size):
        temp_poly = []

        # select sequence type for polymerf
        poly_seq = sequence_list[random.randint(0, len(sequence_list) - 1)]
        temp_poly.append(poly_seq)

        # select monomer types for polymer
        for num in range(num_type_mono):
            # randomly select a monomer index
            poly_monomer = random.randint(0, len(smiles_list) - 1)
            temp_poly.append(poly_monomer)
            # increase frequency count for monomer in mono_list
            mono_list[poly_monomer][1] += 1

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
        # print("population:", population)
        poly_property_list = find_poly_dipole(population, poly_size, smiles_list)
    else:
        print("Error: opt_property not recognized. trace:main:initial pop properties")

    # set initial values for min, max, and avg polymer weights
    min_test = min(poly_property_list)
    max_test = max(poly_property_list)
    avg_test = (max_test - min_test) / 2.0

    # create new output files
    analysis_file = open('gens_analysis.txt', 'w+')
    population_file = open('gens_population.txt', 'w+')
    values_file = open('gens_values', 'w+')

    # write files headers
    analysis_file.write('min, max, avg, spearman, \n')
    population_file.write('polymer populations \n')
    values_file.write('%s values \n' % (opt_property))

    #capture initial population data
    analysis_file.write('%f, %f, %f, n/a, \n' % (min_test, max_test, avg_test))

    for polymer in population:
        # capture monomer indexes and numerical sequence as strings for population file
        mono1 = str(polymer[1])
        mono2 = str(polymer[2])
        seq = polymer[0]
        # create string of actual length sequence (e.g. hexamer when poly_size = 6)
        triple_seq = seq + seq + seq
        true_length_seq = ''.join(str(triple_seq[i]) for i in range(poly_size))

        population_file.write('%s_%s_%s, ' % (mono1, mono2, true_length_seq))

    population_file.write('\n')

    for value in poly_property_list:
        values_file.write('%f, ' % (value))
    values_file.write('\n')

    # Loop

    # check for convergence among top 30% (or top 8, whichever is larger) candidates between 5 generations
    perc = 0.3
    n = int(mono_list_size*perc)
    print("n=", n)
    #n = int(1234*perc)
    #if n < 8:
        #n = 8

    method_file = open('spear_method_test.txt', 'a')
    #method_file.write('pop_size  %s\npoly_size   %s\nnum_type_mono   %s\nperc    %f\n' % (pop_size, poly_size, num_type_mono, perc))

    # initialize generation counter
    gen_counter = 0

    # initialize convergence counters for each of 3 speaman methods
    spear_counter = 0
    lam_counter = 0
    three_counter = 0

    # initialize trips to stop printing after convergence
    spear_trip = 0
    lam_trip = 0
    three_trip = 0

    #while spear_counter < 10:
    for x in range(300):
        gen_counter += 1

        # create sorted monomer list with most freq first
        gen1 = make_mono_indicies_list(mono_list)

        # make file to test Ilana's spearman method
        spear_file = open('spear_input.txt', 'w')
        spear_file.write('Run\tSMILES\tCount\n')
        for x in range(len(mono_list)):
            spear_file.write('%i\t%s\t%i\n' % (gen_counter-1, str(mono_list[x][0]), mono_list[x][1]))

        # Selection - select heaviest (best) 50% of polymers as parents
        population = parent_select(opt_property, population, poly_property_list)

        # Crossover & Mutation - create children to repopulate bottom 50% of polymers in population
        population = crossover_mutate(population, pop_size, num_type_mono, sequence_list, smiles_list, mono_list)

        # calculate desired polymer property
        if opt_property == "mw":
            poly_property_list = find_poly_mw(population, poly_size, mw_list)
        elif opt_property == "dip":
            poly_property_list = find_poly_dipole(population, poly_size, smiles_list)
        else:
            print("Error: opt_property not recognized. trace:main:loop pop properties")

        # record representative generation properties
        min_test = min(poly_property_list)
        max_test = max(poly_property_list)
        avg_test = mean(poly_property_list)

        # create sorted monomer list with most freq first
        gen2 = make_mono_indicies_list(mono_list)

        # calculate monomer idx method Spearman correlation coefficient
        spear = stats.spearmanr(gen1[:n], gen2[:n])[0]

        # make file and calculate Ilana's method Spearman correlation coefficient
        for x in range(len(mono_list)):
            spear_file.write('%i\t%s\t%i\n' % (gen_counter, str(mono_list[x][0]), mono_list[x][1]))
        lamarck_spear = analysis('spear_input.txt', 1-perc)

        # make second rank list for method 3
        spearman_ranks = get_spearman_ranks(mono_list, gen1)
        three_spear = stats.spearmanr(spearman_ranks[0][:n], spearman_ranks[1][:n])[0]

        # capture monomer indexes and numerical sequence as strings for population file
        analysis_file.write('%f, %f, %f, %f, \n' % (min_test, max_test, avg_test, spear))

        for polymer in population:
            # capture monomer indexes and numerical sequence as strings for population file
            mono1 = str(polymer[1])
            mono2 = str(polymer[2])
            seq = polymer[0]
            # create string of actual length sequence (e.g. hexamer when poly_size = 6)
            triple_seq = seq + seq + seq
            true_length_seq = ''.join(str(triple_seq[i]) for i in range(poly_size))

            population_file.write('%s_%s_%s, ' % (mono1, mono2, true_length_seq))

        population_file.write('\n')

        for value in poly_property_list:
            values_file.write('%f, ' % (value))

        values_file.write('\n')


        # keep track of number of successive generations meeting Spearman criterion
        if spear > 0.92:
            spear_counter += 1
        else:
            spear_counter = 0

        if lamarck_spear > 0.92:
            lam_counter += 1
        else:
            lam_counter = 0


        if three_spear > 0.92:
            three_counter += 1
        else:
            three_counter = 0

        # check for convergence
        if spear_counter == 10 and spear_trip == 0:
            spear_gen = gen_counter
            spear_max = max_test
            spear_coeff = spear

            print('method 1 convergence reached')
            print('generation:', gen_counter)
            print('top values:', poly_property_list[:4])
            #print(population[:4])
            print('spearman coeff:', spear)
            print('\n')
            spear_trip = 1

        if lam_counter == 10 and lam_trip == 0:
            lam_gen = gen_counter
            lam_max = max_test
            lam_coeff = lamarck_spear

            print('method 2 convergence reached')
            print('generation:', gen_counter)
            print('top values:', poly_property_list[:4])
            #print(population[:4])
            print('spearman coeff:', lamarck_spear)
            print('\n')
            lam_trip = 1

        if three_counter == 10 and three_trip == 0:
            three_gen = gen_counter
            three_max = max_test
            three_coeff = three_spear

            print('method 3 convergence reached')
            print('generation:', gen_counter)
            print('top values:', poly_property_list[:4])
            #print(population[:4])
            print('spearman coeff:', three_spear)
            print('\n')
            three_trip = 1

    #method_file.write('convergence generation\t\t\ttop MW\t\t\tspearman coeff\t\t')
    method_file.write('%s\t%s\t%s\t\t%f\t%f\t%f\t\t%f\t%f\t%f\n' % (spear_gen, lam_gen, three_gen, spear_max, lam_max, three_max, spear_coeff, lam_coeff, three_coeff))



    # close all output files
    analysis_file.close()
    population_file.close()
    values_file.close()
    method_file.close()

    '''
    # print("LAMARCK SPEAR CONVERGED\n generation: %i\n propery value: %i" % (gen_counter, max_test))
    #for x in range(len(converge_generations_lam)):
        #print(converge_generations_lam[x])


    #Print out SMILES strings meeting MW criteria
    polymer_smiles_list = []
    for polymer in population:
        polymer_smiles_list.append(construct_polymer_string(polymer, smiles_list, poly_size))
    print(polymer_smiles_list)
    '''

if __name__ == '__main__':
    main()
