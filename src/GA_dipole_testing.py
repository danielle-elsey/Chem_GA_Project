'''
Notes:
-each polymer is defined in the format [(numerical sequence), monomer1 index, ... monomerX index]
-population is list of these polymers
'''

import os
import shutil
import subprocess
import random
import pybel
import numpy as np
from scipy import stats
from statistics import mean
from copy import deepcopy
import pickle

from scipy.constants import Boltzmann

import utils
import scoring


def find_poly_mw(population, poly_size, smiles_list):
    '''
    Calculates molecular weight of polymers

    Parameters
    ---------
    population: list
        list of polymers in population
    poly_size: int
        number of monomers per polymer
    smiles_list: list
        list of all possible monomer SMILES

    Returns
    -------
    poly_mw_list: list
        list of molecular weights of polymers in population
    '''
    poly_mw_list = []
    for polymer in population:
        # make polymer into SMILES string
        poly_smiles = utils.make_polymer_str(polymer, smiles_list, poly_size)
        # make polymer string into pybel molecule object
        mol = pybel.readstring('smi', poly_smiles)

        # add mw of polymer to list
        poly_mw = mol.molwt
        poly_mw_list.append(poly_mw)

    return poly_mw_list


def run_geo_opt(polymer, poly_size, smiles_list):
    '''
    Runs geometry optimization calculation on given polymer

    Parameters
    ---------
    polymer: list (specific format)
        [(#,#,#,#), A, B]
    poly_size: int
        number of monomers per polymer
    smiles_list: list
        list of all possible monomer SMILES

    '''
    # make file name string w/ convention monoIdx1_monoIdx2_fullNumerSequence
    file_name = utils.make_file_name(polymer, poly_size)

    # if output file already exists, skip xTB
    exists = os.path.isfile('output/%s.out' % (file_name))
    if exists:
        print("output file existed")
        return

    # make polymer into SMILES string
    poly_smiles = utils.make_polymer_str(polymer, smiles_list, poly_size)

    # make polymer string into pybel molecule object
    mol = pybel.readstring('smi', poly_smiles)
    utils.make3D(mol)

    # write polymer .xyz file to containing folder
    mol.write('xyz', 'input/%s.xyz' % (file_name), overwrite=True)

    # make directory to run xtb in for the polymer
    mkdir_poly = subprocess.call('(mkdir %s)' % (file_name), shell=True)

    # run xTB geometry optimization
    xtb = subprocess.call('(cd %s && /ihome/ghutchison/geoffh/xtb/xtb ../input/%s.xyz --opt >../output/%s.out)' %
                          (file_name, file_name, file_name), shell=True)

    save_opt_file = subprocess.call(
        '(cp %s/xtbopt.xyz opt/%s_opt.xyz)' % (file_name, file_name), shell=True)

    # delete xtb run directory for the polymer
    del_polydir = subprocess.call('(rm -r %s)' % (file_name), shell=True)


def find_elec_prop(population, poly_size, smiles_list):
    '''
    Calculates dipole moment and polarizability of each polymer in population
    TODO: add dipole tensor functionality
    TODO: update parser to catch failures/errors in output file

    Parameters
    ---------
    population: list
        list of polymers in population
    poly_size: int
        number of monomers per polymer
    smiles_list: list
        list of all possible monomer SMILES

    Returns
    -------
    elec_prop_lists: list
        nested list of [list of polymer dipole moments, list of polymer polarizabilities]
    '''

    poly_polar_list = []
    poly_dipole_list = []

    # run xTB geometry optimization
    #nproc = 8
    for polymer in population:
        run_geo_opt(polymer, poly_size, smiles_list)

    # parse xTB output files
    for polymer in population:
         # make file name string w/ convention monoIdx1_monoIdx2_fullNumerSequence
        file_name = utils.make_file_name(polymer, poly_size)
        # count number of successful parsed properties
        num_succ_reads = 0

        # parse output file for static polarizability and dipole moment, check for xTB success before recording values
        read_output = open('output/%s.out' % (file_name), 'r')
        for line in read_output:
            # create list of tokens in line
            tokens = line.split()

            if line.startswith(" Mol. α(0)"):
                temp_polar = float(tokens[4])
                num_succ_reads += 1
            elif line.startswith("molecular dipole:"):
                # iterate down to line with dipole value - necessary b/c dip value on non unique line
                next(read_output)
                next(read_output)
                new_tokens = next(read_output).split()

                # dipole tensor - STILL A LIST OF STRINGS (not floats)
                # TODO: add tensor functionality later
                temp_dipole = float(new_tokens[4])
                num_succ_reads += 1
            if line.startswith("#  - GEOMETRY OPTIMIZATION FAILED!"):
                num_succ_reads = 0

        read_output.close()

        if num_succ_reads == 2:
            poly_dipole_list.append(temp_dipole)
            poly_polar_list.append(temp_polar)
        # This is a catch to move files to the failed folder if polarizability
        # .. and dipole are not parsed.
        else:
            poly_polar_list.append(-10)
            poly_dipole_list.append(-10)
            move_fail_file = subprocess.call(
                '(mv output/%s.out failed/%s.out)' % (file_name, file_name), shell=True)
        # reset the number of successful reads
        num_succ_reads = 0

    # make nested list of dipole moment and polarizability lists
    elec_prop_lists = [poly_dipole_list, poly_polar_list]

    return elec_prop_lists


def fitness_fn(opt_property, poly_property_list1, poly_property_list2, population, poly_size):
#     properties = []
#     for list in poly_property_list:
#         properties.append(list)

    if opt_property in ('mw', 'dip', 'pol'):
        ranked_indicies = scoring.simple_descending(poly_property_list1)
    elif opt_property == 'dip_pol':
        ranked_indicies = scoring.comb_dip_pol(poly_property_list1, poly_property_list2, 1, 1, population, poly_size)
    else:
        print('Error: opt_property not recognized; Traceback: fitnesss_fn')
    return(ranked_indicies)


# def fitness_fn(opt_property, poly_property_list):
#     '''
#     Ranks polymers in population based on associated property values
#     Rank 0 = best
#
#     Parameters
#     ---------
#     opt_property: str
#         property being optimized
#     poly_prop_list: list
#         list of property values corresponding to polymers in population
#
#     Returns
#     -------
#     ranked_indicies: list
#         list of ranked polymer indicies with first entry being best
#     '''
#     # rank polymers based on highest property value = best
#     if opt_property in ('mw', 'dip', 'pol') :
#         # make list of indicies of polymers in population, sorted based on property values
#         ranked_indicies = list(np.argsort(poly_property_list))
#         # reverse list so highest property value = 0th
#         ranked_indicies.reverse()
#     else:
#         print("Error: opt_property not recognized. trace:fitness_fn")
#
#     return ranked_indicies
#
#
# def fitness_fn_multi(opt_property_1, prop_list_1, opt_property_2, prop_list_2):
#     if opt_property_1 in ('mw', 'dip', 'pol') :
#         # find ranks of properties
#         get_ranks = list(stats.rankdata(prop_list_1, method = 'ordinal'))
#         # convert ranks so that 1st = highest value
#         ranks_1 = []
#         for x in range (len(get_ranks)):
#             ranks_1.append(len(get_ranks)-get_ranks[x])
#     else:
#         print("Error: opt_property not recognized. trace:fitness_fn")
#
#     # make sorted list of polymer indicies based on second property
#     if opt_property_2 in ('mw', 'dip', 'pol') :
#         # find ranks of properties
#         get_ranks = list(stats.rankdata(prop_list_2, method = 'ordinal'))
#         # convert ranks so that 1st = highest value
#         ranks_2 = []
#         for x in range (len(get_ranks)):
#             ranks_2.append(len(get_ranks)-get_ranks[x])
#     else:
#         print("Error: opt_property not recognized. trace:fitness_fn")
#
#     # average ranks from both properties for each polymer
#     avg_ranks= []
#     for x in range(len(ranks_1)):
#         temp_index = mean([float(ranks_1[x]), float(ranks_2[x])])
#         avg_ranks.append(temp_index)
#
#     # make list of indicies of polymers in population, sorted based on averaged ranks (0th = highest/best)
#     ranked_indicies = list(np.argsort(avg_ranks))
#
#     return ranked_indicies
#

def parent_select(opt_property, population, poly_property_list, poly_size):
    '''
    Finds top half of population

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
    # find number of parents (half of population)
    parent_count = int(len(population) / 2)

    # make list of ranked polymer indicies
#     if opt_property != 'dip_pol':
#         fitness_list = fitness_fn(opt_property, poly_property_list)
#     else:
#         fitness_list = fitness_fn_multi('dip', poly_property_list[0], 'pol', poly_property_list[1])

    fitness_list = fitness_fn('dip_pol', poly_property_list[0], poly_property_list[1], population, poly_size)

    # make list of top polymers
    parent_list = []
    for x in range(parent_count):
        parent_list.append(population[fitness_list[x]])

    return parent_list


def mutate(polymer, sequence_list, smiles_list, mono_list):
    '''
    Creates point mutation in given polymer if polymer is selected for mutation

    Parameters
    ---------
    polymer: list (specific format)
        [(#,#,#,#), A, B]
    sequence_list: list
        list of sequences
    smiles_list: list
        list of monomer SMILES

    Returns
    -------
    polymer: list (specific format)
        [(#,#,#,#), A, B]
    '''

    # set mutation rate
    mut_rate = 0.4

    # determine whether to mutate based on mutation rate
    rand = random.randint(1, 10)
    if rand <= (mut_rate * 10):
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
        #new_mono = utils.weighted_random_pick(mono_list)
        polymer[point] = new_mono
        # increase frequency count for monomer in mono_list
        mono_list[new_mono][1] += 1

    return polymer


def crossover_mutate(parent_list, pop_size, poly_size, num_mono_species, sequence_list, smiles_list, mono_list):
    '''
    Performs crossover and mutation functions on given population
    TODO: fix possible duplication problem after mutation

    Parameters
    ---------
    parent_list: list
        list of parent polymers
    pop_size: int
        number of polymers in each generation
    num_mono_species: int
        number of monomer species in each polymer (e.g. copolymer = 2)
    sequence_list: list
        list of sequences
    smiles_list: list
        list of all possible monomer SMILES

    Returns
    -------
    new_pop: list
        population list after crossover and mutation
    '''

    # initialize new population with parents
    new_pop = deepcopy(parent_list)
    new_pop_str = []
    for parent in new_pop:
        parent_str = utils.make_polymer_str(parent, smiles_list, poly_size)
        new_pop_str.append(parent_str)

    # loop until enough children have been added to reach population size
    while len(new_pop) < pop_size:

        # randomly select two parents (as indexes from parent list) to cross
        parent_a = random.randint(0, len(parent_list) - 1)
        parent_b = random.randint(0, len(parent_list) - 1)

        # ensure parents are unique indiviudals
        if len(parent_list) > 1:
            while parent_b == parent_a:
                parent_b = random.randint(0, len(parent_list) - 1)

        # determine number of monomers taken from parent A
        num_mono_a = random.randint(1, num_mono_species)

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
        if num_mono_a < num_mono_species:
            for monomer in range(num_mono_a + 1, num_mono_species + 1):
                temp_child.append(parent_list[parent_b][monomer])

        # give child opportunity for mutation
        temp_child = mutate(temp_child, sequence_list, smiles_list, mono_list)

        temp_child_str = utils.make_polymer_str(
            temp_child, smiles_list, poly_size)

        # try to avoid duplicates in population, but prevent infinite loop if unique individual not found after so many attempts
        # TODO: fix possible duplication problem after mutation
        if temp_child_str in new_pop_str:
            pass
        else:
            new_pop.append(temp_child)
            new_pop_str.append(temp_child_str)

    return new_pop


def init_gen(pop_size, poly_size, num_mono_species, opt_property, perc, smiles_list):
    '''
    Initializes parameter, creates population, and runs initial generation

    Parameters
    ----------
    pop_size: int
        number of polymers in each generation
    poly_size: int
        number of monomers per polymer
    num_mono_species: int
        number of monomer species in each polymer (e.g. copolymer = 2)
    opt_property: str
        property being optimized
    perc: float
        percentage of number of monomers to compare with Spearman calculation
    smiles_list: list
        list of all possible monomer SMILES

    Returns
    -------
    params: list (specific format)
        [pop_size, poly_size, num_mono_species, opt_property, smiles_list, sequence_list, mono_list, population, poly_property_list, n, gen_counter, spear_counter, prop_value_counter]
    '''

    # create all possible numerical sequences for given number of monomer types
    sequence_list = utils.find_sequences(num_mono_species)

    n = int(len(smiles_list) * perc)
    # n_05 = int(len(smiles_list) * .05)
    # n_10 = int(len(smiles_list) * .10)
    # n_15 = int(len(smiles_list) * .15)

    # initialize generation counter
    gen_counter = 1

    # initialize convergence counter
    spear_counter = 0
    prop_value_counter = 0

    # create monomer frequency list [(mono index 1, frequency), (mono index 2, frequency),...]
    mono_list = []
    for x in range(len(smiles_list)):
        mono_list.append([x, 0])

    # create inital population as list of polymers
    population = []
    population_str = []
    counter = 0
    while counter < pop_size:
        # for polymer in range(pop_size):
        temp_poly = []

        # select sequence type for polymer
        poly_seq = sequence_list[random.randint(0, len(sequence_list) - 1)]
        temp_poly.append(poly_seq)

        # select monomer types for polymer
        for num in range(num_mono_species):
            # randomly select a monomer index
            poly_monomer = random.randint(0, len(smiles_list) - 1)
            temp_poly.append(poly_monomer)
            # increase frequency count for monomer in mono_list
            mono_list[poly_monomer][1] += 1

        # make SMILES string of polymer
        temp_poly_str = utils.make_polymer_str(
            temp_poly, smiles_list, poly_size)

        # add polymer to population
        # check for duplication - use str for comparison to avoid homopolymer, etc. type duplicates
        if temp_poly_str in population_str:
            pass
        else:
            population.append(temp_poly)
            population_str.append(temp_poly_str)
            counter += 1

    # find initial population properties
    if opt_property == 'mw':
        # calculate polymer molecular weights
        poly_property_list = find_poly_mw(population, poly_size, smiles_list)
    elif opt_property == 'dip':
        # initialize list of polarizabilities
        polar_list = []
        # calculate electronic properties for each polymer
        elec_prop_list = find_elec_prop(population, poly_size, smiles_list)
        poly_property_list = elec_prop_list[0]
        polar_list = elec_prop_list[1]
    elif opt_property == 'pol':
        # initialize list of dipole moments
        dip_list = []
        # calculate electronic properties for each polymer
        elec_prop_list = find_elec_prop(population, poly_size, smiles_list)
        poly_property_list = elec_prop_list[1]
        dip_list = elec_prop_list[0]
    elif opt_property == 'dip_pol':
        # calculate electronic properties for each polymer
        elec_prop_list = find_elec_prop(population, poly_size, smiles_list)
        # store list of dipole moments and list of polarizabilities
        poly_property_list = []
        poly_property_list.append(elec_prop_list[0])
        poly_property_list.append(elec_prop_list[1])
    else:
        print("Error: opt_property not recognized. trace:main:initial pop properties")

    # set initial values for min, max, and avg polymer weights
    if opt_property != 'dip_pol':
        min_test = min(poly_property_list)
        max_test = max(poly_property_list)
        avg_test = mean(poly_property_list)
    else:
        min_test = [min(poly_property_list[0]), min(poly_property_list[1])]
        max_test = [max(poly_property_list[0]), max(poly_property_list[1])]
        avg_test = [mean(poly_property_list[0]), mean(poly_property_list[1])]

    if opt_property == 'dip':
        compound = utils.make_file_name(
            population[poly_property_list.index(max_test)], poly_size)
        polar_val = polar_list[poly_property_list.index(max_test)]

    if opt_property == 'pol':
        compound = utils.make_file_name(
            population[poly_property_list.index(max_test)], poly_size)
        dip_val = dip_list[poly_property_list.index(max_test)]

    # create new output files
    analysis_file = open('gens_analysis.txt', 'w+')
    population_file = open('gens_population.txt', 'w+')
    values_file = open('gens_values.txt', 'w+')
    if opt_property == 'dip':
        dip_polar_file = open('gens_dip_polar.txt', 'w+')
    if opt_property == 'pol':
        polar_dip_file = open('gens_polar_dip.txt', 'w+')
    if opt_property == 'dip_pol':
        multi_file = open('gens_multi.txt', 'w+')
    #spear_file = open('gens_spear.txt', 'w+')


    # write files headers
    if opt_property != 'dip_pol':
        analysis_file.write('min, max, avg, spearman, \n')
        population_file.write('polymer populations \n')
        values_file.write('%s values \n' % (opt_property))
    else:
        analysis_file.write('dip_min, dip_max, dip_avg, pol_min, pol_max, pol_avg, spearman, \n')
        population_file.write('polymer populations \n')
        values_file.write('dip values -break- pol values \n')
    if opt_property == 'dip':
        dip_polar_file.write('compound, gen, dipole, polar \n')
    if opt_property == 'pol':
        polar_dip_file.write('compound, gen, polar, dip \n')
    if opt_property == 'dip_pol':
        multi_file.write('min_dip_term, max_dip_term, med_dip_term, min_pol_term, max_pol_term, med_pol_term, \n')
    #spear_file.write('gen, spear_05, spear_10, spear_15 \n')

    # capture initial population data
    if opt_property != 'dip_pol':
        analysis_file.write('%f, %f, %f, n/a, \n' % (min_test, max_test, avg_test))
    else:
        analysis_file.write('%f, %f, %f, %f, %f, %f, n/a, \n' % (min_test[0], max_test[0], avg_test[0], min_test[1], max_test[1], avg_test[1]))

    if opt_property == 'dip':
        dip_polar_file.write('%s, %d, %f, %f, \n' %
                             (compound, 1, max_test, polar_val))
    if opt_property == 'pol':
        polar_dip_file.write('%s, %d, %f, %f, \n' %
                             (compound, 1, max_test, dip_val))
    if opt_property == 'dip_pol':
        # determine and write to file properties of best polymer, worst polymer and "median" polymer (i.e. max dip and pol values are for SAME "best" polymer)
        fitness_list = fitness_fn('dip_pol', poly_property_list[0], poly_property_list[1], population, poly_size)

        min_polymer = population[fitness_list[len(fitness_list)-1]]
        max_polymer = population[fitness_list[0]]
        median = int((len(fitness_list)-1)/2)
        med_polymer = population[fitness_list[median]]

        # find volumes of specified polymers
        file_name_min = utils.make_file_name(min_polymer, poly_size)
        vol_min = utils.find_volume(file_name_min)

        file_name_max = utils.make_file_name(max_polymer, poly_size)
        vol_max = utils.find_volume(file_name_max)

        file_name_med = utils.make_file_name(med_polymer, poly_size)
        vol_med = utils.find_volume(file_name_med)

        # find dipole moments of specified polymers
        mu_min = poly_property_list[0][population.index(min_polymer)]
        mu_max = poly_property_list[0][population.index(max_polymer)]
        mu_med = poly_property_list[0][population.index(med_polymer)]

        # calculate dipole terms for specified polymers
        min_dip_term = mu_min**2/(3*Boltzmann*298*vol_min)
        max_dip_term = mu_max**2/(3*Boltzmann*298*vol_max)
        med_dip_term = mu_med**2/(3*Boltzmann*298*vol_med)

        # find polarizabilities for specified polymers
        alpha_min = poly_property_list[1][population.index(min_polymer)]
        alpha_max = poly_property_list[1][population.index(max_polymer)]
        alpha_med = poly_property_list[1][population.index(med_polymer)]

        # calculate polarizability terms for specified polymers
        min_pol_term = alpha_min/vol_min
        max_pol_term = alpha_max/vol_max
        med_pol_term = alpha_med/vol_med

        multi_file.write('%f, %f, %f, %f, %f, %f, \n' % (min_dip_term, max_dip_term, med_dip_term, min_pol_term, max_pol_term, med_pol_term))

    #spear_file.write('1, n/a, n/a, n/a, \n')

    # write polymer population to file
    for polymer in population:
        poly_name = utils.make_file_name(polymer, poly_size)
        population_file.write('%s, ' % (poly_name))
    population_file.write('\n')


    if opt_property != 'dip_pol':
        for value in poly_property_list:
            values_file.write('%f, ' % (value))
        values_file.write('\n')
    else:
        for item in poly_property_list:
            for value in item:
                values_file.write('%f, ' % (value))
            values_file.write('-break-, ')
        values_file.write('\n')

    # close all output files
    analysis_file.close()
    population_file.close()
    values_file.close()
    if opt_property == 'dip':
        dip_polar_file.close()
    if opt_property == 'pol':
        polar_dip_file.close()
    if opt_property == 'dip_pol':
        multi_file.close()
    #spear_file.close()

    # make backup copies of output files
    shutil.copy('gens_analysis.txt', 'gens_analysis_copy.txt')
    shutil.copy('gens_population.txt', 'gens_population_copy.txt')
    shutil.copy('gens_values.txt', 'gens_values_copy.txt')
    if opt_property == 'dip':
        shutil.copy('gens_dip_polar.txt', 'gens_dip_polar_copy.txt')
    if opt_property == 'pol':
        shutil.copy('gens_polar_dip.txt', 'gens_polar_dip_copy.txt')
    if opt_property == 'dip_pol':
        shutil.copy('gens_multi.txt', 'gens_multi_copy.txt')
    #shutil.copy('gens_spear.txt', 'gens_spear_copy.txt')

    params = [pop_size, poly_size, num_mono_species, opt_property, smiles_list, sequence_list,
              mono_list, population, poly_property_list, n, gen_counter, spear_counter, prop_value_counter]
    return(params)


def next_gen(params):
    '''
    Runs one post-initial generation

    Parameters
    ---------
    params: list (specific format)
        [pop_size, poly_size, num_mono_species, opt_property, smiles_list, sequence_list, mono_list, population, poly_property_list, n, gen_counter, spear_counter, prop_value_counter]

    Returns
    -------
    params: list (specific format)
        [pop_size, poly_size, num_mono_species, opt_property, smiles_list, sequence_list, mono_list, population, poly_property_list, n, gen_counter, spear_counter, prop_value_counter]
    '''
    # unpack parameters
    pop_size = params[0]
    poly_size = params[1]
    num_mono_species = params[2]
    opt_property = params[3]
    smiles_list = params[4]
    sequence_list = params[5]
    mono_list = params[6]
    population = params[7]
    poly_property_list = params[8]
    n = params[9]
    gen_counter = params[10]
    spear_counter = params[11]
    prop_value_counter = params[12]

    # open output files
    analysis_file = open('gens_analysis.txt', 'a+')
    population_file = open('gens_population.txt', 'a+')
    values_file = open('gens_values.txt', 'a+')
    if opt_property == 'dip':
        dip_polar_file = open('gens_dip_polar.txt', 'a+')
    if opt_property == 'pol':
        polar_dip_file = open('gens_polar_dip.txt', 'a+')
    if opt_property == 'dip_pol':
        multi_file = open('gens_multi.txt', 'a+')
    #spear_file = open('gens_spear.txt', 'a+')

    gen_counter += 1

    if opt_property != 'dip_pol':
        max_init = max(poly_property_list)

    # create sorted monomer list with most freq first
    gen1 = utils.sort_mono_indicies_list(mono_list)

    # Selection - select heaviest (best) 50% of polymers as parents
    population = parent_select(
        opt_property, population, poly_property_list, poly_size)

    # Crossover & Mutation - create children to repopulate bottom 50% of polymers in population
    population = crossover_mutate(
        population, pop_size, poly_size, num_mono_species, sequence_list, smiles_list, mono_list)

    # calculate desired polymer property
    if opt_property == "mw":
        poly_property_list = find_poly_mw(
            population, poly_size, smiles_list)
    elif opt_property == "dip":
        elec_prop_list = find_elec_prop(population, poly_size, smiles_list)
        poly_property_list = elec_prop_list[0]
        polar_list = elec_prop_list[1]
    elif opt_property == 'pol':
        elec_prop_list = find_elec_prop(population, poly_size, smiles_list)
        poly_property_list = elec_prop_list[1]
        dip_list = elec_prop_list[0]
    elif opt_property == 'dip_pol':
        # calculate electronic properties for each polymer
        elec_prop_list = find_elec_prop(population, poly_size, smiles_list)
        # store list of dipole moments and list of polarizabilities
        poly_property_list = []
        poly_property_list.append(elec_prop_list[0])
        poly_property_list.append(elec_prop_list[1])
    else:
        print("Error: opt_property not recognized. trace:main:loop pop properties")

    # record representative generation properties
    if opt_property != 'dip_pol':
        min_test = min(poly_property_list)
        max_test = max(poly_property_list)
        avg_test = mean(poly_property_list)
    else:
        min_test = [min(poly_property_list[0]), min(poly_property_list[1])]
        max_test = [max(poly_property_list[0]), max(poly_property_list[1])]
        avg_test = [mean(poly_property_list[0]), mean(poly_property_list[1])]



    if opt_property == 'dip':
        compound = utils.make_file_name(
            population[poly_property_list.index(max_test)], poly_size)
        polar_val = polar_list[poly_property_list.index(max_test)]

    if opt_property == 'pol':
        compound = utils.make_file_name(
            population[poly_property_list.index(max_test)], poly_size)
        dip_val = dip_list[poly_property_list.index(max_test)]

    # create sorted monomer list with most freq first
    gen2 = utils.sort_mono_indicies_list(mono_list)

    # calculate Spearman correlation coefficient for begin and end sorted monomer lists
    spear = stats.spearmanr(gen1[:n], gen2[:n])[0]

    # spear_05 = stats.spearmanr(gen1[:n_05], gen2[:n_05])[0]
    # spear_10 = stats.spearmanr(gen1[:n_10], gen2[:n_10])[0]
    # spear_15 = stats.spearmanr(gen1[:n_15], gen2[:n_15])[0]


    # capture monomer indexes and numerical sequence as strings for population file
    if opt_property != 'dip_pol':
        analysis_file.write('%f, %f, %f, n/a, \n' % (min_test, max_test, avg_test))
    else:
        analysis_file.write('%f, %f, %f, %f, %f, %f, n/a, \n' % (min_test[0], max_test[0], avg_test[0], min_test[1], max_test[1], avg_test[1]))

    if opt_property == 'dip':
        dip_polar_file.write('%s, %d, %f, %f, \n' %
                             (compound, 1, max_test, polar_val))
    if opt_property == 'pol':
        polar_dip_file.write('%s, %d, %f, %f, \n' %
                             (compound, 1, max_test, dip_val))
    if opt_property == 'dip_pol':
        # determine and write to file properties of best polymer, worst polymer and "median" polymer (i.e. max dip and pol values are for SAME "best" polymer)
        fitness_list = fitness_fn('dip_pol', poly_property_list[0], poly_property_list[1], population, poly_size)

        min_polymer = population[fitness_list[len(fitness_list)-1]]
        max_polymer = population[fitness_list[0]]
        median = int((len(fitness_list)-1)/2)
        med_polymer = population[fitness_list[median]]

        # find volumes of specified polymers
        file_name_min = utils.make_file_name(min_polymer, poly_size)
        vol_min = utils.find_volume(file_name_min)

        file_name_max = utils.make_file_name(max_polymer, poly_size)
        vol_max = utils.find_volume(file_name_max)

        file_name_med = utils.make_file_name(med_polymer, poly_size)
        vol_med = utils.find_volume(file_name_med)

        # find dipole moments of specified polymers
        mu_min = poly_property_list[0][population.index(min_polymer)]
        mu_max = poly_property_list[0][population.index(max_polymer)]
        mu_med = poly_property_list[0][population.index(med_polymer)]

        # calculate dipole terms for specified polymers
        min_dip_term = mu_min**2/(3*Boltzmann*298*vol_min)
        max_dip_term = mu_max**2/(3*Boltzmann*298*vol_max)
        med_dip_term = mu_med**2/(3*Boltzmann*298*vol_med)

        # find polarizabilities for specified polymers
        alpha_min = poly_property_list[1][population.index(min_polymer)]
        alpha_max = poly_property_list[1][population.index(max_polymer)]
        alpha_med = poly_property_list[1][population.index(med_polymer)]

        # calculate polarizability terms for specified polymers
        min_pol_term = alpha_min/vol_min
        max_pol_term = alpha_max/vol_max
        med_pol_term = alpha_med/vol_med

        multi_file.write('%f, %f, %f, %f, %f, %f, \n' % (min_dip_term, max_dip_term, med_dip_term, min_pol_term, max_pol_term, med_pol_term))

    # spear_file.write('%d, %f, %f, %f, \n' %
        # (gen_counter, spear_05, spear_10, spear_15))

    # write polymer population to file
    for polymer in population:
        poly_name = utils.make_file_name(polymer, poly_size)
        population_file.write('%s, ' % (poly_name))
    population_file.write('\n')


    if opt_property != 'dip_pol':
        for value in poly_property_list:
            values_file.write('%f, ' % (value))
        values_file.write('\n')
    else:
        for item in poly_property_list:
            for value in item:
                values_file.write('%f, ' % (value))
            values_file.write('-break-, ')
        values_file.write('\n')

    # keep track of number of successive generations meeting Spearman criterion
    if spear > 0.92:
        spear_counter += 1
    else:
        spear_counter = 0

    # keep track of number of successive generations meeting property value convergence criterion
    if opt_property != 'dip_pol':
        if max_test >= (max_init - max_init * 0.05) and max_test <= (max_init + max_init * 0.05):
            prop_value_counter += 1
        else:
            prop_value_counter = 0

    # close all output files
    analysis_file.close()
    population_file.close()
    values_file.close()
    if opt_property == 'dip':
        dip_polar_file.close()
    if opt_property == 'pol':
        polar_dip_file.close()
    if opt_property == 'dip_pol':
        multi_file.close()
    #spear_file.close()

    # make backup copies of output files
    shutil.copy('gens_analysis.txt', 'gens_analysis_copy.txt')
    shutil.copy('gens_population.txt', 'gens_population_copy.txt')
    shutil.copy('gens_values.txt', 'gens_values_copy.txt')
    if opt_property == 'dip':
        shutil.copy('gens_dip_polar.txt', 'gens_dip_polar_copy.txt')
    if opt_property == 'pol':
        shutil.copy('gens_polar_dip.txt', 'gens_polar_dip_copy.txt')
    if opt_property == 'dip_pol':
        shutil.copy('gens_multi.txt', 'gens_multi_copy.txt')
    #shutil.copy('gens_spear.txt', 'gens_spear_copy.txt')

    params = [pop_size, poly_size, num_mono_species, opt_property, smiles_list, sequence_list,
              mono_list, population, poly_property_list, n, gen_counter, spear_counter, prop_value_counter]
    return(params)


def main():
    # flag for restart from save
    restart = 'y'

    # number of polymers in population
    pop_size = 32
    # number of monomers per polymer
    poly_size = 6
    # number of species of monomers in each polymer
    num_mono_species = 2
    # property of interest (options: molecular weight 'mw', dipole moment 'dip', polarizability 'pol', dipole+polar 'dip_pol')
    opt_property = 'dip_pol'

    # check for convergence among top 30% (or top 8, whichever is larger) candidates between 5 generations
    perc = 0.1

    # Read in monomers from input file
    # read_file = open('../input_files/1235MonomerList.txt', 'r')
    read_file = open(
        '/ihome/ghutchison/dch45/Chem_GA_Project/input_files/1235MonomerList.txt', 'r')

    # create list of monomer SMILES strings
    # assumes input file has one monomer per line
    smiles_list = []
    for line in read_file:
        # grab first token from each line
        temp = line.split()[0]
        smiles_list.append(temp)
    read_file.close()

    if restart == 'y':
        # reload parameters and random state from restart file
        open_params = open('last_gen.p', 'rb')
        params = pickle.load(open_params)
        open_params.close()

        open_rand = open('randstate.p', 'rb')
        randstate = pickle.load(open_rand)
        random.setstate(randstate)
        open_rand.close()

    else:
        # run initial generation if NOT loading from restart file
        params = init_gen(pop_size, poly_size, num_mono_species, opt_property, perc, smiles_list)

    # set convergence counters
    spear_counter = params[11]
    prop_value_counter = params[12]

    # while spear_counter < 10 or prop_value_counter < 10:
    for x in range(100):
        # run next generation of GA
        params = next_gen(params)

        # pickle parameters needed for restart
        params_file = open('last_gen.p', 'wb')
        pickle.dump(params, params_file)
        params_file.close()

        # pickle random state for restart
        randstate = random.getstate()
        rand_file = open('randstate.p', 'wb')
        pickle.dump(randstate, rand_file)
        rand_file.close()

        # update convergence counters
        spear_counter = params[11]
        prop_value_counter = params[12]

    # remove unnecessary copies
    del_copies = subprocess.call('(rm -f *_copy.txt)', shell=True)


if __name__ == '__main__':
    main()
