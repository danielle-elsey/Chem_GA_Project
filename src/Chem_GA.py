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


def fitness_fn(opt_property, poly_property_list, population, poly_size):

    if opt_property in ('mw', 'dip', 'pol'):
        ranked_indicies = scoring.simple_descending(poly_property_list)
    elif opt_property == 'dip_pol':
        ranked_indicies = scoring.comb_dip_pol(poly_property_list[0], poly_property_list[1], 1, 1, population, poly_size)
    else:
        print('Error: opt_property not recognized; Traceback: fitnesss_fn')
    return(ranked_indicies)


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

    fitness_list = fitness_fn(opt_property, poly_property_list, population, poly_size)

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
    Initializes parameters, creates population, and runs initial generation

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


    # create new analysis output files
    quick_file = open('quick_analysis_data.txt', 'w+')
    analysis_file = open('full_analysis_data.txt', 'w+')  
    
    # write analysis file headers
    if opt_property == 'mw':
        quick_file.write('gen,min_mw,max_mw,avg_mw,spearman\n')
        analysis_file.write('run,gen,filename,mw\n')
    elif opt_property == 'pol':
        quick_file.write('gen,min_alpha,max_alpha,avg_alpha,spearman\n')
        analysis_file.write('run,gen,filename,POL,dip\n')
    elif opt_property == 'dip':
        quick_file.write('gen,min_mu,max_mu,avg_mu,spearman\n')
        analysis_file.write('run,gen,filename,pol,DIP\n')
    elif opt_property == 'dip_pol':
        quick_file.write('gen,min_alpha,max_alpha,med_alpha,min_mu,max_mu,med_mu,spearman\n')
        analysis_file.write('run,gen,filename,alpha,mu,vol\n')

    # find values for quick analysis file
    if opt_property != 'dip_pol':
        min_test = min(poly_property_list)
        max_test = max(poly_property_list)
        avg_test = mean(poly_property_list)
    else:
        fitness_list = fitness_fn('dip_pol', poly_property_list, population, poly_size)

        min_polymer = population[fitness_list[len(fitness_list)-1]]
        max_polymer = population[fitness_list[0]]
        median = int((len(fitness_list)-1)/2)
        med_polymer = population[fitness_list[median]]

        min_test_a = poly_property_list[1][population.index(min_polymer)]
        max_test_a = poly_property_list[1][population.index(max_polymer)]
        med_test_a = poly_property_list[1][population.index(med_polymer)]

        min_test_m = poly_property_list[0][population.index(min_polymer)]
        max_test_m = poly_property_list[0][population.index(max_polymer)]
        med_test_m = poly_property_list[0][population.index(med_polymer)]


    # write to quick analysis file
    if opt_property != 'dip_pol':
        quick_file.write('1,%f,%f,%f,n/a\n' % (min_test, max_test, avg_test))
    else:
        quick_file.write('1,%f,%f,%f,%f,%f,%f,n/a\n' % (min_test_a, max_test_a, med_test_a, min_test_m, max_test_m, med_test_m))


    # find values and write to full analysis file
    fitness_list = fitness_fn(opt_property, poly_property_list, population, poly_size)
    # loop over every polymer in population
    for x in range(pop_size):
        #use ranks so that polymers appear in order in output file
        poly = population[fitness_list[x]]

        file_name = utils.make_file_name(poly, poly_size)
        prop = poly_property_list[population.index(poly)]

        if opt_property == 'mw':
            analysis_file.write('1,1,%s,%f\n' % (file_name, prop))

        elif opt_property == 'pol':
            dip_val = dip_list[population.index(poly)]
            analysis_file.write('1,1,%s,%f,%f\n' % (file_name, prop, dip_val))
        
        elif opt_property == 'dip':
            pol_val = polar_list[population.index(poly)]
            analysis_file.write('1,1,%s,%f,%f\n' % (file_name, pol_val, prop))
        
        elif opt_property == 'dip_pol':
            alpha = poly_property_list[1][population.index(poly)]
            mu = poly_property_list[0][population.index(poly)]
            vol = utils.find_volume('opt/%s_opt.xyz' % (file_name))
            analysis_file.write('1,1,%s,%f,%f,%f\n' % (file_name, alpha, mu, vol))

    # close output files
    quick_file.close()
    analysis_file.close()

    # make backup copies of output files
    shutil.copy('quick_analysis_data.txt', 'quick_analysis_data_copy.txt')
    shutil.copy('full_analysis_data.txt', 'full_analysis_data_copy.txt')
    
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

    # create sorted monomer list with most freq first
    gen2 = utils.sort_mono_indicies_list(mono_list)

    # calculate Spearman correlation coefficient for begin and end sorted monomer lists
    spear = stats.spearmanr(gen1[:n], gen2[:n])[0]


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

    
    # # keep track of number of successive generations meeting Spearman criterion
    # if spear > 0.92:
    #     spear_counter += 1
    # else:
    #     spear_counter = 0

    # # keep track of number of successive generations meeting property value convergence criterion
    # if opt_property != 'dip_pol':
    #     if max_test >= (max_init - max_init * 0.05) and max_test <= (max_init + max_init * 0.05):
    #         prop_value_counter += 1
    #     else:
    #         prop_value_counter = 0


    # open analysis output files
    quick_file = open('quick_analysis_data.txt', 'a+')
    analysis_file = open('full_analysis_data.txt', 'a+')  
    
    # find values for quick analysis file
    if opt_property != 'dip_pol':
        min_test = min(poly_property_list)
        max_test = max(poly_property_list)
        avg_test = mean(poly_property_list)
    else:
        fitness_list = fitness_fn('dip_pol', poly_property_list, population, poly_size)

        min_polymer = population[fitness_list[len(fitness_list)-1]]
        max_polymer = population[fitness_list[0]]
        median = int((len(fitness_list)-1)/2)
        med_polymer = population[fitness_list[median]]

        min_test_a = poly_property_list[1][population.index(min_polymer)]
        max_test_a = poly_property_list[1][population.index(max_polymer)]
        med_test_a = poly_property_list[1][population.index(med_polymer)]

        min_test_m = poly_property_list[0][population.index(min_polymer)]
        max_test_m = poly_property_list[0][population.index(max_polymer)]
        med_test_m = poly_property_list[0][population.index(med_polymer)]


    # write to quick analysis file
    if opt_property != 'dip_pol':
        quick_file.write('%d,%f,%f,%f,%f\n' % (gen_counter, min_test, max_test, avg_test, spear))
    else:
        quick_file.write('%d,%f,%f,%f,%f,%f,%f,%f\n' % (gen_counter, min_test_a, max_test_a, med_test_a, min_test_m, max_test_m, med_test_m, spear))

    
    # find values and write to full analysis file
    fitness_list = fitness_fn(opt_property, poly_property_list, population, poly_size)
    # loop over every polymer in population
    for x in range(pop_size):
        #use ranks so that polymers appear in order in output file
        poly = population[fitness_list[x]]

        file_name = utils.make_file_name(poly, poly_size)
        prop = poly_property_list[population.index(poly)]

        if opt_property == 'mw':
            analysis_file.write('1,%d,%s,%f\n' % (gen_counter, file_name, prop))

        elif opt_property == 'pol':
            dip_val = dip_list[population.index(poly)]
            analysis_file.write('1,%d,%s,%f,%f\n' % (gen_counter, file_name, prop, dip_val))
        
        elif opt_property == 'dip':
            pol_value = polar_list[population.index(poly)]
            analysis_file.write('1,%d,%s,%f,%f\n' % (gen_counter, file_name, pol_value, prop))
        
        elif opt_property == 'dip_pol':
            alpha = poly_property_list[1][population.index(poly)]
            mu = poly_property_list[0][population.index(poly)]
            vol = utils.find_volume('opt/%s_opt.xyz' % (file_name))
            analysis_file.write('1,%d,%s,%f,%f,%f\n' % (gen_counter, file_name, alpha, mu, vol))

    
    # close output files
    quick_file.close()
    analysis_file.close()

    # make backup copies of output files
    shutil.copy('quick_analysis_data.txt', 'quick_analysis_data_copy.txt')
    shutil.copy('full_analysis_data.txt', 'full_analysis_data_copy.txt')

    params = [pop_size, poly_size, num_mono_species, opt_property, smiles_list, sequence_list,
              mono_list, population, poly_property_list, n, gen_counter, spear_counter, prop_value_counter]
    return(params)


def main():
    # flag for restart from save
    restart = 'n'

    # number of polymers in population
    pop_size = 2
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
    # spear_counter = params[11]
    # prop_value_counter = params[12]

    # while spear_counter < 10 or prop_value_counter < 10:
    for x in range(2):
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
        # spear_counter = params[11]
        # prop_value_counter = params[12]

    
    # delete final newline characters from files
    quick_file = open('quick_analysis_data.txt', 'a+')
    analysis_file = open('full_analysis_data.txt', 'a+')

    with open('quick_analysis_data.txt', 'rb+') as filehandle:
        filehandle.seek(-1, os.SEEK_END)
        filehandle.truncate()

    with open('full_analysis_data.txt', 'rb+') as filehandle:
        filehandle.seek(-1, os.SEEK_END)
        filehandle.truncate()

    # remove unnecessary copies
    del_copies = subprocess.call('(rm -f *_copy.txt)', shell=True)


if __name__ == '__main__':
    main()
