'''
Notes:
-each polymer is defined in the format [(numerical sequence), monomer1 index, ... monomerX index]
-population is list of these polymers
'''

import string
import random
import pybel as pb
from statistics import mean
from itertools import product
from copy import deepcopy

#given the number of types of monomers in each polymer, find all possible sequences [combinations]
def find_sequences(num_type_mono):
    #find length of each sequence
    seq_len = num_type_mono**2

    #find all possible sequences as numerical lists [start index 0]
    #(use cartesian product of each type of monomer over sequence length)
    numer_seqs = list(product(range(num_type_mono), repeat=seq_len))

    return numer_seqs

#given a population, polymer size, and list of monomer mw's, calculate polymer mw's
def find_poly_mw(population, poly_size, mw_list):
    poly_mw_list = []
    for polymer in population:
        temp_mw = 0
        #counter for cycling over sequence
        cycle = 0
        #loop over total monomers in each polymer
        for i in range(poly_size):
            #cycle over sequence until total number of monomers is reached
            seq_cycle_location = i-cycle*(len(polymer[0]))
            seq_monomer_index = polymer[0][seq_cycle_location]
            if seq_cycle_location == (len(polymer[0])-1):
                cycle += 1
            #find index in polymer of monomer index (from smiles_list) from given sequence value
            monomer_index = polymer[seq_monomer_index+1]
            temp_mw += mw_list[monomer_index]
        poly_mw_list.append(temp_mw)
    return poly_mw_list

#given list of sequences, list of molecular weights,
#and a population with polymers in index format: [sequence index, monomer1 index, ..., monomerX index]
#return population with polymers in human readable format: [[alpha sequence], 'monomer SMILES string 1', ..., 'monomer SMILES string X']
def make_pop_readable(population, sequence_list, smiles_list):

    #create dictionary mapping numbers [start index 0] to uppercase English alphabet letters
    alpha_dict = dict(zip(range(26), string.ascii_uppercase))

    readable_poly = []
    for polymer in population:
        temp_poly = []

        #find numerical sequence from index of sequence_list
        temp_numer_seq = sequence_list[polymer[0]]
        #translate numerical sequence into alphabetic sequence
        temp_alpha_seq = []
        for number in temp_numer_seq:
            temp_alpha_seq.append(alph_dict[number])
        temp_poly.append(temp_alpha_seq)

        for monomer_index in polymer[1:]:
            #find monomer SMILES string from index of smiles_list
            temp_mono_smiles = smiles_list[monomer_index]
            temp_poly.append(temp_mono_smiles)

        readable_poly.append(temp_poly)

    return readable_poly

#given list of polymer molecular weights and polymer population list, return list of heaviest 50% of polymers
def parent_select(polymer_mw_list, population):
    parent_count = int(len(polymer_mw_list)/2)
    parent_list = []
    temp_mw_list = deepcopy(polymer_mw_list)

    for x in range (parent_count):
        temp_max_mw = max(temp_mw_list)
        #add heaviest polymer to parent_list
        parent_list.append(population[polymer_mw_list.index(temp_max_mw)])
        #replace heaviest mw with 0 in mw list (to preserve indexes)
        temp_mw_list[temp_mw_list.index(temp_max_mw)] = 0
    return parent_list

#given a polymer and sequence and smiles lists, creates random point mutation within polymer
def mutate(polymer, sequence_list, smiles_list):
    #choose point of mutation (sequence or specific monomer)
    point = random.randint(0,len(polymer)-1)
    #replace sequence
    if point == 0:
        polymer[point] = sequence_list[random.randint(0,len(sequence_list)-1)]
    #or replace specific monomer
    else:
        polymer[point] = random.randint(0, len(smiles_list)-1)
    return polymer

#given parent population, number of monomers per polymer, and number of polymers in population,
#create equal number of children via crossover and return new population of parents and children
def crossover(parent_list, pop_size):
    #calculate number of children to generate
    num_children = pop_size-len(parent_list)
    #initialize new population with parents
    new_pop = deepcopy(parent_list)

    for child in range(num_children):
        #randomly select two parents (as indexes from parent list) to cross
        parent_a = random.randint(0, len(parent_list)-1)
        parent_b = random.randint(0, len(parent_list)-1)

        #ensure parents are unique indiviudals
        while parent_b == parent_a:
            parent_b = random.randint(0, len(parent_list)-1)

        #create child with sequence from A, monomer types from B
        temp_child = []
        temp_child.append(parent_list[parent_a][0])
        for monomer in parent_list[parent_b][1:]:
            temp_child.append(monomer)

        new_pop.append(temp_child)

    return new_pop


'''
#given a polymer as a list of monomer indicies from population, list of all monomers, and number of monomers per polymer, construct SMILES string of polymer
def construct_polymer_string(polymer, smiles_list, poly_size):
    poly_string = ''
    for monomer in range(poly_size):
        poly_string = poly_string + smiles_list[polymer[monomer]]
    return poly_string
'''

def main():
    #number of polymers in population
    pop_size = 10
    #number of monomers per polymer
    poly_size = 5
    #number of types of monomers in each polymer
    num_type_mono = 2
    #threshold for minimum MW
    min_mw_std = 3434

    #Read in monomers from input file
    read_file = open('../input_files/1235MonomerList.txt','r')

    #create list of monomer SMILES strings
    #assumes input file has one monomer per line
    smiles_list = []
    for line in read_file:
        #grab first token from each line
        temp = line.split()[0]
        smiles_list.append(temp)

    read_file.close()

    #calculate MW for each monomer
    mw_list = []
    for monomer in smiles_list:
        temp_mono = pb.readstring("smi", monomer).molwt
        mw_list.append(temp_mono)

    print("max", max(mw_list))
    #create all possible numerical sequences for given number of monomer types
    sequence_list = find_sequences(num_type_mono)

    #create inital population as list of polymers
    population = []
    for polymer in range(pop_size):
        temp_poly=[]

        #select sequence type for polymer
        poly_seq = sequence_list[random.randint(0,len(sequence_list)-1)]
        temp_poly.append(poly_seq)

        #select monomer types for polymer
        for num in range(num_type_mono):
            poly_monomer = random.randint(0, len(smiles_list)-1)
            temp_poly.append(poly_monomer)

        #add polymer to population
        population.append(temp_poly)

    #Calculate polymer molecular weights
    poly_mw_list = find_poly_mw(population, poly_size, mw_list)

    #set initial values for min, max, and avg polymer weights
    min_mw_test = min(poly_mw_list)
    max_mw_test = max(poly_mw_list)
    avg_mw_test = (max_mw_test-min_mw_test)/2.0

    #create new output read_file
    write_file = open("ga_polymer_output.txt","w+")
    write_file.write("min_wt avg_wt max_wt\n")

    #Loop
    while (min_mw_test < min_mw_std):

        #Selection - select best 50% as parents

        #select heaviest (best) 50% of polymers as parents
        population = parent_select(poly_mw_list, population)

        #Mutation - one point mutation per parent
        for polymer in population:
            polymer = mutate(polymer, sequence_list, smiles_list)

        #Crossover - create children to repopulate bottom 50% of polymers in population
        population = crossover(population, pop_size)


        poly_mw_list = find_poly_mw(population, poly_size, mw_list)

        #print(poly_mw_list)

        #find minimum polymer weight
        min_mw_test = min(poly_mw_list)
        max_mw_test = max(poly_mw_list)
        avg_mw_test = mean(poly_mw_list)

        write_file.write("{} {} {}\n".format(min_mw_test, avg_mw_test, max_mw_test))

        print(min_mw_test)


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
