import random
import pybel as pb
from statistics import mean


#given a polymer population, list of monomer mw, and number of monomer units per polymer, calculate polymer molecular weights
def find_poly_mw(population, mw_list, poly_size):
    poly_mw_list = []
    for polymer_index in range(len(population)):
        temp_mw = 0
        for monomer_index in range(poly_size):
            temp_mw += mw_list[population[polymer_index][monomer_index]]
        poly_mw_list.append(temp_mw)
    return poly_mw_list

#given list of polymer molecular weights and polymer population list, return list of heaviest 50% of polymers
def parent_select(polymer_mw_list, population):
    parent_count = int(len(polymer_mw_list)/2)
    parent_list = []

    for x in range (parent_count):
        temp_max_mw = max(polymer_mw_list)
        #add heaviest polymer to parent_list
        parent_list.append(population[polymer_mw_list.index(temp_max_mw)])
        #replace heaviest mw with 0 in mw list (to preserve indexes)
        polymer_mw_list[polymer_mw_list.index(temp_max_mw)] = 0
    return parent_list

#given a polymer, and the number of monomers per polymer, creates random point mutation within polymer
def mutate(polymer, poly_size):
    #create point for mutation
    point = random.randint(0, poly_size-1)
    #mutate: change monomer at specified point to another monomer in same polymer
    #ensures mutant monomer is of an allowed type for given polymer
    mutant_mono = random.randint(0, poly_size-1)
    #ensure mutant monomer is not from mutation point
    while mutant_mono == point:
        mutant_mono = random.randint(0, poly_size-1)
    polymer[point] = polymer[mutant_mono]
    #not sure if return is correct, as python seems to be pass by value (?)
    return polymer

#given parent population, number of monomers per polymer, and number of polymers in population,
#create equal number of children via crossover and return new population of parents and children
def crossover(parent_list, poly_size, pop_size):
    #calculate number of children to generate
    num_children = pop_size-len(parent_list)
    #initialize new population with parents
    new_pop = parent_list

    for child in range(num_children):
        #randomly select two parents (as indexes from parent list) to cross
        parent_a = random.randint(0, len(parent_list)-1)
        parent_b = random.randint(0, len(parent_list)-1)

        #ensure parents are indiviudals
        while parent_b == parent_a:
            parent_b = random.randint(0, len(parent_list)-1)

        #number of monomers taken from parent A
        num_mono_a = int(poly_size/2)
        #number of monomers taken from parent B
        num_mono_b = poly_size-num_mono_a

        #create child with first half monomers from A, second half from B
        temp_child = []
        for monomer in range(num_mono_a):
            temp_child.append(parent_list[parent_a][monomer])
        for monomer in range(num_mono_b):
            temp_child.append(parent_list[parent_b][num_mono_a-1+monomer])

        new_pop.append(temp_child)

    return new_pop

#given a polymer as a list of monomer indicies from population, list of all monomers, and number of monomers per polymer, construct SMILES string of polymer
def construct_polymer_string(polymer, smiles_list, poly_size):
    poly_string = ''
    for monomer in range(poly_size):
        poly_string = poly_string + smiles_list[polymer[monomer]]
    return poly_string

def main():
    #number of polymers in population
    pop_size = 10
    #number of monomers per polymer
    poly_size = 5
    #number of types of monomers in each polymers
    num_type_mono = 2
    #threshold for minimum MW
    min_mw_std = 1500

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

    #create inital population as list of polymers
    population = []
    for polymer in range(pop_size):
        #select monomer types for polymer
        mono_types = []
        for num in range(num_type_mono):
            temp_mono_type = random.randint(0, len(smiles_list)-1)
            mono_types.append(temp_mono_type)

        #create polymer as a list of monomer indexes from smiles_list
        temp_poly = []
        for monomer in range(poly_size):
            temp_mono_index = random.randint(0,len(mono_types)-1)
            temp_poly.append(mono_types[temp_mono_index])
        tuple(temp_poly)

        #add polymer to population
        population.append(temp_poly)

    #Calculate polymer molecular weights
    poly_mw_list = find_poly_mw(population, mw_list, poly_size)

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
            polymer = mutate(polymer, poly_size)

        #Crossover - create children to repopulate bottom 50% of polymers in population
        population = crossover(population, poly_size, pop_size)


        poly_mw_list = find_poly_mw(population, mw_list, poly_size)

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
