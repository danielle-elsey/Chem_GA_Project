import random
import pybel as pb
from statistics import mean

'''
#given a SMILES string, calculate and return molecular weight
def find_mw(smiles_fragment):
    mw_counter = 0

    for character in smiles_fragment:
        if character.isupper():
            #add element's MW
            temp_element = pt.elements.isotope(character)
            mw_counter += temp_element.mass

    return mw_counter
'''

#given a polymer population, list of monomer mw, and number of monomer units per polymer, calculate polymer molecular weights
def find_poly_mw(population, mw_list, num_mono):
    poly_mw_list = []
    for polymer_index in range(len(population)):
        temp_mw = 0
        for monomer_index in range(num_mono):
            temp_mw += mw_list[population[polymer_index][monomer_index]]
        poly_mw_list.append(temp_mw)
    return poly_mw_list

#given list of polymer molecular weights and polymer population list, return list of heaviest 50% of polymers
def parent_select(polymer_mw_list, population):
    parent_count = int(len(polymer_mw_list)/2)
    parent_list = []
    print_mw_list = []

    for x in range (parent_count):
        temp_max_mw = max(polymer_mw_list)
        print_mw_list.append(temp_max_mw)
        #add heaviest polymer to parent_list
        parent_list.append(population[polymer_mw_list.index(temp_max_mw)])
        #replace heaviest mw with 0 in mw list (to preserve indexes)
        polymer_mw_list[polymer_mw_list.index(temp_max_mw)] = 0
    return parent_list

#given a polymer, number of monomers per polymer, and the number of types of monomers, creates random point mutation within polymer
def mutate(polymer, num_mono, num_type_mono):
    #create point for mutation
    point = random.randint(0, num_mono-1)
    #mutate polymer at specified monomer unit
    polymer[point] = random.randint(0, num_type_mono-1)
    #not sure if return is correct, as python seems to be pass by value (?)
    return polymer

#given parent population, number of monomers per polymer, and number of polymers in population,
#create equal number of children via crossover and return new population of parents and children
def crossover(parent_population, num_mono, num_poly):
    temp_p_list = parent_population
    #calculate number of children to generate (cannot assume same as number of parents, in case of odd number of total polymers)
    num_children = num_poly-len(temp_p_list)
    total_pop = temp_p_list

    for child in range(num_children):
        temp_p = temp_p_list[child]
        temp_c = []
        for monomer in range(num_mono):
            temp_c.append(temp_p[random.randint(0, num_mono-1)])
        total_pop.append(temp_c)
    return total_pop

#given a polymer as a list of monomer indicies from population, list of all monomers, and number of monomers per polymer, construct SMILES string of polymer
def construct_polymer_string(polymer, smiles_list, num_mono):
    poly_string = ''
    for monomer in range(num_mono):
        poly_string = poly_string + smiles_list[polymer[monomer]]
    return poly_string

def main():
    #number of polymers in population
    num_poly = 10
    #number of monomers per polymer
    num_mono = 5
    #number of types of monomers in each polymers
    num_type_mono = 2
    #threshold for minimum MW
    min_mw_std = 2500

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

    #create inital population of polymers
    population = []
    for x in range(num_poly):
        temp_poly = []
        for y in range(num_mono):
            temp_mono_index = random.randint(0,len(smiles_list)-1)
            temp_poly.append(temp_mono_index)
        tuple(temp_poly)
        population.append(temp_poly)

    #Calculate polymer molecular weights
    poly_mw_list = find_poly_mw(population, mw_list, num_mono)

    #print(poly_mw_list)

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
            polymer = mutate(polymer, num_mono, len(smiles_list))

        #Crossover - create children to repopulate bottom 50% of polymers in population
        population = crossover(population, num_mono, num_poly)


        poly_mw_list = find_poly_mw(population, mw_list, num_mono)

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
        polymer_smiles_list.append(construct_polymer_string(polymer, smiles_list, num_mono))
    print(polymer_smiles_list)
    '''

if __name__ == '__main__':
    main()
