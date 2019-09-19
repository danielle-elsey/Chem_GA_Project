import numpy as np
from scipy import stats
from statistics import mean


def fitness_fn_multi(opt_property_1, prop_list_1, opt_property_2, prop_list_2):
    if opt_property_1 in ('mw', 'dip', 'pol') :
        # find ranks of properties
        get_ranks = list(stats.rankdata(prop_list_1, method = 'ordinal'))
        # convert ranks so that 1st = highest value
        ranks_1 = []
        for x in range (len(get_ranks)):
            ranks_1.append(len(get_ranks)-get_ranks[x])
    else:
        print("Error: opt_property not recognized. trace:fitness_fn")

    # make sorted list of polymer indicies based on second property
    if opt_property_2 in ('mw', 'dip', 'pol') :
        # find ranks of properties
        get_ranks = list(stats.rankdata(prop_list_2, method = 'ordinal'))
        # convert ranks so that 1st = highest value
        ranks_2 = []
        for x in range (len(get_ranks)):
            ranks_2.append(len(get_ranks)-get_ranks[x])
    else:
        print("Error: opt_property not recognized. trace:fitness_fn")

    # average ranks from both properties for each polymer
    avg_ranks= []
    for x in range(len(ranks_1)):
        temp_index = mean([float(ranks_1[x]), float(ranks_2[x])])
        avg_ranks.append(temp_index)

    # make list of indicies of polymers in population, sorted based on averaged ranks
    ranked_indicies = list(np.argsort(avg_ranks))
    
    return ranked_indicies


def main():
    population = ['A', 'B', 'C', 'D', 'E']

    dip_list = [5, 12, 10, 9, 19]
    pol_list = [501, 2009, 1012, 334, 987]

    poly_property_list = []
    poly_property_list.append(dip_list)
    poly_property_list.append(pol_list)

    fitness_list = fitness_fn_multi('dip', poly_property_list[0], 'pol', poly_property_list[1])

    ordered_list = []
    for x in range(5):
        ordered_list.append(population[fitness_list[x]])
    print('ordered_list', ordered_list)

    max_polymer = population[fitness_list[0]]
    min_polymer = population[fitness_list[len(fitness_list)-1]]
    median = int((len(fitness_list)-1)/2)
    med_polymer = population[fitness_list[median]]
    print('max_polymer', max_polymer)
    print('min_polymer', min_polymer)
    print('median', median)
    print('med_polymer', med_polymer)

    max_test = [poly_property_list[0][population.index(max_polymer)], poly_property_list[1][population.index(max_polymer)]]
    min_test = [poly_property_list[0][population.index(min_polymer)], poly_property_list[1][population.index(min_polymer)]]
    med_test = [poly_property_list[0][population.index(med_polymer)], poly_property_list[1][population.index(med_polymer)]]
    print('max_test', max_test)
    print('min_test', min_test)
    print('med_test', med_test)


if __name__ == '__main__':
    main()
