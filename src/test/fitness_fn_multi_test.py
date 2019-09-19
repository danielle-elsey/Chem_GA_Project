
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

    print('avg_ranks', avg_ranks)

    # make list of indicies (0th = best) of polymers in population, sorted based on averaged ranks
    ranked_indicies = list(np.argsort(avg_ranks))
    # reverse list so highest property value = 0th
    #ranked_indicies.reverse()

    return ranked_indicies


def test():
    opt_property_1 = 'dip'
    # prop_list_1 = [8, 16.9, 16.52, 16.53]
    prop_list_1 = [39.215000, 38.663000, 38.213000, 35.150000, 41.555000, 38.409000, 39.771000, 33.835000, 34.705000, 40.828000, 36.014000, 42.352000, 35.445000, 48.179000, 44.298000, 35.424000, 34.401000, 30.154000, 29.354000, 14.965000, 62.926000, 30.127000, 34.179000, 24.825000, 39.266000, 12.539000, 34.188000, 34.828000, 18.548000, 19.052000, 33.447000, 31.689000]

    opt_property_2 = 'pol'
    # prop_list_2 = [600, 695, 888, 353]
    prop_list_2 = [1426.386009, 1394.098348, 1393.923648, 1513.365354, 1346.154773, 1366.956205, 1361.362124, 1513.368000, 1486.598373, 1328.821168, 1369.866550, 1266.239109, 1365.870890, 1239.820860, 1244.918252, 1361.403948, 1128.732168, 1190.357489, 1600.275502, 1861.170007, 1193.336960, 1296.401774, 1338.207891, 1359.166817, 1393.884111, 1227.588191, 1456.472057, 1120.638671, 964.288807, 875.430378, 1392.751754, 1319.436531]


    population = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'AA', 'BB', 'CC', 'DD', 'EE', 'FF']

    output = fitness_fn_multi(opt_property_1, prop_list_1, opt_property_2, prop_list_2)

    print(output)


    parent_list = []
    for x in range(32):
        parent_list.append(population[output[x]])

    print(parent_list)




def main():
    test()




if __name__ == '__main__':
    main()
