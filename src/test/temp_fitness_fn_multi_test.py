
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
    # reverse list so highest property value = 0th
    ranked_indicies.reverse()

    return ranked_indicies


def parent_select(opt_property, population, poly_property_list):
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
    if opt_property != 'dip_pol':
        fitness_list = fitness_fn(opt_property, poly_property_list)
    else:
        fitness_list = fitness_fn_multi('dip', poly_property_list[0], 'pol', poly_property_list[1])

    # make list of top polymers
    parent_list = []
    for x in range(parent_count):
        parent_list.append(population[fitness_list[x]])

    return parent_list


def main():

    dip_list = [1.282000, 2.445000, 2.554000, 1.326000, 3.594000, 3.800000, 3.845000, 5.754000, 1.638000, 2.935000, 6.049000, 2.450000, 7.302000, 2.737000, 4.913000, 5.009000, 7.194000, 30.759000, 22.975000, 11.034000, 10.859000, 15.442000, 4.585000, 18.540000, 4.307000, 7.481000, 13.107000, 14.059000, 3.075000, -10.000000, -10.000000, 4.298000, 3.202000]
    pol_list = [245.512110, 298.087766, 319.914303, 406.160803, 313.726667, 294.072621, 297.983325, 235.643561, 485.655179, 370.380665, 277.839281, 548.130880, 245.566080, 526.476321, 357.408870, 360.620156, 560.257659, 1026.243772, 285.451194, 790.346375, 402.264543, 638.055904, 419.925084, 882.302953, 485.605424, 255.568479, 581.089518, 275.436257, 610.444738, -10.000000, -10.000000, 909.161144, 275.127103]

    mw_list = [3594.984920, 3004.660540, 2741.678480, 2473.458920, 2062.495680, 1970.493180, 1770.500740, 1681.095860, 1641.411660, 1615.612000, 1599.259840, 1587.394338, 1584.219080, 1572.292960, 1539.045880, 1487.074280, 2593.878000, 1592.110080, 1295.466080, 1473.304180, 2733.745600, 1774.014438, 1528.089620, 903.125720, 1403.431928, 1299.593480, 1800.536060, 1271.490960, 4015.583360, 1191.715840, 1397.941600, 1571.468560,
]

# ***Check on len = 33 problem, then test parent_select/fitness_fn_multi with -10s***
    print(len(dip_list), len(pol_list), len(mw_list))

    # fitness_list = fitness_fn_multi('dip', dip_list, 'pol', pol_list)
    # print('fitness_list', fitness_list)
    #
    # parent_list = parent_select('dip_pol', )

if __name__ == '__main__':
    main()
