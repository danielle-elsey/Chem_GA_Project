import scoring

def fitness_fn(opt_property, *poly_property_list):
    properties = []
    for list in poly_property_list:
        properties.append(list)

    if opt_property in ('mw', 'dip', 'pol'):
        ranked_indicies = scoring.simple_descending(properties[0])
    elif opt_property == 'dip_pol':
        ranked_indicies = scoring.comb_dip_pol(properties[0], properties[1], 1, 1)
    else:
        print('Error: opt_property not recognized; Traceback: fitnesss_fn')
    return(ranked_indicies)


def main():
    dip_prop_list = [5, 17, 99, 12, 11]
    pol_prop_list = [777, 234, 1009, 238, 561]
    dip_coeff = 1
    pol_coeff = 1

    test1 = scoring.comb_dip_pol(dip_prop_list, pol_prop_list, dip_coeff, pol_coeff)

    test2 = fitness_fn('dip_pol', dip_prop_list, pol_prop_list)

    print('test1 \n', test1)
    print('\ntest2 \n', test2)


if __name__ == '__main__':
    main()
