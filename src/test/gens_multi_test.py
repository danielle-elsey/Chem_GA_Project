
def multi_snippet():
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

