import numpy as np
import scipy.constants as sc

import utils

# rank polymers based on highest property value = best
def simple_descending(score_list):
    # make list of indicies of polymers in population, sorted based on property values
    ranked_indicies = list(np.argsort(score_list))
    # reverse list so highest property value = 0th
    ranked_indicies.reverse()

    return ranked_indicies

def comb_dip_pol(dip_prop_list, pol_prop_list, dip_coeff, pol_coeff, population, poly_size):

    # find polymer volumes
    vol_list = utils.make_vol_list(population, poly_size)

    score_list = []
    for x in range(len(dip_prop_list)):
        
        # alpha in bohr^3, mu in debye, vol in A^3
        alpha = pol_prop_list[x]
        mu = dip_prop_list[x]
        vol = vol_list[x]

        # convert polarizability (volume) to SI units (bohr^3 to m^3)    
        si_alpha_prime = alpha * (5.29177E-11)**3

        # convert dipole moment to SI units (debye to C*m)
        # 1 debye = 3.3E-30 C*m
        si_mu = mu * 3.3356E-30

        # convert volume to SI units (A^3 to m^3)
        si_vol = vol * 1E-30

        # define molecular density N/V (assume packing of 68%)
        n_v = 0.68/si_vol    

        # calculate Debye equation term values (assume temperature of 300K)
        a_term = n_v * (si_alpha_prime * 4 * sc.pi * sc.epsilon_0)/(3 * sc.epsilon_0)
        m_term = n_v * (si_mu**2)/(9  * sc.epsilon_0 * sc.Boltzmann * 300)

        # calculate RHS of Debye equation
        debye_sum = a_term + m_term 

        temp_score = (debye_sum)
        score_list.append(temp_score)

    ranked_indicies = simple_descending(score_list)

    return(ranked_indicies)
