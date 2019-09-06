import numpy as np
from scipy.constants import Boltzmann

# rank polymers based on highest property value = best
def simple_descending(score_list):

    # make list of indicies of polymers in population, sorted based on property values
    ranked_indicies = list(np.argsort(score_list))
    # reverse list so highest property value = 0th
    ranked_indicies.reverse()

    return ranked_indicies

def comb_dip_pol(dip_prop_list, pol_prop_list, dip_coeff, pol_coeff):
    score_list = []
    for x in range(len(dip_prop_list)):
        alpha = pol_prop_list[x]
        mu_2 = dip_prop_list[x]**2
        vol = 1

        temp_score = ((alpha/vol)**pol_coeff + (mu_2/(3*Boltzmann*298*vol))**dip_coeff)

        score_list.append(temp_score)

    ranked_indicies = simple_descending(score_list)

    return(ranked_indicies)
