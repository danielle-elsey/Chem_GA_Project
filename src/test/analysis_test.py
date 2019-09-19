import pandas as pd
import numpy as np

def analysis(data_file, perc):

    df = pd.read_csv(data_file, sep='\t', header=0)

    # Makes a pivot table with SMILES as the rows and columns of different row numbers
    pt = pd.pivot_table(df, index=["SMILES"], values=["Run"], columns=["Run"], fill_value=0)
    # Add across each row to get the total number of monomers and sort from largest to smallest
    pt_sum = pt.sum(1).sort_values(ascending=False)

    # Takes the values that are larger than the quartile being analyzed
    a = pt_sum[pt_sum > pt_sum.quantile(perc)]

    # Identifies the number of monomers within that top percentage which will be used in the analysis
    num_monomers = len(a.index)

    # Generates smaller pivot table which includes only the top monomers
    pt_subset = pt.loc[a.index.values]
    print(pt_subset)

    # Pairwise spearman correlation of all runs
    corr = pt_subset.corr(method='spearman')

    # Removes duplicate Spearman correlations and substituting with nan
    c = corr.copy()
    c.values[np.tril_indices_from(c)] = np.nan

    s = c.unstack().mean()

    # Generate ndarray of values of the Spearman correlations
    s1 = s.values

    # Removes nan values leaving just the desired  Spearman correlations
    sv = s1[~np.isnan(s1)]

    # Converts numpy array to list of values which can then be saved individually
    sv_values = list(sv)

    # Saves the values in a test file in the format percentage, Generations, Spearman Value which
    # can then be plotted. This file will just be appended to which allows many data sets to be combined into
    # one long list with all data which can then be plotted, analyzed, etc.

    for i in sv_values:
        #generations = int(sys.argv[3])
        generations = 2

        #monomers = int(sys.argv[4])
        monomers = 1234

        #with open(sys.argv[2], "a") as output_file:
        with open('spear_output.txt', 'a') as output_file:
            output_file.write("%i\t%i\t%s\t%f\n" % (monomers, generations, perc, i))

    return sv_values[0]

def main():
    analysis('spear_input.txt', 0.7)



if __name__ == '__main__':
    main()
