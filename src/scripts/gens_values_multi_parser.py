import numpy as np
import pandas as pd



def main():

    read_file = open('gens_values.txt', 'r')

    dip_values = []
    pol_values = []

    # remove header
    next(read_file)
    # parse data
    for line in read_file:
        lists = line.split('-break-, ')
        dip_list = lists[0].split(', ')
        pol_list = lists[1].split(', ')

        del dip_list[-1]
        del pol_list[-1]

        for x in range(len(dip_list)):
            dip_list[x] = float(dip_list[x])
            pol_list[x] = float(pol_list[x])

        # sort values for each generation highest to lowest
        dip_list.sort(reverse=True)
        pol_list.sort(reverse=True)

        dip_values.append(dip_list)
        pol_values.append(pol_list)

    read_file.close()

    # create dataframes where each column is a generation
    dip_df = pd.DataFrame(np.array(dip_values)).transpose()
    pol_df = pd.DataFrame(np.array(pol_values)).transpose()

    dip_df.to_csv(r'dip_df.csv')
    pol_df.to_csv(r'pol_df.csv')




if __name__ == '__main__':
    main()
