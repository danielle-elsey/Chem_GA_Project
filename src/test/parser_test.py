import glob
import pandas as pd

def main():
    # initialize list of all monomer data dictionaries
    data = []
    dataset = '/home/user/Dev/Mini_GA_Project/src/output/*.out'
    print('here1')
    for file in glob.iglob(dataset):
        print('here2')
        fname = file.split('/')[-1].split('.')[0]
        print(fname)

        # parse out property values
        read_file = open(file, 'r')
        for line in read_file:
            # create list of tokens in line
            tokens = line.split()
            if line.startswith(" Mol. α(0)"):
                polar = float(tokens[4])
            elif line.startswith("   full:"):
            # dipole tensor - STILL A LIST OF STRINGS (not floats)
            # TODO: add tensor functionality later
                dipole_line = tokens
                dipole = float(tokens[4])
                # break inner for loop to avoid overwriting with other lines starting with "full"
                break
        read_file.close()

        # create dictionary of data values for monomer
        d = {}
        d.update({'mono':fname})
        d.update({'mu':dipole})
        d.update({'alpha':polar})
        data.append(d)

    print('here3')
    # create df and csv of all monomer data
    df = pd.DataFrame(data, columns = ['mono', 'mu', 'alpha'])
    df.to_csv('mono.csv', index = False)

if __name__ == '__main__':
    main()
