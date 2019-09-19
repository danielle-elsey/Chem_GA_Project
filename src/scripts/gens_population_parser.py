def make_all_monos_list(file, generation_num):
    '''
    Make list of all monomers

    Parameters
    ---------
    file: str
        name of file to be parsed
    generation_num: int
        line number of first generation of consecutive pair to be analyzed

    Returns
    -------
    all_monos: list
        list containing two sublists of strings of monomer indexes used in each generation
    '''


    read_file = open(file, 'r')

    # read over header (line 0) and previous generations
    i = 0
    while i != generation_num:
        read_file.readline()
        i += 1

    all_monos = []

    for x in range(2):
        gen_list = read_file.readline().split()
        gen_monos = []
        for polymer in gen_list:
            temp = polymer.split('_')
            #add both monomer indexes from polymer string to list
            gen_monos.append(temp[0])
            gen_monos.append(temp[1])
        all_monos.append(gen_monos)

    read_file.close()

    return all_monos


def count_monos(all_monos):
    '''
    Count number of occurences of each monomer and write data to output file in the form:
    Run SMILES Count

    Parameters
    ---------
    all_monos: list
        list containing two sublists of strings of monomer indexes used in each generation

    Returns
    -------
    N/A
    '''


    write_file = open('run_test_mono_counts.txt', 'w+')
    # write headers for output file
    write_file.write('Run\tSMILES\tCount\n')

    # loop over each generation in list
    gen_counter = 0
    for generation in all_monos:
        gen_counter += 1
        counted_monos = []
        # loop over each monomer in given generation
        for mono in generation:
            count = generation.count('%s' % mono)
            # skip over types of monomers that have already been counted
            if mono not in counted_monos:
                write_file.write('%i\t%s\t%i\n' % (gen_counter, mono, count))


    write_file.close()


def main():

    all_monos = make_all_monos_list('gens_population.txt', 2)

    count_monos(all_monos)


if __name__ == '__main__':
    main()
