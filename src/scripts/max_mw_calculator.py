import pybel as pb
from copy import deepcopy

def find_poly_mw(polymer, poly_size, mw_list):
    temp_mw = 0
    # counter for cycling over sequence
    cycle = 0
    # loop over total monomers in each polymer
    for i in range(poly_size):
        #print("cycle:", cycle)
        # cycle over sequence until total number of monomers is reached
        seq_cycle_location = i-cycle*(len(polymer[0]))
        seq_monomer_index = polymer[0][seq_cycle_location]
        #print("sequence index:", i-cycle*(len(polymer[0])))
        if seq_cycle_location == (len(polymer[0])-1):
            cycle += 1
        # find index in polymer of monomer index (from smiles_list) from given sequence value
        monomer_index = polymer[seq_monomer_index+1]
        temp_mw += mw_list[monomer_index]
    return temp_mw

def main():
    #number of monomers per polymer
    poly_size = 6

    #Read in monomers from input file
    read_file = open('../input_files/1235MonomerList.txt','r')

    #create list of monomer SMILES strings
    #assumes input file has one monomer per line
    smiles_list = []
    for line in read_file:
        #grab first token from each line
        temp = line.split()[0]
        smiles_list.append(temp)

    read_file.close()

    # create list of mw's
    mw_list = []
    for monomer in smiles_list:
        # calculate MW for each monomer
        temp_mw= pb.readstring("smi", monomer).molwt
        mw_list.append(temp_mw)

    # create sorted list of mw's
    sort_mw_list = deepcopy(mw_list)
    sort_mw_list = sorted(sort_mw_list, reverse = True)

    # find highest possible polymer mw
    # test_high_mw = sort_mw_list[0]*poly_size
    temp_polymer = [(0,0,0,0), mw_list.index(sort_mw_list[0]), mw_list.index(sort_mw_list[0])]
    mw0 = find_poly_mw(temp_polymer, poly_size, mw_list)

    # find second highest polymer mw
    temp = [(0,0,0,1), mw_list.index(sort_mw_list[0]), mw_list.index(sort_mw_list[1])]
    mw1 = find_poly_mw(temp, poly_size, mw_list)

    # find third highest polymer mw
    temp = [(0,0,0,1), mw_list.index(sort_mw_list[1]), mw_list.index(sort_mw_list[0])]
    mw2 = find_poly_mw(temp, poly_size, mw_list)

    # find fourth highest polymer mw
    temp = [(0,0,0,0), mw_list.index(sort_mw_list[1]), mw_list.index(sort_mw_list[1])]
    mw3 = find_poly_mw(temp, poly_size, mw_list)

    print(mw0, mw1, mw2, mw3)

    '''
    max_mono_mw = max(mw_list)
    max_poly_mw = poly_size*max(mw_list)

    print(max(mw_list))

    test_polymer = [(0,0,0,0), mw_list.index(max_mono_mw), mw_list.index(max_mono_mw)]
    test_poly_mw = find_poly_mw(test_polymer, poly_size, mw_list)
    print("test_poly_mw", test_poly_mw)

    print("actual max poly mw", max_poly_mw)
    '''

if __name__ == '__main__':
    main()
