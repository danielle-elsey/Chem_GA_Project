import string
from itertools import product

'''
#p = combinations_with_replacement('AB', 4)
p = product(range(2), repeat=4)
for i in list(p):
    print(i)

#test_dict = dict.fromkeys(range(10), string.ascii_uppercase[])
test_dict = dict(zip(range(1, 27), string.ascii_uppercase))

print(test_dict)
'''


#given the number of types of monomers in each polymer, find all possible sequences [combinations]
def find_sequences(num_type_mono):
    #find length of each sequence
    seq_len = num_type_mono**2

    #find all possible sequences as numerical lists [start index 0]
    #(use cartesian product of each type of monomer over sequence length)
    numer_seqs = product(range(num_type_mono), repeat=seq_len)


    #create dictionary mapping numbers [start index 0] to uppercase English alphabet letters
    alpha_dict = dict(zip(range(26), string.ascii_uppercase))
    #translate numerical sequences into alphabetic find_sequences
    alpha_seqs = []
    for sequence in numer_seqs:
        temp_seq = []
        for number in sequence:
            temp_seq.append(alpha_dict[number])
        alpha_seqs.append(temp_seq)

    return alpha_seqs

print(find_sequences(3))
