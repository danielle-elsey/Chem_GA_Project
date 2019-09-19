poly_dipole_list = []
poly_polar_list = []
num_succ_reads = 0

# parse output file for static polarizability and dipole moment, check for xTB success before recording values
read_output = open('output/%s.out' % (file_name), 'r')

for line in read_output:
    # create list of tokens in line
    tokens = line.split()

    if line.startswith(" Mol. α(0)"):
        temp_polar = float(tokens[4])
        num_succ_reads += 1
    elif line.startswith("molecular dipole:"):
        # iterate down to line with dipole value - necessary b/c dip value on non unique line
        next(read_output)
        next(read_output)
        new_tokens = next(read_output).split()

        # dipole tensor - STILL A LIST OF STRINGS (not floats)
        # TODO: add tensor functionality later
        temp_dipole = float(new_tokens[4])
        num_succ_reads += 1
    if line.startswith("#  - GEOMETRY OPTIMIZATION FAILED!"):
        num_succ_reads = 0

read_output.close()

if num_succ_reads == 2:
    poly_dipole_list.append(temp_dipole)
    poly_polar_list.append(temp_polar)
