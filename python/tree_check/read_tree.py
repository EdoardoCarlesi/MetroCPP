'''
        This script reads a merger tree from an ASCII file and returns an array containing the number of particles per step
'''

import numpy as np

'''
        file_ahf = open(file_name, 'r')

        line = file_ahf.readline()
        halos_ahf = []
        count = 0

        while line:
                full_line = file_ahf.readline()
                line = full_line.strip()
                column = line.split()
                n_col = len(column)

                if n_col > 1 and line[0] != "#":
'''

def read_ahf_tree(f_name, n_steps):
    n_steps = int(n_steps) -2
    ahf_tree = np.zeros((n_steps))

    f_ahf = open(f_name, 'r')
    line = f_ahf.readline()
    i_line = 0

    while line:
        this_line = f_ahf.readline()
        line = this_line.strip()
        cols = line.split()
        n_col = len(cols)

        # Particle number is in column number 6
        if n_col > 1 and line[0] != "#":
            ahf_tree[i_line] = float(cols[5])
            i_line += 1

    return ahf_tree



def read_mcpp_tree(f_name, n_steps):
    n_steps = int(n_steps)-1
    mcpp_tree = np.zeros((n_steps))
    f_mcpp = open(f_name, 'r')

    line = f_mcpp.readline()
    i_line = 0

    while line:
        this_line = f_mcpp.readline()
        line = this_line.strip()
        cols = line.split()
        n_col = len(cols)

        # Particle number is in column number 5
        if n_col > 1 and line[0] != "#":
            mcpp_tree[i_line] = float(cols[1])

            i_line += 1


    return mcpp_tree



def read_darf_list(f_dwarfs):

	dwarf_list = []

	return dwarf_list
