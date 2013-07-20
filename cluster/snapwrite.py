'''

Execution: ./snapwrite.py FOLDER OUTPUT

FOLDER: folder containing all the input files, which are:

    header.txt position.txt velocity.txt id.txt masses.txt (always)
    energy.txt density.txt smoothing.txt (in case there is gas)

The columns of the positions and velocities files should be separated with
tabs, as well as the columns of the header. Among the input files, the one
corresponding to the header is the only one that accepts comments, which
should be preceded by a #.

OUTPUT: name of the file which will be the output snapshot.

'''

import numpy as np
import struct
import sys

n_gas = 0
output = sys.argv[1]

# To create a list containing the non-commented, non-empty lines
# of the input file for the header
def process_input(file_):
    h_file = []
    input_ = open(file_, 'r')
    for line in input_:
        if line.find("#") != -1:
            continue
        elif line.find("\n") == 0:
			continue
        else:
            h_file.append(line)
    return h_file

def read_header(folder):
    h_file = process_input(folder + "header.txt")
    h_data = []
    global n_gas
    n_gas = int(h_file[0].split('\t')[0])
    for j in h_file[0].split('\t')[0:6]: # n_part
        h_data.append(int(j))
    for j in h_file[1].split('\t')[0:6]: # mass
        h_data.append(float(j))
    h_data.append(float(h_file[2].split('\t')[0])) # time
    h_data.append(float(h_file[3].split('\t')[0])) # redshift
    h_data.append(int(h_file[4].split('\t')[0])) # flag_sfr
    h_data.append(int(h_file[5].split('\t')[0])) # flag_feedback
    for j in h_file[6].split('\t')[0:6]:
        h_data.append(int(j)) # n_part_total
    h_data.append(int(h_file[7].split('\t')[0])) # flag_coooling
    h_data.append(int(h_file[8].split('\t')[0])) # num_files
    h_data.append(float(h_file[9].split('\t')[0])) # box_size
    h_data.append(float(h_file[10].split('\t')[0])) # omega0
    h_data.append(float(h_file[11].split('\t')[0])) # omega_lambda
    h_data.append(float(h_file[12].split('\t')[0])) # hubble_param

    # blank, present in the header
    for i in np.arange(96):
        h_data.append('\0')
    s = struct.Struct(' iiiiii dddddd d d i i iiiiii i i dddd cccccccccccc\
                        cccccccccccccccccccccccccccccccccccccccccccccccccc\
                        cccccccccccccccccccccccccccccccccc')
    packed_data = s.pack(*h_data)
    return packed_data

def write_dummy(f, n_dummies):
    dummy = [256]
    s = struct.Struct('i')
    d = s.pack(*dummy)
    for i in np.arange(n_dummies):
        f.write(d)

def write_block(f, block_data, data_type):
    write_dummy(f, 1)
    fmt = data_type * len(block_data)
    f.write(struct.pack(fmt, *block_data))
    write_dummy(f, 1)

def write_snapshot(folder):

	# Erasing the input file before opening it
    open(output, 'w').close()
    f = file(output, 'a')
    header_data = read_header(folder)
    write_dummy(f, 1)
    f.write(header_data)
    write_dummy(f, 1)
    
    pos_file = np.fromfile(folder + "position.txt", dtype=float, count=-1,
                                                                 sep='\t')
    vel_file = np.fromfile(folder + "velocity.txt", dtype=float, count=-1,
                                                                 sep='\t')
    ID_file = np.fromfile(folder + "id.txt", dtype=int, count=-1, sep='\t')
    mass_file = np.fromfile(folder + "masses.txt", dtype=float, count=-1,
                                                                sep='\t')
    if(n_gas > 0):
        U_file = np.fromfile(folder + "energy.txt", dtype=float, count=-1,
                                                                 sep='\t')
        rho_file = np.fromfile(folder + "density.txt", dtype=float,
                                                count=-1, sep='\t')
    smoothing_file = np.fromfile(folder + "smoothing.txt", dtype=float,
                                                        count=-1, sep='\t')

	# writing
    write_block(f, pos_file, 'f')
    write_block(f, vel_file, 'f')
    write_block(f, ID_file, 'i')
    write_block(f, mass_file, 'f')
    if(n_gas > 0):
        write_block(f, U_file, 'f')
        write_block(f, rho_file, 'f')
        write_block(f, smoothing_file, 'f')
    f.close()
