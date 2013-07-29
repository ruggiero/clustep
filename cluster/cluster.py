'''
Use: python cluster.py SNAPSHOT, where SNAPSHOT is the name of the output file

Two files should be present in the script's execution folder: cluster_param.txt
and header.txt (see examples)
'''

from snapwrite import *
import rejection as rej # source found in source/rejection.pyx
import numpy as np
import numpy.random as nprand
import os
import sys

# Extracting the values of the global variables from cluster_param.txt
def init():
    global Mh, a, N
    if not (os.path.isfile("header.txt") and os.path.isfile(
            "cluster_param.txt")):
        print ("The following parameter files are required: header.txt and\n"
               "cluster_param.txt. Please check.")
        sys.exit(0)
    vars_ = process_input("cluster_param.txt") # From snapwrite.py
    Mh, a, N = float(vars_[0][0]), float(vars_[1][0]), float(vars_[2][0])

def inverse_cumulative(Mc):
    return (a * ((Mc * Mh)**0.5 + Mc)) / (Mh - Mc)

def set_positions():
    
    # This factor Mh 200^2 / 201^2 is used for restricting the radius to 200a
    radii = inverse_cumulative(nprand.sample(N) * ((Mh * 40000) / 40401))
    thetas = np.arccos(nprand.sample(N) * 2 - 1)
    phis = 2 * np.pi * nprand.sample(N)
    xs = radii * np.sin(thetas) * np.cos(phis)
    ys = radii * np.sin(thetas) * np.sin(phis)
    zs = radii * np.cos(thetas)
    coords = np.column_stack((xs, ys, zs))

    # Older NumPy versions freak out without this line
    coords = np.array(coords, order='C')
    coords.shape = (1, -1) # linearizing the array

    # coords is of the form [[a, b, c, ...]]; so coords[0] is the relevant part 
    return coords[0], radii 

def set_velocities(radii):
    vels = []
    for i in np.arange(len(radii)):
        vels.append(rej.set_velocity(radii[i], Mh, a))
        if(i % 1000 == 0): 
            print 'set velocity', i, 'of', N
    vels = np.array(vels, order='C')
    vels.shape = (1, -1)
    return vels[0]

def write_input_file(coords, vels):
    length = len(coords) / 3
    ids = np.arange(1, length + 1, 1)
    masses = np.empty(length)
    masses.fill(float(Mh) / N)
    smooths = np.zeros(length)
    write_snapshot(from_text=False, data_list=[coords, vels, ids, masses,
                                               smooths])

def main():
    init()
    coords, radii = set_positions()
    print "done with the positions"
    vels = set_velocities(radii)
    print "writing output file..."
    write_input_file(coords, vels)

if __name__ == '__main__':
    main()
