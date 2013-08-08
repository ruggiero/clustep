'''
Use: python cluster.py SNAPSHOT, where SNAPSHOT is the name of the output file

Two files should be present in the script's execution folder: cluster_param.txt
and header.txt (see examples)
'''

from snapwrite import * 
import optimized_functions as opt
import numpy as np
import numpy.random as nprand
from scipy import integrate
import os
import sys
from bisect import *

G = 43007.1

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

def potential(radius, Mh, a):
    return -(G * Mh) / (radius + a)

def DF_numerical(E, Mh, a):
    epsilon = -E
    if(epsilon <= 0):
        return 0
    else:
        cte = a / (8**0.5 * (np.pi * G)**3 * Mh**2)
        integral = integrate.quad(opt.aux, 0, epsilon, args = (Mh, a, epsilon), full_output=-1)
        return cte * integral[0]

def interpolate(E, DF_tabulated):
    index = bisect_left(DF_tabulated[:, 0], E)
    if(index >= 999):
        return 0
    else:
        delta_y = DF_tabulated[index + 1][1] - DF_tabulated[index][1]
        delta_x = DF_tabulated[index + 1][0] - DF_tabulated[index][0]
        y = DF_tabulated[index][1] + delta_y / delta_x * (DF_tabulated[index][0] - E)
        return y

def sample_velocity(radius, Mh, a, DF_tabulated):
    phi = potential(radius, Mh, a)
    fmax = interpolate(phi, DF_tabulated)
    vesc = (-2.0 * phi)**0.5
    while(True):
        vx, vy, vz, v2 = opt.random_velocity(vesc)
        y = fmax * nprand.rand()
        if(y < interpolate(phi + 0.5 * v2, DF_tabulated)):
            break
    return [vx, vy, vz]

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

    # coords is in the form [[a, b, c, ...]]; so coords[0] is the relevant part 
    return coords[0], radii 

def set_velocities(radii):
    vels = []
    DF_tabulated = []
    for i in np.linspace(-(G * Mh) / a, 0, 1000):
        DF_tabulated.append([i, DF_numerical(i, Mh, a)])
    DF_tabulated = np.array(DF_tabulated)
    print "done with tabulation"
    for i in np.arange(len(radii)):
        vels.append(sample_velocity(radii[i], Mh, a, DF_tabulated))
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
    write_snapshot(N, from_text=False, data_list=[coords, vels, ids, masses,
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
