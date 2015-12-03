# -*- coding: utf-8 -*-


from sys import path as syspath
from os import path
from bisect import bisect_left
from argparse import ArgumentParser as parser

import numpy as np
from numpy import pi, cos, sin, arctan
from scipy import integrate
import matplotlib
matplotlib.use('Agg') # To be able to plot under an SSH session.
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from pygadgetreader import *
import centers


# These variables must be updated manually.
a_halo = 200
M_halo = 100000
a_gas = 200
M_gas = 1e5*0.15/0.85
time = 0.0


def main():
    input_ = init()
    print "reading..."
    npart = readheader(input_, 'npartTotal')
    if(npart[0] > 0):
      data_gas = process_data(input_, 'gas')
      part_gas, aux_gas = log_partition(data_gas, 1.3)
      density_plot(input_, data_gas, part_gas, aux_gas, 'gas')
    if(npart[1] > 0):
      data_dm = process_data(input_, 'dm')
      part_dm, aux_dm = log_partition(data_dm, 1.3)
      density_plot(input_, data_dm, part_dm, aux_dm, 'dm')


def init():
    global dm_core, gas_core
    flags = parser(description="Plots stuff.")
    flags.add_argument('--dm-core', help='Sets the density profile for the\
                       dark matter to have a core.', action='store_true')
    flags.add_argument('--gas-core', help='The same, but for the bulge.',
                       action='store_true')
    flags.add_argument('i', help='The name of the input file.',
                       metavar="file.dat")
    args = flags.parse_args()
    dm_core = args.dm_core
    gas_core = args.gas_core
    input_ = args.i
    return input_


def process_data(input_, component):
  coords = readsnap(input_, 'pos', component)
  data = np.array([np.linalg.norm(i) for i in coords])
  data = sorted(data)
  return data
 

def density(r, M, a, core=False):
  if(core):
    return (3*M*a) / (4*np.pi*(r+a)**4)
  else:
    if(r == 0):
      return 0
    else:
      return (M*a) / (2*np.pi*r*(r+a)**3)


# Given a data vector, in which each element represents a different
# particle by a list of the form [radius, radial_velocity^2], ordered
# according to the radii; and a multiplication factor, returns the right
# indexes of a log partition of the vector. Also returns an auxiliary
# vector, which will be useful in the functions that calculate the
# distribution functions.
def log_partition(data, factor):
    limits = []
    auxiliary = []
    left_limit = 0
    right_limit = 0.01
    left_index = 0
    while(right_limit < 200 * a_halo):
        # Before right_index, everybody is smaller than right_limit.
        right_index = left_index + bisect_left(data[left_index:], right_limit)
        limits.append(right_index)
        auxiliary.append([right_index - left_index, (right_limit + left_limit) /
                          2])
        left_limit = right_limit
        left_index = right_index
        right_limit *= factor
    return limits, auxiliary


# Returns a list containing elements of the form [radius, density].
def density_distribution(data, partition, aux, M):
    N = len(data)
    distribution = []
    left = 0
    cte = (10**10*3*M) / (4*np.pi*N)
    for j in np.arange(len(partition)):
        right = partition[j]
        if(right >= len(data)):
            break
        count = aux[j][0]
        middle_radius = aux[j][1]
        if(count > 0):
            density = (cte * count) / (data[right]**3 - data[left]**3)
            distribution.append([middle_radius, density])
        else:
            distribution.append([middle_radius, 0])
        left = right
    return distribution


def density_plot(input_, data, part, aux, name):
    if(name == 'gas'):
      M = M_gas
      a = a_gas
      core = gas_core
    else:
      M = M_dm
      a = a_dm
      core = dm_core
    dist = density_distribution(data, part, aux, M)
    x_axis = np.logspace(np.log10(dist[0][0]), np.log10(dist[-1][0]), num=500)
    p1, = plt.plot([i[0] for i in dist], [i[1] for i in dist], '-', color='blue')
    p2, = plt.plot(x_axis, [10**10 * density(i, M, a, core) for i in x_axis], color='black')
    plt.legend([p1, p2], [u"Data", "Model"], loc=1)
    plt.xlabel("Radius (kpc)")
    plt.ylabel("Density ( M$_{\odot}$/kpc$^3$)")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1, 3e3])
    plt.ylim([1, 10**9])
    plt.title(name + ", t = %1.2f Gyr" % time)
    plt.savefig(input_ + "-" + name + "-density.png", bbox_inches='tight')
    print "Done with " + name + " density for " + input_
    plt.close()


if __name__ == '__main__':
    main()
