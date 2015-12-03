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
M_gas = 1e5*0.15/0.85
time = 0.0


def main():
    input_ = init()
    print "reading..."
    data_halo, data_disk = process_data(input_)
    part_halo, aux_halo = log_partition(data_halo, 1.3)
    density_plot(input_, data_halo, part_halo, aux_halo)
    #density_plot(input_, data_bulge, part_bulge, aux_bulge, bulge=True)
    #velocity_plot(input_, data_disk, 'c')
    #velocity_plot(input_, data_disk, 'r')
    #velocity_plot(input_, data_disk, 'z')
    #gas_velocity_plot(input_, data_gas, 'c')
    #gas_velocity_plot(input_, data_gas, 'r')


def init():
    global halo_core, bulge_core
    flags = parser(description="Plots stuff.")
    flags.add_argument('--halo-core', help='Sets the density profile for the\
                       halo to have a core.', action='store_true')
    flags.add_argument('--bulge-core', help='The same, but for the bulge.',
                       action='store_true')
    flags.add_argument('i', help='The name of the input file.',
                       metavar="file.dat")
    args = flags.parse_args()
    halo_core = args.halo_core
    bulge_core = args.bulge_core
    input_ = args.i
    return input_


def process_data(input_):
    coords_halo = readsnap(input_, 'pos', 'dm')
    coords_bulge = readsnap(input_, 'pos', 'bulge')
    coords_disk = readsnap(input_, 'pos', 'disk')
#    coords_gas = readsnap(input_, 'pos', 'gas')
    vels_disk = readsnap(input_, 'vel', 'disk')
    COD = centers.COD(coords_bulge)
    data_halo = np.array([np.linalg.norm(i-COD) for i in coords_halo])
    data_halo = sorted(data_halo)
    
    coords_disk -= COD
    rhos_disk = (coords_disk[:,0]**2+coords_disk[:,1]**2)**0.5
    phis = np.arctan2(coords_disk[:,1], coords_disk[:,0])
    vphis_disk = vels_disk[:,1]*np.cos(phis) - vels_disk[:,0]*np.sin(phis)
    vrs_disk = vels_disk[:,0]*np.cos(phis) + vels_disk[:,1]*np.sin(phis)
    vzs_disk = vels_disk[:,2]
    data_disk = np.column_stack((rhos_disk, vphis_disk, vrs_disk, vzs_disk))
    #data_disk = np.sort(data_disk)
    for i in range(len(data_disk)):
      print data_disk[i][0], data_disk[i][1]

    return data_halo, data_disk
 

def density(r, core=False, bulge=False):
    if(bulge):
        if(core):
            return (3*M_bulge*a_bulge) / (4*np.pi*(r+a_bulge)**4)
        else:
            if(r == 0):
                return 0
            else:
                return (M_bulge*a_bulge) / (2*np.pi*r*(r+a_bulge)**3)
    else:
        if(core):
            return (3*M_halo*a_halo) / (4*np.pi*(r+a_halo)**4)
        else:
            if(r == 0):
                return 0
            else:
                return (M_halo*a_halo) / (2*np.pi*r*(r+a_halo)**3)




# Given a data vector, in which each element represents a different
# PARTICLe by a list of the form [radius, radial_velocity^2], ordered
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
def density_distribution(data, partition, aux, bulge=False):
    N = len(data)
    distribution = []
    left = 0
    if(bulge):
        cte = (10**10*3*M_bulge) / (4*np.pi*N)
    else:
        cte = (10**10*3*M_halo) / (4*np.pi*N)
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



def density_plot(input_, data, part, aux, bulge=False):
    dist = density_distribution(data, part, aux, bulge)
    x_axis = np.logspace(np.log10(dist[0][0]), np.log10(dist[-1][0]), num=500)
    p1, = plt.plot([i[0] for i in dist], [i[1] for i in dist], '-', color='blue')
    if(bulge):
        p2, = plt.plot(x_axis, [10**10 * density(i, core=bulge_core, bulge=True) for i in x_axis])
    else:
        p2, = plt.plot(x_axis, [10**10 * density(i, core=halo_core) for i in x_axis], color='black')
    plt.legend([p1, p2], [u"Simulação", "Modelo"], loc=1)
    plt.xlabel("Raio (kpc)")
    plt.ylabel("Densidade ( M$_{\odot}$/kpc$^3$)")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1, 10**3])
    plt.ylim([1, 10**9])
    if(bulge):
        plt.title("Bulge, t = %1.2f Gyr" % time)
        plt.savefig(input_ + "-bulge-density.png")
        print "Done with bulge density for " + input_
    else:
        plt.title("Halo, t = %1.2f Gyr" % time)
        plt.savefig(input_ + "-halo-density.png", bbox_inches='tight')
        print "Done with halo density for " + input_
    plt.close()


def velocity_plot(input_, data_disk, comp):
    #formatter = FuncFormatter(lambda x, pos : "%1.2f" % (x / 10**6))
    #ax = plt.subplot(111)
    #ax.yaxis.set_major_formatter(formatter)
    if(comp == 'c'):

#        data = numpy.loadtxt('saida')
#        bins = np.linspace(0, 20, N)
#        digitized = numpy.digitize(data[:,0], bins)
#        bin_means = [data[digitized == i][:,1].mean() for i in range(len(bins))]

        p1, = plt.plot(data_disk[:,0], data_disk[:,1], '.', color='black')
        #plt.legend('Circular velocity', loc=1)
        plt.ylabel("$v_\phi$ (km/s)")
    elif(comp == 'r'):
        p1, = plt.plot(data_disk[:,0], data_disk[:,2], '.', color='black')
        #plt.legend('Radial velocity', loc=1)
        plt.ylabel("$v_R$ (km/s)")
    else:
        p1, = plt.plot(data_disk[:,0], data_disk[:,3], '.', color='black')
        #plt.legend('Radial velocity', loc=1)
        plt.ylabel("$v_z$ (km/s)")
    plt.xlabel("Raio (kpc)")
    plt.xlim([0, 20])
    #plt.title("t = %1.2f Gyr" % time)
    #plt.gcf().subplots_adjust(left=0.17)
    #plt.yscale('log')
    #plt.ylim([0, 100])
    if(comp == 'c'):
        plt.savefig(input_ + "-circular-velocity.png", bbox_inches='tight')
        print "Done with circular velocity for " + input_
    elif(comp == 'r'):
        plt.savefig(input_ + "-radial-velocity.png", bbox_inches='tight')
        print "Done with radial velocity for " + input_
    else:
        plt.savefig(input_ + "-vertical-velocity.png", bbox_inches='tight')
        print "Done with vertical velocity for " + input_
    plt.clf()

def gas_velocity_plot(input_, data_disk, comp):
    #formatter = FuncFormatter(lambda x, pos : "%1.2f" % (x / 10**6))
    #ax = plt.subplot(111)
    #ax.yaxis.set_major_formatter(formatter)
    if(comp == 'c'):
        p1, = plt.plot([i[0] for i in data_disk], [i[1] for i in data_disk], '.', color='black')
        #plt.legend('Circular velocity', loc=1)
        plt.ylabel("$v_\phi$ (km/s)")
    else:
        p1, = plt.plot([i[0] for i in data_disk], [i[2] for i in data_disk], '.')
        #plt.legend('Radial velocity', loc=1)
        plt.ylabel("$v_R$ (km/s)")
    plt.xlabel("Raio (kpc)")
    #plt.xscale('log')
    #plt.title("t = %1.2f Gyr" % time)
    #plt.gcf().subplots_adjust(left=0.17)
    #plt.yscale('log')
    plt.xlim([1, 40])
    #plt.ylim([0, 100])
    if(comp == 'c'):
        plt.savefig(input_ + "-gas-circular-velocity.png", bbox_inches='tight')
        print "Done with gas circular velocity for " + input_
    else:
        plt.savefig(input_ + "-gas-radial-velocity.png")
        print "Done with gas radial velocity for " + input_
    plt.clf()


if __name__ == '__main__':
    main()
