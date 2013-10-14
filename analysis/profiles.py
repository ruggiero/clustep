from sys import path as syspath
from os import path
from bisect import bisect_left
from argparse import ArgumentParser as parser

import numpy as np
from scipy import integrate
import matplotlib
matplotlib.use('Agg') # To be able to plot under an SSH session.
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

syspath.append(path.join(path.dirname(__file__), '..', 'misc'))
syspath.append(path.join(path.dirname(__file__), '..', 'cluster'))
from units import internal_energy_to_temp, temp_to_internal_energy
import optimized_functions as opt
from snapread import read_data, header
import centers


G = 43007.1

a_dm = 200
a_gas = 200
M_dm = 100000
M_gas = 100000
time = 0.0


def density(r, gas=False):
    if(r == 0):
        return 0
    else:
        if(gas):
            if(gas_core):
                return (3 * M_gas * a_gas) / (4 * np.pi * (r+a_gas)**4)
            else:
                return (M_gas * a_gas) / (2 * np.pi * r * (r + a_gas)**3)
        else:
            if(dm_core):
                return (3 * M_dm * a_dm) / (4 * np.pi * (r+a_dm)**4)
            else:
                return (M_dm * a_dm) / (2 * np.pi * r * (r + a_dm)**3)


# Only applies when there is only dark matter, and it follows a Hernquist
# density profile
def radial_velocity(r):
    cte = (G*M_dm) / (12*a_dm)
    t1 = ((12 * r * (r + a_dm)**3) / a_dm**4) * np.log((r + a_dm) / r)
    t2 = -r/(r+a_dm) * (25 + 52*r/a_dm + 42 * (r/a_dm)**2 + 12 * (r/a_dm)**3)
    return (cte * (t1+t2))**0.5


def temperature(r):
    MP_OVER_KB = 121.148
    HYDROGEN_MASSFRAC = 0.76
    meanweight_n = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC)
    meanweight_i = 4.0 / (3 + 5 * HYDROGEN_MASSFRAC)

    integral = integrate.quad(opt.T_integrand,
        r, np.inf, args=(M_gas, a_gas, M_dm, a_dm, int(gas_core), int(dm_core)), full_output=-1)
    result = integral[0] / opt.gas_density(r, M_gas, a_gas, int(gas_core))


    temp_i = MP_OVER_KB * meanweight_i * result
    temp_n = MP_OVER_KB * meanweight_n * result

    if(temp_i > 1.0e4):
        return temp_i
    else:
        return temp_n



# Given a data vector, in which each element represents a different
# PARTICLe by a list of the form [radius, radial_velocity^2], ordered
# according to the radii; and a multiplication factor, returns the right
# indexes of a log partition of the vector. Also returns an auxiliary
# vector, which will be useful in the functions that calculate the
# distribution functions.
def log_partition(data, factor):
    limits = []
    auxiliary = []
    radii = [i[0] for i in data]
    left_limit = 0
    right_limit = 0.01
    left_index = 0
    while(right_limit < 200 * a_gas):
        # Before right_index, everybody is smaller than right_limit.
        right_index = left_index + bisect_left(radii[left_index:], right_limit)
        limits.append(right_index)
        auxiliary.append([right_index - left_index, (right_limit + left_limit) /
                          2])
        left_limit = right_limit
        left_index = right_index
        right_limit *= factor
    return limits, auxiliary


# Returns a list containing elements of the form [radius, density].
def density_distribution(data, partition, aux, gas=False):
    distribution = []
    left = 0
    if(gas):
        cte = (10**10 * 3 * M_gas) / (4 * np.pi * N_gas)
    else:
        cte = (10**10 * 3 * M_dm) / (4 * np.pi * N_dm)
    for j in np.arange(len(partition)):
        right = partition[j]
        count = aux[j][0]
        middle_radius = aux[j][1]
        if(count > 0):
            density = (cte * count) / (data[right][0]**3 - data[left][0]**3)
            distribution.append([middle_radius, density])
        else:
            distribution.append([middle_radius, 0])
        left = right
    return distribution


# Returns a list containing elements of the form [radius, vr].
def radial_velocity_distribution(data, partition, aux):
    distribution = []
    left = 0
    for j in np.arange(len(partition)):
        right = partition[j]
        count = aux[j][0]
        middle_radius = aux[j][1]
        if(count > 0):
            sum_ = sum(i[1] for i in data[left:right])
            distribution.append([middle_radius, (sum_ / count)**0.5])
        else:
            distribution.append([middle_radius, 0])
        left = right
    return distribution


# Returns a list containing elements of the form [radius, vr].
def temperature_distribution(gas_data, partition, aux):
    distribution = []
    left = 0
    for j in np.arange(len(partition)):
        right = partition[j]
        count = aux[j][0]
        middle_radius = aux[j][1]
        if(count > 0):
            sum_ = sum(i[1] for i in gas_data[left:right])
            distribution.append([middle_radius, sum_ / count])
        else:
            distribution.append([middle_radius, 0])
        left = right
    return distribution


def density_plot(input_, data, part, aux, gas=False):
    dist = density_distribution(data, part, aux, gas)
    x_axis = np.logspace(np.log10(dist[0][0]), np.log10(dist[-1][0]), num=500)
    p1, = plt.plot([i[0] for i in dist], [i[1] for i in dist], 'o')
    if(gas):
        p2, = plt.plot(x_axis, [10**10 * density(i, gas=True) for i in x_axis])
    else:
        p2, = plt.plot(x_axis, [10**10 * density(i) for i in x_axis])
    plt.legend([p1, p2], ["Simulation", "Theoretical value"], loc=1)
    plt.xlabel("Radius (kpc)")
    plt.ylabel("Density ( M$_{\odot}$/kpc$^3$)")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1, 10**4])
    plt.ylim([1, 10**10])
    if(gas):
        plt.title("Gas, t = %1.2f Gyr" % time)
        plt.savefig(input_ + "-gas-density.png")
        print "Done with gas density for " + input_
    else:
        plt.title("Dark matter, t = %1.2f Gyr" % time)
        plt.savefig(input_ + "-dm-density.png")
        print "Done with dark matter density for " + input_
    plt.close()


def radial_velocity_plot(input_, data, part, aux):
    dist = radial_velocity_distribution(data, part, aux)
    x_axis = np.logspace(np.log10(dist[0][0]), np.log10(dist[-1][0]), num=500)
    #formatter = FuncFormatter(lambda x, pos : "%1.2f" % (x / 10**6))
    #ax = plt.subplot(111)
    #ax.yaxis.set_major_formatter(formatter)
    p1, = plt.plot([i[0] for i in dist], [i[1] for i in dist], 'o')
    p2, = plt.plot(x_axis, radial_velocity(x_axis))
    plt.legend([p1, p2], ["Simulation", "Theoretical value"], loc=1)
    plt.xlabel("Radius (kpc)")
    plt.ylabel("$\sigma_r$ (km/s)")
    plt.xscale('log')
    plt.title("t = %1.2f Gyr" % time)
    #plt.gcf().subplots_adjust(left=0.17)
    #plt.yscale('log')
    plt.xlim([1, 10**4])
    #plt.ylim([0, 2.5 * 10**6])
    plt.ylim([0, 1650])
    plt.savefig(input_ + "-velocity.png")
    plt.close()
    print "Done with velocity for " + input_


def temperature_plot(input_, gas_data, part, aux):
    dist = temperature_distribution(gas_data, part, aux)
    x_axis = np.logspace(np.log10(dist[0][0]), np.log10(dist[-1][0]), num=500)
    formatter = FuncFormatter(lambda x, pos : "%1.2f" % (x / 10**8))
    ax = plt.subplot(111)
    ax.yaxis.set_major_formatter(formatter)
    p1, = plt.plot([i[0] for i in dist], [i[1] for i in dist], 'o')
    p2, = plt.plot(x_axis, [temperature(i) for i in x_axis])
    plt.legend([p1, p2], ["Simulation", "Theoretical value"], loc=1)
    plt.xlabel("Radius (kpc)")
    plt.ylabel("T (10$^8$K)")
    plt.xscale('log')
    plt.title("t = %1.2f Gyr" % time)
    #plt.gcf().subplots_adjust(left=0.17)
    #plt.yscale('log')
    plt.xlim([1, 10**4])
    #plt.ylim([0, 2.5 * 10**6])
    plt.ylim([0, 3.5e8])
    plt.savefig(input_ + "-temperature.png")
    plt.close()
    print "Done with temperature for " + input_



def init():
    global gas_core, dm_core
    flags = parser(description="Plots stuff.")
    flags.add_argument('--gas-core', help='Sets the density profile for the\
                       gas to have a core.', action='store_true')
    flags.add_argument('--dm-core', help='The same, but for the dark matter.',
                       action='store_true')
    flags.add_argument('-i', help='The name of the input file.',
                       metavar="file.dat", required=True)
    args = flags.parse_args()
    gas_core = args.gas_core
    dm_core = args.dm_core
    input_ = args.i
    return input_


def process_data(input_):
    global gas, time, N_dm, N_gas
    snapshot = open(input_, 'r')
    h = header(snapshot)
    time = h.time
    N_dm = h.n_part_total[1]
    N_gas = h.n_part_total[0]
    if N_gas: 
        gas = True
    else:
        gas = False
    p_list = read_data(snapshot, h)
    snapshot.close()
    data_gas = []
    data_dm = []
    COD = centers.COD(p_list)
    for i in p_list:
        i.pos -= COD
        r = np.linalg.norm(i.pos)
        if(i.ID <= h.n_part_total[0]):
            data_gas.append([r, internal_energy_to_temp(i.U)])
        else:
            vr = np.dot(i.pos, i.vel) / r
            data_dm.append([r, vr**2])
            continue
    del(p_list)
    data_dm = sorted(data_dm)
    data_gas = sorted(data_gas)
    return data_dm, data_gas
 


def main():
    input_ = init()
    print "reading..."
    data_dm, data_gas = process_data(input_)
    part_dm, aux_dm = log_partition(data_dm, 1.3)
    part_gas, aux_gas = log_partition(data_gas, 1.3)
    if gas:
        density_plot(input_, data_dm, part_dm, aux_dm)
        density_plot(input_, data_gas, part_gas, aux_gas, gas=True)
        temperature_plot(input_, data_gas, part_gas, aux_gas)
    else:
        density_plot(input_, data_dm, part_dm, aux_dm, gas=False)
        if not dm_core:
            radial_velocity_plot(input_, data_dm, part_dm, aux_dm)



if __name__ == '__main__':
    main()
