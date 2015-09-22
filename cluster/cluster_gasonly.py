from sys import exit
from sys import path as syspath
from bisect import bisect_left
from os import path
from argparse import ArgumentParser as parser

import numpy as np
import numpy.random as nprand
from scipy import integrate
from scipy.optimize import brentq

from snapwrite import process_input, write_snapshot
import optimized_functions as opt
syspath.append(path.join(path.dirname(__file__), '..', 'misc'))
from units import temp_to_internal_energy


G = 43007.1


def main():
    init()
    if(gas):
        cluster_data = generate_cluster_with_gas()
    else:
        cluster_data = generate_cluster_without_gas()
    write_input_file(cluster_data)


def generate_cluster_with_gas():
    print "defining the positions"
    coords, radii_gas = set_positions()
    print "calculating the temperatures"
    U = set_temperatures(radii_gas)
    rho = np.zeros(N_gas)
    vels = np.zeros(3*N_gas)
    print "writing output file..."
    return [coords, vels, U, rho]


def init():
    global gas, gas_core, dm_core
    global M_dm, a_dm, N_dm, M_gas, a_gas, N_gas
    flags = parser(description="Generates an initial conditions file\
                                for a galaxy cluster halo simulation.")
    flags.add_argument('--dm-core', help='The same, but for the dark matter.',
                       action='store_true')
    flags.add_argument('--gas-core', help='Sets the density profile for the\
                       gas to have a core.', action='store_true')
    flags.add_argument('-o', help='The name of the output file.',
                       metavar="init.dat", default="init.dat")
    args = flags.parse_args()
    gas_core = args.gas_core
    dm_core = args.dm_core
    if not (path.isfile("header.txt") and path.isfile("cluster_param.txt")):
        print "header.txt or cluster_param.txt missing."
        exit(0)
    gas = True
    vars_ = process_input("cluster_param.txt")
    M_dm, a_dm, N_dm = (float(i[0]) for i in vars_[0:3])
    if(gas):
        M_gas, a_gas, N_gas = (float(i[0]) for i in vars_[3:6])


# Inverse cumulative mass function. Depends on both the parameters M and
# a, in the Dehnen density profile. Mc is a number between 0 and M.
def inverse_cumulative(Mc, M, a, core):
    if(core):
        return ((a * (Mc**(2/3.)*M**(4/3.) + Mc*M + Mc**(4/3.)*M**(2/3.))) /
                   (Mc**(1/3.) * M**(2/3.) * (M-Mc)))
    else:
        return (a * ((Mc*M)**0.5 + Mc)) / (M-Mc)


def gas_density(r):
    if(gas_core): 
        return (3*M_gas*a_gas) / (4*np.pi*(r+a_gas)**4)
    else:  
        return (M_gas*a_gas) / (2*np.pi*r*(r+a_gas)**3)


# Positions are restricted to the radius where 90% of the mass is
# at, so particles don't go too far
def set_positions():
    radii_gas = inverse_cumulative(nprand.sample(N_gas) *
                                   (M_gas*0.9), M_gas, a_gas, gas_core)
    thetas = np.arccos(nprand.sample(N_gas) * 2 - 1)
    phis = 2 * np.pi * nprand.sample(N_gas)
    xs = radii_gas * np.sin(thetas) * np.cos(phis)
    ys = radii_gas * np.sin(thetas) * np.sin(phis)
    zs = radii_gas * np.cos(thetas)
    coords_gas = np.column_stack((xs, ys, zs))

    coords_gas = np.array(coords_gas, order='C')
    coords_gas.shape = (1, -1) 
    return coords_gas[0], radii_gas



def interpolate(E, DF_tabulated):
    index = bisect_left(DF_tabulated[:, 0], E)
    if(index >= len(DF_tabulated) - 1):
        return 0.0
    elif(index == 0):
        return DF_tabulated[0][1]
    else:
        dy = DF_tabulated[index][1] - DF_tabulated[index - 1][1]
        dx = DF_tabulated[index][0] - DF_tabulated[index - 1][0]
        y = (DF_tabulated[index - 1][1] + dy / dx *
             (DF_tabulated[index - 1][0] - E))
        return y


def temperature(r):
    MP_OVER_KB = 121.148
    HYDROGEN_MASSFRAC = 0.76
    meanweight_n = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC)
    meanweight_i = 4.0 / (3 + 5 * HYDROGEN_MASSFRAC)

    integral = integrate.quad(opt.T_integrand, r, np.inf,
        args=(M_gas, a_gas, M_dm, a_dm, int(gas_core), int(dm_core)),
        full_output=-1)
    result = integral[0] / opt.gas_density(r, M_gas, a_gas, int(gas_core))

    temp_i = MP_OVER_KB * meanweight_i * result
    temp_n = MP_OVER_KB * meanweight_n * result

    if(temp_i > 1.0e4):
        return temp_to_internal_energy(temp_i)
    else:
        return temp_to_internal_energy(temp_n)


def set_temperatures(radii_gas):
    temps = np.zeros(N_gas)

    T_tabulated = []
    # This 0.99 avoids numerical problems.
    for r in np.logspace(-1, np.log10(200*a_gas), 1000):
        T_tabulated.append([r, temperature(r)])
    T_tabulated = np.array(T_tabulated)
    for i, r in enumerate(radii_gas):
        temps[i] = interpolate(r, T_tabulated)
    return temps


def write_input_file(cluster_data):
    coords = cluster_data[0]
    vels = cluster_data[1]
    ids = np.arange(1, N_gas + 1, 1)
    U = cluster_data[2]
    rho = cluster_data[3]
    masses_gas = np.empty(N_gas)
    masses_gas.fill(M_gas / N_gas)
    smooths = np.zeros(N_gas)
    write_snapshot(n_part=[N_gas, 0, 0, 0, 0, 0], from_text=False,
                   data_list=[coords, vels, ids, masses_gas, U, rho, smooths])


if __name__ == '__main__':
    main()
