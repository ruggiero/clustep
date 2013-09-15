"""

DESCRIPTION:

Script that generates a snapshot containing a galaxy cluster halo, with
dark matter and gas, both following a Hernquist density profile. The
number of particles of each of these, as well as the parameters a and
M in the density profiles, are defined in the file cluster_param.txt
(see example).

Here, the value for the gravitational constant G is such that the unit
for length is 1.0 kpc, for mass 1.0e10 solar masses, and for velocity
1.0 km/s.


USE:

python cluster.py SNAPSHOT,

where SNAPSHOT is the name of the output file.

"""


from sys import exit, argv, path
from os import path
from bisect import bisect_left

import numpy as np
import numpy.random as nprand
from scipy import integrate

from snapwrite import process_input, write_snapshot
import optimized_functions as opt
path.append(os.path.join(os.path.dirname(__file__), '..', 'misc'))
from units import temp_to_internal_energy


G = 43007.1
gas = False


def main():
    init()
    if(gas):
        cluster_data = generate_cluster_with_gas()
    else:
        cluster_data = generate_cluster_without_gas()
    write_input_file(cluster_data)


def generate_cluster_with_gas():
    print "defining the positions"
    coords, radii_gas, radii_dm = set_positions()
    print "defining the velocities"
    vels = set_velocities(radii_dm)
    print "calculating the temperatures"
    U = set_temperatures(radii_gas)
    print "calculating the densities"
    rho = set_densities(radii_gas)
    print "writing output file..."
    return[coords, vels, U, rho]


def generate_cluster_without_gas():
    print "defining the positions"
    coords, radii_dm = set_positions()
    print "defining the velocities"
    vels = set_velocities(radii_dm)
    print "writing output file..."
    return [coords, vels]



def init():
    global gas
    global M_dm, a_dm, N_dm, M_gas, a_gas, N_gas
    args = argv[1:]
    if not (path.isfile("header.txt") and isfile("cluster_param.txt")):
        print "header.txt or cluster_param.txt missing."
        exit(0)
    if("--gas" in args):
        gas = True
    elif("--dm" in args):
        gas = False
    else:
        print "Please tell whether you want a cluster with or without gas."
        exit(0)
    vars_ = process_input("cluster_param.txt")
    M_dm, a_dm, N_dm = (float(i[0]) for i in vars_[0:3])
    if(gas):
        M_gas, a_gas, N_gas = (float(i[0]) for i in vars_[3:6])


# Inverse cumulative mass function. Depends on both the parameters M and
# a, in the hernquist density profile. Mc is a number between 0 and 1.
def inverse_cumulative(Mc, M, a):
    return (a * ((Mc * M)**0.5 + Mc)) / (M - Mc)


def potential(r):
    phi = - (G * M_dm) / (r + a_dm)
    if(gas):
        phi += -(G * M_gas) / (r + a_gas) 
    return phi


def gas_density(r):
    return (M_gas * a_gas) / (2 * np.pi * r * (r + a_gas)**3)


def cumulative_mass(r):
    gas_mass = (M_gas * r**2) / (r + a_gas)**2
    dm_mass = (M_dm * r**2) / (r + a_dm)**2
    return gas_mass + dm_mass



def set_positions():
    if(gas):

        # The factor M * 200^2 / 201^2 restricts the radius to 200 * a.
        radii_gas = inverse_cumulative(nprand.sample(N_gas) *
                                       ((M_gas * 40000) / 40401), M_gas, a_gas)
        thetas = np.arccos(nprand.sample(N_gas) * 2 - 1)
        phis = 2 * np.pi * nprand.sample(N_gas)
        xs = radii_gas * np.sin(thetas) * np.cos(phis)
        ys = radii_gas * np.sin(thetas) * np.sin(phis)
        zs = radii_gas * np.cos(thetas)
        coords_gas = np.column_stack((xs, ys, zs))

        # Older NumPy versions freak out without this line.
        coords_gas = np.array(coords_gas, order='C')
        coords_gas.shape = (1, -1) # Linearizing the array

    radii_dm = inverse_cumulative(nprand.sample(N_dm) *
                                  ((M_dm * 40000) / 40401), M_dm, a_dm)
    thetas = np.arccos(nprand.sample(N_dm) * 2 - 1)
    phis = 2 * np.pi * nprand.sample(N_dm)
    xs = radii_dm * np.sin(thetas) * np.cos(phis)
    ys = radii_dm * np.sin(thetas) * np.sin(phis)
    zs = radii_dm * np.cos(thetas)
    coords_dm = np.column_stack((xs, ys, zs))
    coords_dm = np.array(coords_dm, order='C')
    coords_dm.shape = (1, -1)
    if(gas):
        coords_total = np.concatenate((coords_gas[0], coords_dm[0]))
        return coords_total, radii_gas, radii_dm
    else:
        return coords_dm[0], radii_dm



def set_velocities(radii_dm):
    vels = []
    DF_tabulated = []
    for i in np.linspace(potential(0), 0, 1000):
        DF_tabulated.append([i, DF_numerical(i)])
    DF_tabulated = np.array(DF_tabulated)
    print "done with tabulation"
    if(gas):
        for i in np.arange(N_gas):
            vels.append([0.0, 0.0, 0.0])
    for i in np.arange(len(radii_dm)):
        vels.append(sample_velocity(radii_dm[i], DF_tabulated))
        if(i % 1000 == 0):
            print 'set velocity', i, 'of', N_dm
    vels = np.array(vels, order='C')
    vels.shape = (1, -1)
    return vels[0]


# Provisory
def DF_numerical(E):
    epsilon = -E
    if(epsilon <= 0):
        return 0
    else:
        if(gas):
            limit1 = (((M_gas**2+2*M_dm*M_gas+M_dm**2)*G**2+epsilon*((2*a_dm-2*a_gas)*M_gas+(2*a_gas-2*a_dm)*M_dm)*G+(a_gas**2-2*a_dm*a_gas+a_dm**2)*epsilon**2)**0.5+(M_gas+M_dm)*G+(-a_gas-a_dm)*epsilon)/(2*epsilon) # This is provisory
            integral = integrate.quad(opt.aux_gas, limit1, np.inf,
                args=(M_gas, a_gas, M_dm, a_dm, epsilon), full_output=-1)
            return -integral[0] / (8**0.5 * np.pi**2)
        else:
            limit1 = (G * M_dm) / epsilon - a_dm
            cte = a_dm / (np.pi**3 * G * 8**0.5)
            integral = integrate.quad(opt.aux_dm, limit1, np.inf,
                args=(M_dm, a_dm, epsilon), full_output=-1)
            return cte * integral[0]


def interpolate(E, DF_tabulated):
    index = bisect_left(DF_tabulated[:, 0], E)
    if(index >= len(DF_tabulated) - 1):
        return 0.0
    elif(index == 0):
        return DF_tabulated[0][1]
    else:
        dy = DF_tabulated[index][1] - DF_tabulated[index - 1][1]
        dx = DF_tabulated[index][0] - DF_tabulated[index - 1][0]
        y = DF_tabulated[index - 1][1] + dy / dx *\
            (DF_tabulated[index - 1][0] - E)
        return y


def sample_velocity(radius, DF_tabulated):
    phi = potential(radius)
    DFmax = interpolate(phi, DF_tabulated)
    vesc = (-2.0 * phi)**0.5
    while(True):
        vx, vy, vz, v2 = opt.random_velocity(vesc)
        y = DFmax * nprand.rand()
        if(y < interpolate(phi + 0.5 * v2, DF_tabulated)):
            break
    return [vx, vy, vz]



def temperature(r):
    mp_over_kb = 121.148
    integral = integrate.quad(opt.T_integrand,
        r, np.inf, args=(M_gas, a_gas, M_dm, a_dm), full_output=-1)
    result = integral[0] / opt.gas_density(r, M_gas, a_gas)

    HYDROGEN_MASSFRAC = 0.76
    meanweight_n = 4.0 / (1 + 3 * HYDROGEN_MASSFRAC)
    meanweight_i = 4.0 / (3 + 5 * HYDROGEN_MASSFRAC)
    
    temp_i = mp_over_kb * meanweight_i * result
    temp_n = mp_over_kb * meanweight_n * result
    if(temp_i > 1.0e4):
        return temp_to_internal_energy(temp_i)
    elif(temp_n < 1.0e4):
        return temp_to_internal_energy(temp_n)
    else:
         return temp_to_internal_energy(0.5 * (temp_i + temp_n)) # ???


def set_temperatures(radii_gas):
    temps = np.zeros(N_gas)
    for i, r in enumerate(radii_gas):
        temps[i] = temperature(r)
    return temps


def set_densities(radii_gas):
    rho = np.zeros(N_gas)
    for i, r in enumerate(radii_gas):
        rho[i] = gas_density(r)
    return rho



def write_input_file(cluster_data):
    coords = cluster_data[0]
    vels = cluster_data[1]
    masses_dm = np.empty(N_dm)
    masses_dm.fill(M_dm / N_dm)
    if(gas):
        ids = np.arange(1, N_gas + N_dm + 1, 1)
        U = cluster_data[2]
        rho = cluster_data[3]
        masses_gas = np.empty(N_gas)
        masses_gas.fill(M_gas / N_gas)
        masses = np.concatenate((masses_gas, masses_dm))
        smooths = np.zeros(N_gas)
        write_snapshot(n_part=[N_gas, N_dm, 0, 0, 0, 0], from_text=False,
                       data_list=[coords, vels, ids, masses, U, rho, smooths])
    else:
        ids = np.arange(1, N_dm + 1, 1)
        masses = masses_dm
        write_snapshot(n_part=[0, N_dm, 0, 0, 0, 0], from_text=False,
                       data_list=[coords, vels, ids, masses])



if __name__ == '__main__':
    main()
