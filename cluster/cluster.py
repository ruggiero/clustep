"""
DESCRIPTION:

Script that generates a snapshot containing a galaxy cluster halo,
with both dark matter and gas, each of which following a Dehnen density
profile with either gamma equals 1 or 2. The number of particles in each
of these, as well as the parameters a and M in the density profiles,
are defined in the file cluster_param.txt (see example).

Here, the value for the gravitational constant G is such that the unit
for length is 1.0 kpc, for mass 1.0e10 solar masses, and for velocity
1.0 km/s.


USAGE:

cluster.py [-h] [--gas-core] [--dm-core] [--dm-only] [-o init.dat]

optional arguments:
  -h, --help   show this help message and exit
  --gas-core   Sets the density profile for the gas to have a core.
  --dm-core    The same, but for the dark matter.
  --dm-only    Generates an initial conditions file containing only dark
               matter.
  -o init.dat  The name of the output file.

"""


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
    cluster_data = generate_cluster()
    write_input_file(cluster_data)


def generate_cluster():
    if(gas):
        if(dm):
            print "defining the positions"
            coords, radii_gas, radii_dm = set_positions()
            print "defining the velocities"
            vels = set_velocities(radii_dm)
            print "calculating the temperatures"
            U = set_temperatures(radii_gas)
            rho = np.zeros(N_gas)
            print "writing output file..."
            return [coords, vels, U, rho]
        else:
            print "defining the positions"
            coords, radii_gas = set_positions()
            vels = set_velocities()
            print "calculating the temperatures"
            U = set_temperatures(radii_gas)
            rho = np.zeros(N_gas)
            print "writing output file..."
            return [coords, vels, U, rho]
    else:
        print "defining the positions"
        coords, radii_dm = set_positions()
        print "defining the velocities"
        vels = set_velocities(radii_dm)
        print "writing output file..."
        return [coords, vels]


def init():
    global gas, dm, gas_core, dm_core, output
    global M_dm, a_dm, N_dm, M_gas, a_gas, N_gas
    flags = parser(description="Generates an initial conditions file\
                                for a galaxy cluster halo simulation.")
    flags.add_argument('--gas-core', help='Sets the density profile for the\
                       gas to have a core.', action='store_true')
    flags.add_argument('--dm-core', help='The same, but for the dark matter.',
                       action='store_true')
    flags.add_argument('--dm-only', help='Generates an initial conditions\
                       file containing only dark matter.', action='store_true')
    flags.add_argument('--gas-only', help='No dark matter, only gas. The potential
                                           of the dark matter profile supplied is
                                           still used when calculating the
                                           temperatures.',
                       action='store_true')
    flags.add_argument('-o', help='The name of the output file.',
                       metavar="init.dat", default="init.dat")
    args = flags.parse_args()
    gas_core = args.gas_core
    dm_core = args.dm_core
    output = args.o
    if not (path.isfile("header.txt") and path.isfile("cluster_param.txt")):
        print "header.txt or cluster_param.txt missing."
        exit(0)
    if args.dm_only:
        if args.gas_only:
            print "Neither gas or dark matter were selected!"
            exit(0)
        else:
            gas = False
            dm = True
    elif args.gas_only:
        gas = True
        dm = False
    else:
        gas = True
        dm = True
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


def potential(r):
    phi = 0
    if(gas):
        if(gas_core):
            phi -= -(G*M_gas) / 2 * (2*r + a_gas) / (r+a_gas)**2 
        else:
            phi -= -(G*M_gas) / (r+a_gas)
    if(dm_core):
        phi -= (G*M_dm) / 2 * (2*r + a_dm) / (r+a_dm)**2
    else:
        phi -= (G*M_dm) / (r + a_dm)
    return phi


def gas_density(r):
    if(gas_core): 
        return (3*M_gas*a_gas) / (4*np.pi*(r+a_gas)**4)
    else:  
        return (M_gas*a_gas) / (2*np.pi*r*(r+a_gas)**3)


def set_positions():
<<<<<<< HEAD
    # The factor 0.9 cuts the profile at 90% of the mass
    if(dm):
        radii_dm = inverse_cumulative(nprand.sample(N_dm) *
                                     (M_dm * 0.9), M_dm, a_dm, dm_core)
        thetas = np.arccos(nprand.sample(N_dm)*2 - 1)
        phis = 2 * np.pi * nprand.sample(N_dm)
        xs = radii_dm * np.sin(thetas) * np.cos(phis)
        ys = radii_dm * np.sin(thetas) * np.sin(phis)
        zs = radii_dm * np.cos(thetas)

        # Older NumPy versions freak out without this line.
        coords_dm = np.column_stack((xs, ys, zs))
        coords_dm = np.array(coords_dm, order='C')
        coords_dm.shape = (1, -1) # Linearizing the array.
        print "maximum radius (dm): %f" % max(radii_dm)
    if(gas):
        radii_gas = inverse_cumulative(nprand.sample(N_gas) * (M_gas * 0.9),
=======
    # The factor M * 200^2 / 201^2 restricts the radius to 200 * a.
    radii_dm = inverse_cumulative(nprand.sample(N_dm) *
                                  (M_dm*0.9), M_dm, a_dm, dm_core)
    thetas = np.arccos(nprand.sample(N_dm)*2 - 1)
    phis = 2 * np.pi * nprand.sample(N_dm)
    xs = radii_dm * np.sin(thetas) * np.cos(phis)
    ys = radii_dm * np.sin(thetas) * np.sin(phis)
    zs = radii_dm * np.cos(thetas)

    # Older NumPy versions freak out without this line.
    coords_dm = np.column_stack((xs, ys, zs))
    coords_dm = np.array(coords_dm, order='C')
    coords_dm.shape = (1, -1) # Linearizing the array.
    print max(radii_dm)
    if(gas):
        radii_gas = inverse_cumulative(nprand.sample(N_gas) *
                                       ((M_gas*0.9)),
>>>>>>> 4e453f9a58bf0cc5ffe706377661196530e5bb8f
                                       M_gas, a_gas, gas_core)
        thetas = np.arccos(nprand.sample(N_gas) * 2 - 1)
        phis = 2 * np.pi * nprand.sample(N_gas)
        xs = radii_gas * np.sin(thetas) * np.cos(phis)
        ys = radii_gas * np.sin(thetas) * np.sin(phis)
        zs = radii_gas * np.cos(thetas)
        coords_gas = np.column_stack((xs, ys, zs))

        coords_gas = np.array(coords_gas, order='C')
        coords_gas.shape = (1, -1) 
        print "maximum radius (gas): %f" % max(radii_gas)
    if(gas):
        if(dm):
            coords_total = np.concatenate((coords_gas[0], coords_dm[0]))
            return coords_total, radii_gas, radii_dm
        else:
            return coords_gas[0], radii_gas
    else:
        return coords_dm[0], radii_dm


def set_velocities(radii_dm=None):
    vels = []
    DF_tabulated = []
    # This 0.99 avoids numerical problems.
    for i in np.linspace(potential(0)*0.99, 0, 1000):
        DF_tabulated.append([i, DF(i)])
    DF_tabulated = np.array(DF_tabulated)
    print "done with DF tabulation"
    if(gas):
        for i in np.arange(N_gas):
            vels.append([0.0, 0.0, 0.0])
    if(dm):
        for i in np.arange(len(radii_dm)):
            vels.append(sample_velocity(radii_dm[i], DF_tabulated))
            if(i % 1000 == 0):
                print 'set velocity', i, 'of', N_dm
    vels = np.array(vels, order='C')
    vel_COM = sum(vels) # The velocity of the center of mass.
    if(gas):
        vel_COM /= (N_gas + N_dm)
    else:
        vel_COM /= (N_dm)
    for v in vels:
        v -= vel_COM
    vels.shape = (1, -1)
    return vels[0]


def DF(E):
    epsilon = -E
    if(epsilon <= 0):
        return 0
    else:
        # This is r(epsilon), where psi(r) - epsilon = 0.
        limit = brentq(lambda r : -potential(r) - epsilon, 0, 1.0e10)
        if(gas):
            if(gas_core):
                if(dm_core):
                    integral = integrate.quad(opt.aux_core_core, limit, np.inf,
                        args=(M_gas, a_gas, M_dm, a_dm, epsilon), 
                        full_output=-1)
                else:
                    integral = integrate.quad(opt.aux_core_cusp, limit, np.inf,
                        args=(M_gas, a_gas, M_dm, a_dm, epsilon),
                        full_output=-1)
            else:
                if(dm_core):
                    integral = integrate.quad(opt.aux_cusp_core, limit, np.inf,
                        args=(M_gas, a_gas, M_dm, a_dm, epsilon),
                        full_output=-1)
                else: 
                    integral = integrate.quad(opt.aux_cusp_cusp, limit, np.inf,
                        args=(M_gas, a_gas, M_dm, a_dm, epsilon),
                        full_output=-1)
        else:
            if(dm_core):
                integral = integrate.quad(opt.aux_core, limit, np.inf,
                    args=(M_dm, a_dm, epsilon), full_output=-1) 
            else:
                return opt.DF_analytical(E, M_dm, a_dm)
        return -integral[0] / (8**0.5 * np.pi**2)


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


def sample_velocity(radius, DF_tabulated):
    phi = potential(radius)
    DFmax = interpolate(phi, DF_tabulated)
    vesc = (-2.0*phi)**0.5
    while(True):
        vx, vy, vz, v2 = opt.random_velocity(vesc)
        y = DFmax * nprand.rand()
        # We don't want particles unbound to the potential
        if(y < interpolate(phi + 0.5*v2, DF_tabulated) and v2**0.5 < 0.95*vesc):
            break
    return [vx, vy, vz]


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
    if(gas):
        U = cluster_data[2]
        rho = cluster_data[3]
        smooths = np.zeros(N_gas)
        masses_gas = np.empty(N_gas)
        masses_gas.fill(M_gas / N_gas)
        if(dm):
            masses_dm = np.empty(N_dm)
            masses_dm.fill(M_dm / N_dm)
            masses = np.concatenate((masses_gas, masses_dm))
            ids = np.arange(1, N_gas + N_dm + 1)
            write_snapshot(n_part=[N_gas, N_dm, 0, 0, 0, 0], from_text=False,
                           data_list=[coords, vels, ids, masses, U, rho, smooths])
        else:
            masses = masses_gas
            ids = np.arange(1, N_gas + 1)
            write_snapshot(n_part=[N_gas, 0, 0, 0, 0, 0], from_text=False,
                           data_list=[coords, vels, ids, masses, U, rho, smooths])
    else:
        ids = np.arange(1, N_dm + 1)
        masses_dm = np.empty(N_dm)
        masses_dm.fill(M_dm / N_dm)
        masses = masses_dm
        write_snapshot(n_part=[0, N_dm, 0, 0, 0, 0], from_text=False,
                       outfile=output, data_list=[coords, vels, ids, masses])


if __name__ == '__main__':
    main()
