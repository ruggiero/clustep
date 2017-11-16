from sys import exit
from sys import path as syspath
from bisect import bisect_left
from os import path
from argparse import ArgumentParser as parser
from ConfigParser import ConfigParser
import time

import numpy as np
import numpy.random as nprand
from scipy import integrate
from scipy.optimize import brentq

from snapwrite import write_snapshot
import optimized_functions as opt
syspath.append(path.join(path.dirname(__file__), '..', 'misc'))
from units import temp_to_internal_energy


G = 44920.0


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
  global gas, dm, output
  global M_dm, a_dm, N_dm, M_gas, a_gas, N_gas, Z
  global truncation_radius, gamma_gas, gamma_dm
  flags = parser(description="Generates an initial conditions file\
                for a galaxy cluster halo simulation.")
  flags.add_argument('--no-dm', help='No dark matter particles in the\
           initial conditions. The dark matter potential is\
           still used when calculating the gas temperatures.',
           action='store_true')
  flags.add_argument('--no-gas', help='Gas is completely ignored, and\
           only dark matter is included.',
           action='store_true')
  flags.add_argument('-o', help='The name of the output file.',
           metavar="init.dat", default="init.dat")
  args = flags.parse_args()
  output = args.o
  if not path.isfile("params_cluster.ini"):
    print "params_cluster.ini missing."
    exit(0)
  if args.no_dm:
    if args.no_gas:
      print "Neither gas or dark matter were selected!"
      exit(0)
    else:
      gas = True
      dm = False
  elif args.no_gas:
    gas = False
    dm = True
  else:
    gas = True
    dm = True
  config = ConfigParser()
  config.read("params_cluster.ini")
  M_dm = config.getfloat('dark_matter', 'M_dm')
  a_dm = config.getfloat('dark_matter', 'a_dm')
  N_dm = config.getint('dark_matter', 'N_dm')
  gamma_dm = config.getfloat('dark_matter', 'gamma_dm')
  if(gas):
    M_gas = config.getfloat('gas', 'M_gas')
    a_gas = config.getfloat('gas', 'a_gas')
    N_gas = config.getint('gas', 'N_gas')
    Z = config.getfloat('gas', 'Z')
    gamma_gas = config.getfloat('gas', 'gamma_gas')
  truncation_radius = config.getfloat('global', 'truncation_radius')


def dehnen_cumulative(r, M, a, gamma):
  return M * (r/(r+float(a)))**(3-gamma)


# Inverse cumulative mass function. Mc is a number between 0 and M.
def dehnen_inverse_cumulative(Mc, M, a, gamma):
  results = []
  for i in Mc:
    results.append(brentq(lambda r: dehnen_cumulative(r, M, a, gamma) - i, 0, 1.0e10))
  return np.array(results)


def potential(r):
  phi = 0
  if(gas):
    if gamma_gas != 2:
      phi += (G*M_gas)/a_gas * (-1.0/(2-gamma_gas)) * (1-(r/(r+float(a_gas)))**(2-gamma_gas))
    else:
      phi += (G*M_gas)/a_gas * np.log(r/(r+float(a_gas)))
  if gamma_dm != 2:
    phi += (G*M_dm)/a_dm * (-1.0/(2-gamma_dm)) * (1-(r/(r+float(a_dm)))**(2-gamma_dm))
  else:
    phi += (G*M_dm)/a_dm * np.log(r/(r+float(a_dm)))
  return phi


def set_positions():
  if(dm):
    # the factor variable restricts the radius to truncation_radius
    factor = dehnen_cumulative(truncation_radius, M_dm, a_dm, gamma_dm)/M_dm
    radii_dm = dehnen_inverse_cumulative(nprand.sample(N_dm) *
                  (M_dm*factor), M_dm, a_dm, gamma_dm)
    thetas = np.arccos(nprand.sample(N_dm)*2 - 1)
    phis = 2 * np.pi * nprand.sample(N_dm)
    xs = radii_dm * np.sin(thetas) * np.cos(phis)
    ys = radii_dm * np.sin(thetas) * np.sin(phis)
    zs = radii_dm * np.cos(thetas)
    coords_dm = np.column_stack((xs, ys, zs))
    coords_dm = np.array(coords_dm, order='C')
    coords_dm.shape = (1, -1) # Linearizing the array.
  if(gas):
    factor = dehnen_cumulative(truncation_radius, M_gas, a_gas, gamma_gas)/M_gas
    radii_gas = dehnen_inverse_cumulative(nprand.sample(N_gas) * (M_gas*factor),
                   M_gas, a_gas, gamma_gas)
    thetas = np.arccos(nprand.sample(N_gas)*2 - 1)
    phis = 2 * np.pi * nprand.sample(N_gas)
    xs = radii_gas * np.sin(thetas) * np.cos(phis)
    ys = radii_gas * np.sin(thetas) * np.sin(phis)
    zs = radii_gas * np.cos(thetas)
    coords_gas = np.column_stack((xs, ys, zs))

    coords_gas = np.array(coords_gas, order='C')
    coords_gas.shape = (1, -1) 
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
  if(gas):
    for i in np.arange(N_gas):
      vels.append([0.0, 0.0, 0.0])
  if(dm):
    DF_tabulated = []
    # This 0.99 avoids numerical problems.
    for i in np.linspace(potential(0)*0.99, 0, 1000):
      DF_tabulated.append([i, DF(i)])
    DF_tabulated = np.array(DF_tabulated)
    print "done with DF tabulation"
    t0 = time.time()
    time_count = 0
    for i in np.arange(len(radii_dm)):
      vels.append(sample_velocity(radii_dm[i], DF_tabulated))
      if(int(time.time()-t0) > time_count):
        print 'set velocity', i, 'of', N_dm
        time_count += 1
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
      if(gamma_gas == 0):
        if(gamma_dm == 0):
          integral = integrate.quad(opt.aux_core_core, limit, np.inf,
            args=(M_gas, a_gas, M_dm, a_dm, epsilon), 
            full_output=-1)
        else:
          integral = integrate.quad(opt.aux_core_cusp, limit, np.inf,
            args=(M_gas, a_gas, M_dm, a_dm, epsilon),
            full_output=-1)
      else:
        if(gamma_dm == 0):
          integral = integrate.quad(opt.aux_cusp_core, limit, np.inf,
            args=(M_gas, a_gas, M_dm, a_dm, epsilon),
            full_output=-1)
        else: 
          integral = integrate.quad(opt.aux_cusp_cusp, limit, np.inf,
            args=(M_gas, a_gas, M_dm, a_dm, epsilon),
            full_output=-1)
    else:
      if(gamma_dm == 0):
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
  meanweight_n = 4.0 / (1 + 3*HYDROGEN_MASSFRAC)
  meanweight_i = 4.0 / (3 + 5*HYDROGEN_MASSFRAC)

  integral = integrate.quad(opt.T_integrand, r, np.inf,
    args=(M_gas, a_gas, M_dm, a_dm, gamma_gas, gamma_dm),
    full_output=-1)
  result = integral[0] / opt.dehnen_density(r, M_gas, a_gas, gamma_gas)

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
    Zs = np.zeros(N_gas)
    Zs.fill(Z)
    if(dm):
      masses_dm = np.empty(N_dm)
      masses_dm.fill(M_dm / N_dm)
      masses = np.concatenate((masses_gas, masses_dm))
      ids = np.arange(1, N_gas + N_dm + 1)
      write_snapshot(n_part=[N_gas, N_dm, 0, 0, 0, 0], outfile=output,
             data_list=[coords, vels, ids, masses, U, rho, smooths, Zs])
    else:
      masses = masses_gas
      ids = np.arange(1, N_gas + 1)
      write_snapshot(n_part=[N_gas, 0, 0, 0, 0, 0], outfile=output,
             data_list=[coords, vels, ids, masses, U, rho, smooths, Zs])
  else:
    ids = np.arange(1, N_dm + 1)
    masses_dm = np.empty(N_dm)
    masses_dm.fill(M_dm / N_dm)
    masses = masses_dm
    write_snapshot(n_part=[0, N_dm, 0, 0, 0, 0], outfile=output, 
            data_list=[coords, vels, ids, masses])


if __name__ == '__main__':
  main()
