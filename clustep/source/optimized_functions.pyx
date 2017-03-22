'''

All these huge formulas are integrands for the calculation of the
distribution function, given by the Eddington Formula.  More precisely,
they are

(dr/dpsi)^2 * (d^2rho/dr^2 - drho/dr * dr/dpsi * d^2psi/dr^2) * dpsi/dr /
(epsilon - psi)^0.5

'''

from libc.math cimport asin, sqrt, pow, M_PI, log

import numpy.random as nprand


cdef double G = 44920.0

# The integrand in the integral for the DF, in case the density profile for the
# the gas is cuspy, while the one for the dark matter has a core. All the following
# five functions follow the same logic.
def aux_cusp_core(double r, double M_gas, double a_gas, double M_dm,
                  double a_dm, double epsilon):
    return -(3*sqrt(2)*a_dm*M_dm*(r+a_gas)*(3*M_gas*r**4+3*M_dm*r**4+5*a_gas*M_gas*r**3+7*a_dm*M_gas*r**3+9*a_gas*M_dm*r**3+a_dm*M_dm*r**3+15*a_dm*a_gas*M_gas*r**2+3*a_dm**2*M_gas*r**2+9*a_gas**2*M_dm*r**2+3*a_dm*a_gas*M_dm*r**2+15*a_dm**2*a_gas*M_gas*r-3*a_dm**3*M_gas*r+3*a_gas**3*M_dm*r+3*a_dm*a_gas**2*M_dm*r+5*a_dm**3*a_gas*M_gas-2*a_dm**4*M_gas+a_dm*a_gas**3*M_dm))/(M_PI*(r+a_dm)*(M_gas*r**3+M_dm*r**3+3*a_dm*M_gas*r**2+2*a_gas*M_dm*r**2+3*a_dm**2*M_gas*r+a_gas**2*M_dm*r+a_dm**3*M_gas)**2*abs(r+a_dm)*G*sqrt(-(((2*M_gas+2*M_dm)*r**2+(4*a_dm*M_gas+(2*a_gas+a_dm)*M_dm)*r+2*a_dm**2*M_gas+a_dm*a_gas*M_dm)*G-2*epsilon*r**3+(-2*a_gas-4*a_dm)*epsilon*r**2+(-4*a_dm*a_gas-2*a_dm**2)*epsilon*r-2*a_dm**2*a_gas*epsilon)/(r+a_gas)))


def aux_cusp_cusp(double r, double M_gas, double a_gas, double M_dm,
                  double a_dm, double epsilon):
    return -(a_dm*M_dm*(r+a_gas)*(6*M_gas*r**5+6*M_dm*r**5+10*a_gas*M_gas*r**4+12*a_dm*M_gas*r**4+18*a_gas*M_dm*r**4+4*a_dm*M_dm*r**4+25*a_dm*a_gas*M_gas*r**3+6*a_dm**2*M_gas*r**3+18*a_gas**2*M_dm*r**3+12*a_dm*a_gas*M_dm*r**3+a_dm**2*M_dm*r**3+21*a_dm**2*a_gas*M_gas*r**2+6*a_gas**3*M_dm*r**2+12*a_dm*a_gas**2*M_dm*r**2+3*a_dm**2*a_gas*M_dm*r**2+7*a_dm**3*a_gas*M_gas*r+4*a_dm*a_gas**3*M_dm*r+3*a_dm**2*a_gas**2*M_dm*r+a_dm**4*a_gas*M_gas+a_dm**2*a_gas**3*M_dm))/(M_PI*r**3*(r+a_dm)**3*(M_gas*r**2+M_dm*r**2+2*a_dm*M_gas*r+2*a_gas*M_dm*r+a_dm**2*M_gas+a_gas**2*M_dm)**2*G*sqrt(-(((M_gas+M_dm)*r+a_dm*M_gas+a_gas*M_dm)*G-epsilon*r**2+(-a_gas-a_dm)*epsilon*r-a_dm*a_gas*epsilon)/((r+a_dm)*(r+a_gas))))


def aux_core_cusp(double r, double M_gas, double a_gas, double M_dm,
                  double a_dm, double epsilon):
    return -(a_dm*M_dm*(r+a_gas)**2*(12*M_gas*r**6+12*M_dm*r**6+24*a_gas*M_gas*r**5+24*a_dm*M_gas*r**5+48*a_gas*M_dm*r**5+8*a_dm*M_dm*r**5+63*a_dm*a_gas*M_gas*r**4+12*a_dm**2*M_gas*r**4+72*a_gas**2*M_dm*r**4+32*a_dm*a_gas*M_dm*r**4+2*a_dm**2*M_dm*r**4+57*a_dm**2*a_gas*M_gas*r**3+48*a_gas**3*M_dm*r**3+48*a_dm*a_gas**2*M_dm*r**3+8*a_dm**2*a_gas*M_dm*r**3+21*a_dm**3*a_gas*M_gas*r**2+12*a_gas**4*M_dm*r**2+32*a_dm*a_gas**3*M_dm*r**2+12*a_dm**2*a_gas**2*M_dm*r**2+3*a_dm**4*a_gas*M_gas*r+8*a_dm*a_gas**4*M_dm*r+8*a_dm**2*a_gas**3*M_dm*r+2*a_dm**2*a_gas**4*M_dm)*abs(r+a_gas))/(sqrt(2)*M_PI*r**3*(r+a_dm)**3*(M_gas*r**3+M_dm*r**3+2*a_dm*M_gas*r**2+3*a_gas*M_dm*r**2+a_dm**2*M_gas*r+3*a_gas**2*M_dm*r+a_gas**3*M_dm)**2*G*sqrt(-(((2*M_gas+2*M_dm)*r**2+((a_gas+2*a_dm)*M_gas+4*a_gas*M_dm)*r+a_dm*a_gas*M_gas+2*a_gas**2*M_dm)*G-2*epsilon*r**3+(-4*a_gas-2*a_dm)*epsilon*r**2+(-2*a_gas**2-4*a_dm*a_gas)*epsilon*r-2*a_dm*a_gas**2*epsilon)/(r+a_dm)))


def aux_core_core(double r, double M_gas, double a_gas, double M_dm,
                  double a_dm, double epsilon):
    return -(3*sqrt(2)*a_dm*M_dm*(r+a_gas)**2*(3*M_gas*r**5+3*M_dm*r**5+6*a_gas*M_gas*r**4+7*a_dm*M_gas*r**4+12*a_gas*M_dm*r**4+a_dm*M_dm*r**4+19*a_dm*a_gas*M_gas*r**3+3*a_dm**2*M_gas*r**3+18*a_gas**2*M_dm*r**3+4*a_dm*a_gas*M_dm*r**3+21*a_dm**2*a_gas*M_gas*r**2-3*a_dm**3*M_gas*r**2+12*a_gas**3*M_dm*r**2+6*a_dm*a_gas**2*M_dm*r**2+9*a_dm**3*a_gas*M_gas*r-2*a_dm**4*M_gas*r+3*a_gas**4*M_dm*r+4*a_dm*a_gas**3*M_dm*r+a_dm**4*a_gas*M_gas+a_dm*a_gas**4*M_dm)*abs(r+a_gas))/(M_PI*r**2*(r+a_dm)*(M_gas*r**3+M_dm*r**3+3*a_dm*M_gas*r**2+3*a_gas*M_dm*r**2+3*a_dm**2*M_gas*r+3*a_gas**2*M_dm*r+a_dm**3*M_gas+a_gas**3*M_dm)**2*abs(r+a_dm)*G*sqrt(-((2*M_gas+2*M_dm)*r**3+((a_gas+4*a_dm)*M_gas+(4*a_gas+a_dm)*M_dm)*r**2+((2*a_dm*a_gas+2*a_dm**2)*M_gas+(2*a_gas**2+2*a_dm*a_gas)*M_dm)*r+a_dm**2*a_gas*M_gas+a_dm*a_gas**2*M_dm)*G+2*epsilon*r**4-(-4*a_gas-4*a_dm)*epsilon*r**3-(-2*a_gas**2-8*a_dm*a_gas-2*a_dm**2)*epsilon*r**2-(-4*a_dm*a_gas**2-4*a_dm**2*a_gas)*epsilon*r+2*a_dm**2*a_gas**2*epsilon))


# The same, in case there is only dark matter, with a cuspy profile.
def aux_cusp(double r, double M, double a, double epsilon):
    cdef double n, d
    n = 6 * r * r + 4 * a * r + a * a
    d = pow(r * (r + a), 3) * sqrt(epsilon - (G * M) / (r + a))
    return -a / (G * M_PI) * n / d


# Now with a core profile.
def aux_core(double r, double M_dm, double a_dm, double epsilon):
    return -(3*sqrt(2)*a_dm*(3*r+a_dm))/(M_PI*r**2*(r+a_dm)*abs(r+a_dm)*G*sqrt(-(2*M_dm*r+a_dm*M_dm)*G+2*epsilon*r**2+4*a_dm*epsilon*r+2*a_dm**2*epsilon))


def gas_density(double r, double M_gas, double a_gas, int gas_core):
    if(gas_core):
        return (3 * M_gas * a_gas) / (4 * M_PI * (r+a_gas)**4)
    else:
        return (M_gas * a_gas) / (2 * M_PI * r * pow(r + a_gas, 3))


cdef double cumulative_mass(double r, double M_gas, double a_gas, double M_dm,
                            double a_dm, int gas_core, int dm_core):
    cdef double gas_mass, dm_mass
    if(gas_core):
        gas_mass = M_gas * pow(r / (r + a_gas), 3)
    else:
        gas_mass = (M_gas * r*r) / pow(r + a_gas, 2)
    if(dm_core):
        dm_mass = M_dm * pow(r / (r + a_dm), 3)
    else:
        dm_mass = (M_dm * r*r) / pow(r + a_dm, 2)
    return gas_mass + dm_mass


def T_integrand(double t, double M_gas, double a_gas, double M_dm,
                double a_dm, int gas_core, int dm_core):
    return (G*cumulative_mass(t, M_gas, a_gas, M_dm, a_dm, gas_core, dm_core) *
            gas_density(t, M_gas, a_gas, gas_core)) / (t*t)


def random_velocity(double vesc):
    cdef double vx, vy, vz, v2, vesc2
    vesc2 = vesc * vesc
    while(True):
        vx = vesc * (nprand.rand() * 2 - 1)
        vy = vesc * (nprand.rand() * 2 - 1)
        vz = vesc * (nprand.rand() * 2 - 1)
        v2 = vx*vx + vy*vy + vz*vz
        if(v2 < vesc2):
            break
    return vx, vy, vz, v2


# Works when there is only dark matter.
def DF_analytical(double E, double Mh, double a):
    cdef double q, vg, cte
    if(E >= 0):
        return 0.0
    else:
        vg = sqrt((G * Mh) / a)
        q = sqrt((-a*E) / (G*Mh))
        cte = Mh / (8 * sqrt(2) * pow(M_PI*a*vg, 3))
        return cte * (3*asin(q) + q*sqrt(1-q*q) * (1-2*q*q) *
               (8*pow(q, 4) - 8*q*q - 3)) / pow(1 - q*q, 2.5)
