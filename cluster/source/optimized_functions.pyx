from libc.math cimport asin, sqrt, pow, M_PI, log

import numpy.random as nprand


cdef double G = 43007.1


# The integrand in the integral for the DF, in case the density profile for both
# the gas and the dark matter is cuspy, in the center.
def aux_cusp_cusp(double r, double M_gas, double a_gas, double M_dm, double a_dm,
            double epsilon):
    cdef double dr_dpsi, d2rho_dr2, drho_dr, d2psi_dr2
    dr_dpsi = 1/(-(M_gas*G)/(r+a_gas)**2-(M_dm*G)/(r+a_dm)**2)
    d2rho_dr2 = (a_dm*M_dm)/(M_PI*r**3*(r+a_dm)**3)+(3*a_dm*M_dm)/(M_PI*r**2*(r+a_dm)**4)+(6*a_dm*M_dm)/(M_PI*r*(r+a_dm)**5)
    drho_dr = -(a_dm*M_dm)/(2*M_PI*r**2*(r+a_dm)**3)-(3*a_dm*M_dm)/(2*M_PI*r*(r+a_dm)**4)
    d2psi_dr2 = -(4*M_gas*G)/(r+a_gas)^3+(2*M_dm*G)/(r+a_dm)^3+(3*M_gas*(2*r+a_gas)*G)/(r+a_gas)^4
    return (dr_dpsi)**2 * (d2rho_dr2 - drho_dr * dr_dpsi * d2psi_dr2) / (epsilon - (G * M_gas/(r + a_gas) + G * M_dm/(r + a_dm)))**0.5 * #tem um psi faltando kkk

# The integrand in the integral for the DF, in case the density profile for the
# gas has a core, and for the dark matter is cuspy, in the center.
def aux_core_cusp(double r, double M_gas, double a_gas, double M_dm, double a_dm,
            double epsilon):
    cdef double dr_dpsi, d2rho_dr2, drho_dr, d2psi_dr2
    dr_dpsi = 1/((M_gas*G)/(r+a_gas)^2-(M_dm*G)/(r+a_dm)^2-(M_gas*(2*r+a_gas)*G)/(r+a_gas)^3)
    d2rho_dr2 = (a_dm*M_dm)/(%Pi*r^3*(r+a_dm)^3)+(3*a_dm*M_dm)/(%Pi*r^2*(r+a_dm)^4)+(6*a_dm*M_dm)/(%Pi*r*(r+a_dm)^5)
    drho_dr = -(a_dm*M_dm)/(2*%Pi*r^2*(r+a_dm)^3)-(3*a_dm*M_dm)/(2*%Pi*r*(r+a_dm)^4)
    d2psi_dr2 = 1/((M_gas*G)/(r+a_gas)^2-(M_dm*G)/(r+a_dm)^2-(M_gas*(2*r+a_gas)*G)/(r+a_gas)^3)
    return (dr_dpsi)**2 * (d2rho_dr2 - drho_dr * dr_dpsi * d2psi_dr2) / (epsilon - (G*M_dm)/(r+a_dm)+(M_gas*(2*r+a_gas)*G)/(2*(r+a_gas)^2))**0.5




# The same, in case there is only dark matter.
def aux_dm(double r, double M, double a, double epsilon):
    cdef double n, d
    n = 6 * r * r + 4 * a * r + a * a
    d = pow(r * (r + a), 3) * sqrt(epsilon - (G * M) / (r + a))
    return n / d


def gas_density(double r, double M_gas, double a_gas):
    return (M_gas * a_gas) / (2 * M_PI * r * pow(r + a_gas, 3))


cdef double cumulative_mass(double r, double M_gas, double a_gas, double M_dm,
                            double a_dm):
    cdef double gas_mass, dm_mass
    gas_mass = (M_gas * r*r) / pow(r + a_gas, 2)
    dm_mass = (M_dm * r*r) / pow(r + a_dm, 2)
    return gas_mass + dm_mass


def T_integrand(double t, double M_gas, double a_gas, double M_dm,
                double a_dm):
    return (G * cumulative_mass(t, M_gas, a_gas, M_dm, a_dm) *
            gas_density(t, M_gas, a_gas)) / (t*t)


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

''' Obsolete stuff

# Works when there is only dark matter.
cdef double DF_analytical(double E, double Mh, double a):
    cdef double q, vg, cte
    if(E >= 0):
        return 0.0
    else:
        vg = sqrt((G * Mh) / a)
        q = sqrt((-a * E) / (G * Mh))
        cte = Mh / (8 * sqrt(2) * pow(M_PI * a * vg, 3))
        return cte * (3 * asin(q) + q * sqrt(1 - q * q) * (1 - 2 * q * q) *
               (8 * pow(q, 4) - 8 * q * q - 3)) / pow(1 - q * q, 2.5)
'''
