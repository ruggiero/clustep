from libc.math cimport asin, sqrt, pow, M_PI
import numpy.random as nprand

cdef double G = 43007.1

'''
# The integrand in the integral for the DF
def aux(double psi, double Mh, double a, double epsilon):
    cdef double n1, n2, n3, d1
    n1 = 6 * pow(G * Mh, 2) 
    n2 = 8 * a * psi * G * Mh
    n3 = 3 * pow(a * psi, 2)
    d1 = pow(G * Mh - a * psi, 3) 
    return (pow(psi, 2) * (n1 - n2 + n3)) / (d1 * sqrt(epsilon - psi))
'''

# The integrand in the integral for the DF
def aux(double r, double M, double a, double epsilon):
    cdef double n, d
    n = 6 * r * r + 4 * a * r + a * a
    d = pow(r * (r + a), 3) * sqrt(epsilon - (G * M) / (r + a))
    return n / d

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

def DF_numerical(double E, double Mh, double a):
    cdef double epsilon, cte
    epsilon = -E
    if(epsilon <= 0):
        return 0
    else:
        cte = a / (sqrt(8) * pow(M_PI * G, 3) * pow(Mh, 2))
        integral = integrate.quad(aux, 0, epsilon, args = (Mh, a, epsilon))
        return cte * integral[0]

def set_velocity(double radius, double Mh, double a,
                 np.ndarray[np.float64_t] DF_tabulated,
                 np.ndarray[np.float64_t] DF_pos):
    cdef double fmax, y, phi
    cdef double vx, vy, vz, v2, vesc
    phi = potential(radius, Mh, a)
    fmax = search_array(phi, DF_tabulated, DF_pos)
    vesc = sqrt(-2.0 * phi)
    while(True):
        vx = vesc * (nprand.rand() * 2 - 1)
        vy = vesc * (nprand.rand() * 2 - 1)
        vz = vesc * (nprand.rand() * 2 - 1)
        v2 = vx*vx + vy*vy + vz*vz
        y = fmax * nprand.rand()
        if(y < search_array(phi + 0.5 * v2, DF_tabulated, DF_pos)):
            break
    return [vx, vy, vz]

cdef double potential(double radius, double Mh, double a):
    return -(G * Mh) / (radius + a)
'''
