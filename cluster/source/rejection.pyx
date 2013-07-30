from libc.math cimport asin, sqrt, pow, M_PI
import numpy.random as nprand

cdef double G = 43007.1

cdef double potential(double radius, double Mh, double a):
    return -(G * Mh) / (radius + a)

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

#cdef double DF_numerical(double radius, double v, doubl q

def set_velocity(double radius, double Mh, double a):
    cdef double fmax, y, phi
    cdef double vx, vy, vz, v2, vesc
    phi = potential(radius, Mh, a)
    fmax = DF_analytical(phi, Mh, a)
    vesc = sqrt(-2.0 * phi)
    while(True):
        vx = vesc * (nprand.rand() * 2 - 1)
        vy = vesc * (nprand.rand() * 2 - 1)
        vz = vesc * (nprand.rand() * 2 - 1)
        v2 = vx*vx + vy*vy + vz*vz
        y = fmax * nprand.rand()
        if(y < DF_analytical(phi + 0.5 * v2, Mh, a)):
            break
    return [vx, vy, vz]
