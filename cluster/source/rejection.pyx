from libc.math cimport asin, sqrt, pow, M_PI

def DF(double E, double G, double Mh, double a):
	cdef double q, vg, cte
	if(E>=0):
		return 0.0
	else:
		vg = sqrt((G * Mh) / a)
		q = sqrt((-a * E) / (G * Mh))
		cte = Mh / (8 * sqrt(2) * pow(M_PI * a * vg, 3))
		return cte * (3 * asin(q) + q * sqrt(1 - q * q) * (1 - 2 * q * q) *
               (8 * pow(q, 4) - 8 * q * q - 3)) / pow(1 - q * q, 2.5)
