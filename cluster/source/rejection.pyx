from libc.math cimport asin, sqrt, pow, M_PI

'''
def DF(double E, double G, double Mh, double a):
	cdef double q, vg, cte, c0
	if(E >= 0):
		return 0.0
	else:
		q = sqrt(-a * E / G / Mh)
		vg = sqrt(G * Mh / a)
		cte = Mh / (8.0 * sqrt(2.0) * pow(M_PI, 3)) / pow(a * vg, 3)
		c0 = sqrt(1 - q*q)
		return cte / pow(c0, 5) * 3 * asin(q) + q * c0 * (1 - 2 * q*q) *\
			   (8 * q*q*q*q - 8 * q**q - 3)
'''
def DF(double E, double G, double Mh, double a):
	cdef double q, vg, cte, c0
	if(E>=0):
		return 0.0
	else:
		q = (-a * E / G / Mh)**0.5
		vg = (G * Mh / a)**0.5
		return (Mh * (3 * asin(q) + q * (1 - q**2)**0.5 * (1 - 2 * q**2) * (8 * q**4 - 8 * q**2 - 3))) / (8 * 2**0.5 * M_PI**3 * a**3 * vg**3 * (1 - q**2)**2.5);
