'''
Use: python cluster.py SNAPSHOT, where SNAPSHOT is the name of the output file

Two files should be present in the script's execution folder: cluster_param.txt
and header.txt (see examples)

'''

from snapwrite import *
import numpy as np
import numpy.random as nprand
import os
import sys
import shutil

G = 43007.1
folder = "/tmp/.cluster_temp/"
saida = sys.argv[1]

# Extracting the values of the global variables from cluster_param.txt,
# creating a temporary folder where the input files for snapwrite.py will be
# created and copying the file header.txt to this folder
def init():
	global Mh, a, N, vg
	if not (os.path.isfile("header.txt") and os.path.isfile(\
			"cluster_param.txt")):
		print "The following parameter files are required: header.txt and\
			   cluster_param.txt. Please check."
		sys.exit(0)
	if not os.path.exists(folder):
		os.makedirs(folder)
	shutil.copyfile("header.txt", folder + "header.txt")
	x = process_input("cluster_param.txt")
	Mh, a, N = float(x[0]), float(x[1]), float(x[2])
	vg = ((G * Mh) / a)**0.5

def inverse_cumulative(Mc):
	return (a * ((Mc * Mh)**0.5 + Mc)) / (Mh - Mc)

def potential(radius):
	return -(G * Mh) / (radius + a)

# The distribution function
def DF(E):
	if(E >= 0):
		return 0
	else:
		q = (-(a * E) / (G * Mh))**0.5
		return Mh * (3 * np.arcsin(q) + q * (1 - q**2)**0.5 * (1 - 2 *
			   q**2) * (8 * q**4 - 8 * q**2 - 3)) / (8 * 2**0.5 *
			   np.pi**3 * a**3 * vg**3 * (1 - q**2)**2.5)

def set_positions():
	
	# This factor Mh 200^2 / 201^2 is for restricting the radius to 200a
	radii = inverse_cumulative(nprand.sample(N) * ((Mh * 200**2) / 201**2))
	thetas = np.arccos(nprand.sample(N) * 2 - 1)
	phis = 2 * np.pi * nprand.sample(N)
	xs = radii * np.sin(thetas) * np.cos(phis)
	ys = radii * np.sin(thetas) * np.sin(phis)
	zs = radii * np.cos(thetas)
	return np.column_stack((xs, ys, zs)), radii

# The most intensive function in this script
def set_velocities(radii):
	pots = potential(radii)
	vels = []
	for i in np.arange(len(pots)):
		fmax = DF(pots[i])

		# random coordinates in the rectangle (pots[i], 0) x (0, fmax)
		# will be rejected if y > DF(x), as the rejection algorithm requires
		x = pots[i] - (nprand.rand() * pots[i])
		y = nprand.rand() * fmax
		while(y > DF(x)):
			x = pots[i] - (nprand.rand() * pots[i])
			y = nprand.rand() * fmax
		vels.append((2 * (x - pots[i]))**0.5)
	thetas = np.arccos(nprand.sample(N) * 2 - 1)
	phis = 2 * np.pi * nprand.sample(N)
	vxs = vels * np.sin(thetas) * np.cos(phis)
	vys = vels * np.sin(thetas) * np.sin(phis)
	vzs = vels * np.cos(thetas)
	return np.column_stack((vxs, vys, vzs))

def write_input_files(coords, vels):
	length = len(coords)
	np.savetxt(folder + "position.txt", coords, delimiter='\t', fmt="%f")
	np.savetxt(folder + "velocity.txt", vels, delimiter='\t', fmt="%f")
	ids = np.arange(1, length + 1, 1)
	np.savetxt(folder + "id.txt", ids, delimiter='\t', fmt="%d")
	masses = np.empty(length)
	masses.fill(float(Mh) / N)
	np.savetxt(folder + "masses.txt", masses, delimiter='\t', fmt="%f")
	smooths = np.zeros(length)
	np.savetxt(folder + "smoothing.txt", smooths, delimiter='\t', fmt="%f")

def main():
	init()
	coords, radii = set_positions()
	print "done with the positions"
	vels = set_velocities(radii)
	print "done with the velocities"
	write_input_files(coords, vels)
	write_snapshot(folder)
	shutil.rmtree(folder)

if __name__ == '__main__':
	main()
