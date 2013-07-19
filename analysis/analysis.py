from snapread import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

a = 200
Mh = 100000
N = 100000
G = 43007.1

def hernquist_density(r):
	return (Mh * a) / (2 * np.pi * (r + 0.000001) * (r + a)**3)

def radial_speed_squared(r):
	return ((G * Mh) / (12 * a)) * ((12 * r * (r + a)**3 * \
		   np.log((r + a) / r)) / a**4 - (r / (r + a)) * (25 + 52 * r / a +\
		   42 * (r / a)**2 + 12 * (r / a)**3))

def density_distribution(p_list):
	radii = []
	distribution = []
	for i in p_list:
		r = np.linalg.norm(i.pos)
		radii.append(r)
	radii.sort()
	left_limit = 0
	right_limit = 0.01
	left_index = 0
	while(right_limit < 200 * a):

		# This is a esoteric, yet fast, way to find the index for the first element in
		# the vector x which is greater than right_limit
		right_index = next(x[0] for x in enumerate(radii[left_index:]) if x[1] > right_limit)
		count = right_index - left_index
		density = 10**10 * (count * Mh / N) / (4 / 3. * np.pi * (right_limit**3 - left_limit**3))
		distribution.append([(10 * right_limit + left_limit) / 11, density])
		left_limit = right_limit
		left_index = right_index
		right_limit *= 2
	return distribution

def radial_speed_distribution(p_list):

	# vector which will hold vectors of the form [radius, vr^2]
	pairs = []
	distribution = []
	for i in p_list:
		r = np.linalg.norm(i.pos)
		vr = np.dot(i.pos, i.vel) / r
		pairs.append([r, vr**2])
	pairs = sorted(pairs)
	left_limit = 0
	right_limit = 0.01
	left_index = 0
	while(right_limit < 200 * a):
		sum_ = 0
		right_index = next(x[0] for x in enumerate(pairs[left_index:]) if x[1][0] > right_limit)
		count = right_index - left_index
		if(count > 0):
			for i in range(left_index, right_index):
				sum_ += pairs[i][1]
			distribution.append([(10 * right_limit + left_limit) / 11, sum_ / count])
		else:
			distribution.append([(10 * right_limit + left_limit) / 11, 0])
		left_limit = right_limit
		left_index = right_index
		right_limit *= 2
	return distribution

def main():
	input_ = sys.argv[1]
	snapshot = open(input_, 'r')
	
	# Functions from snapread.py
	h = header(snapshot) 
	p_list = read_data(snapshot, h)	
	snapshot.close()

	# Dealing with the density distribution plot
	x = density_distribution(p_list)
	x_axis = np.logspace(np.log10(x[0][0]), np.log10(x[len(x) - 1][0]), num=1000)
	p1, = plt.plot([i[0] for i in x], [i[1] for i in x], 'o')
	p2, = plt.plot(x_axis, [10**10 * hernquist_density(i) for i in x_axis])
	plt.legend([p1, p2], ["Simulation", "Theoretical value"], loc=1)
	plt.xlabel("Radius (kpc)")
	plt.ylabel("Density (M${_\circ}$/kpc$^3$)")
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim([1, 10**4])
	plt.ylim([1, 10**10])
	plt.title("t = %1.2f Gyr" % h.time)
	plt.savefig(input_ + "-density.png")
	plt.close()
	print "Done with density for " + input_

	# Now the radial speed distribution plot
	x = radial_speed_distribution(p_list)
	x_axis = np.logspace(np.log10(x[0][0]), np.log10(x[len(x) - 1][0]), num=1000)
	formatter = FuncFormatter(lambda x, pos : "%1.2f" % (x / 10**6))
	ax = plt.subplot(111)
	ax.yaxis.set_major_formatter(formatter)

	# 0.9574 is a conversion factor from kpc^2/Gyr^2 to km^2/s^2
	p1, = plt.plot([i[0] for i in x], [0.9574 * i[1] for i in x], 'o')
	p2, = plt.plot(x_axis, 0.9574 * radial_speed_squared(x_axis))
	plt.legend([p1, p2], ["Simulation", "Theoretical value"], loc=1)
	plt.xlabel("Radius (kpc)")
	plt.ylabel("$v_{r}^{2}$ (10$^{6}$ km$^{2}$/s$^{2}$)")
	plt.xscale('log')
	plt.title("t = %1.2f Gyr" % h.time)
	#plt.gcf().subplots_adjust(left=0.17)
	#plt.yscale('log')
	plt.xlim([1, 10**4])
	plt.ylim([0, 2.5 * 10**6])
	plt.savefig(input_ + "-speed.png")
	print "Done with speed for " + input_

if __name__ == '__main__':
	main()
