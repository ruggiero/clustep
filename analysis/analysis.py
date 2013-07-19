from snapread import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from bisect import *

a = 200
Mh = 100000
N = 100000
G = 43007.1
time = 0.0

def hernquist_density(r):
	if(r == 0):
		return 0
	else:
		return (Mh * a) / (2 * np.pi * r * (r + a)**3)

def radial_speed_squared(r):
	return ((G * Mh) / (12 * a)) * ((12 * r * (r + a)**3 *\
				np.log((r + a) / r)) / a**4 - (r / (r + a)) * (25 + 52 * r / a\
				+ 42 * (r / a)**2 + 12 * (r / a)**3))

# Given a data vector, in which each element represents a different particle
# by a list of the form [radius, radial_speed^2], ordered according to the
# radii; and a multiplication factor, returns the right indexes of a log
# partition of the vector. Also returns an auxiliary vector, which will be
# useful in the functions that calculate the distrubution functions.
def log_partition(data, factor):
	limits = []
	auxiliary = []
	radii = [i[0] for i in data]
	left_limit = 0
	right_limit = 0.01
	left_index = 0
	while(right_limit < 200 * a):
		right_index = left_index + bisect_left(radii[left_index:], right_limit)
		limits.append(right_index)
		auxiliary.append([right_index - left_index, (right_limit + left_limit) /
			2])
		left_limit = right_limit
		left_index = right_index
		right_limit *= factor
	return limits, auxiliary

# Returns a list containing elements of the form [radius, density]
def density_distribution(data, partition, aux):
	distribution = []
	left = 0
	for j in range(len(partition)):
		right = partition[j]
		count = aux[j][0]
		middle_radius = aux[j][1]
		if(count > 0):
			density = 10**10 * (3.0 * count * Mh) / (4 * N * np.pi *
				(data[right][0]**3 - data[left][0]**3))
			distribution.append([middle_radius, density])
		else:
			distribution.append([middle_radius, 0])
		left = right
	return distribution

# Returns a list containing elements of the form [radius, density]
def radial_speed_distribution(data, partition, aux):
	distribution = []
	left = 0
	for j in range(len(partition)):
		right = partition[j]
		count = aux[j][0]
		middle_radius = aux[j][1]
		if(count > 0):
			sum_ = 0
			for i in range(left, right):
				sum_ += data[i][1]
			distribution.append([middle_radius, sum_ / count])
		else:
			distribution.append([middle_radius, 0])
		left = right
	return distribution

def density_plot(input_, data, part, aux):
	dist = density_distribution(data, part, aux)
	x_axis = np.logspace(np.log10(dist[0][0]), np.log10(dist[-1][0]),\
		num=1000)
	p1, = plt.plot([i[0] for i in dist], [i[1] for i in dist], 'o')
	p2, = plt.plot(x_axis, [10**10 * hernquist_density(i) for i in x_axis])
	plt.legend([p1, p2], ["Simulation", "Theoretical value"], loc=1)
	plt.xlabel("Radius (kpc)")
	plt.ylabel("Density (M${_\circ}$/kpc$^3$)")
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim([1, 10**4])
	plt.ylim([1, 10**10])
	plt.title("t = %1.2f Gyr" % time)
	plt.savefig(input_ + "-density.png")
	plt.close()
	print "Done with density for " + input_

def radial_speed_plot(input_, data, part, aux):
	dist = radial_speed_distribution(data, part, aux)
	x_axis = np.logspace(np.log10(dist[0][0]), np.log10(dist[-1][0]),\
		num=1000)
	formatter = FuncFormatter(lambda x, pos : "%1.2f" % (x / 10**6))
	ax = plt.subplot(111)
	ax.yaxis.set_major_formatter(formatter)

	# 0.9574 is a conversion factor from kpc^2/Gyr^2 to km^2/s^2
	p1, = plt.plot([i[0] for i in dist], [0.9574 * i[1] for i in dist], 'o')
	p2, = plt.plot(x_axis, 0.9574 * radial_speed_squared(x_axis))
	plt.legend([p1, p2], ["Simulation", "Theoretical value"], loc=1)
	plt.xlabel("Radius (kpc)")
	plt.ylabel("$v_{r}^{2}$ (10$^{6}$ km$^{2}$/s$^{2}$)")
	plt.xscale('log')
	plt.title("t = %1.2f Gyr" % time)
	#plt.gcf().subplots_adjust(left=0.17)
	#plt.yscale('log')
	plt.xlim([1, 10**4])
	plt.ylim([0, 2.5 * 10**6])
	plt.savefig(input_ + "-speed.png")
	plt.close()
	print "Done with speed for " + input_

def main():
	input_ = sys.argv[1]
	snapshot = open(input_, 'r')
	
	# Functions from snapread.py
	h = header(snapshot) 
	global time
	time = h.time
	p_list = read_data(snapshot, h)	
	snapshot.close()
	data = []
	for i in p_list:
		r = np.linalg.norm(i.pos)
		vr = np.dot(i.pos, i.vel) / r
		data.append([r, vr**2])
	data = sorted(data)
	part, aux = log_partition(data, 2)
	density_plot(input_, data, part, aux)
	radial_speed_plot(input_, data, part, aux)

if __name__ == '__main__':
	main()
