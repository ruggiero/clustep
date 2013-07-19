'''
Execution: ./snapread.py SNAPSHOT
'''

import sys
import numpy as np

# For reading and storing the header
class header:
	def __init__(self, snapshot):
		read_dummy(snapshot, 1)
		self.n_part = np.fromfile(snapshot, 'int32', 6)
		self.mass = np.fromfile(snapshot, 'float64', 6)
		self.time = np.fromfile(snapshot, 'float64', 1)
		self.redshift = np.fromfile(snapshot, 'float64', 1)
		self.flag_sfr = np.fromfile(snapshot, 'int32', 1)
		self.flag_feedback = np.fromfile(snapshot, 'int32', 1)
		self.n_part_total = np.fromfile(snapshot, 'int32', 6)
		self.flag_cooling = np.fromfile(snapshot, 'int32', 1)
		self.num_files = np.fromfile(snapshot, 'int32', 1)
		self.box_size = np.fromfile(snapshot, 'float64', 1)
		self.omega0 = np.fromfile(snapshot, 'float64', 1)
		self.omega_lambda = np.fromfile(snapshot, 'float64', 1)
		self.hubble_param = np.fromfile(snapshot, 'float64', 1)
		self.fill = np.fromfile(snapshot, 'int8', 96)
		read_dummy(snapshot, 1)

# These declared variables are for clarity purposes
class particle:
	pos = None
	vel = None
	mass = 0.0
	U = None
	rho = None
	smoothing = None

# For reading the ints that are present in the beginning and in the end
# of each block
def read_dummy(snapshot, n_dummies):
	for i in range(n_dummies):
		dummy = np.fromfile(snapshot, 'int32', 1)

def read_data(snapshot, h):
	p_list = []
	n_part = sum(h.n_part)
	read_dummy(snapshot, 1)

	# positions
	for i in range(n_part):
		p = particle()
		p.pos = np.fromfile(snapshot, 'float32', 3)
		p_list.append(p)
	read_dummy(snapshot, 2)
	
	# velocities
	for i in range(n_part):
		p_list[i].vel = np.fromfile(snapshot, 'float32', 3)
	read_dummy(snapshot, 2)

	# IDs
	for i in range(n_part):
		p_list[i].ID = np.fromfile(snapshot, 'int32', 1)
	read_dummy(snapshot, 2)

	# Variable masses, which are read in case the mass of the
	# particle of type 'i' is declared as 0, in the header
	cur = 0
	read_something = 0
	for i in range(6):
		if(h.mass[i] != 0):
			cur += h.n_part[i]
			continue
		else:
			read_something = 1
			for j in range(h.n_part[i]):
				p_list[cur].mass = np.fromfile(snapshot, 'float32', 1)
				cur += 1

	# The variable masses block might not exist
	if(read_something):
		read_dummy(snapshot, 2)
	
	# Blocks related to the internal energies, densities and smoothing
	# lengths of the gas particles, in case there is any
	if(h.n_part[0] > 0):

		# First the energies
		for i in range(h.n_part[0]):
			p_list[i].U = np.fromfile(snapshot, 'float32', 1)
		read_dummy(snapshot, 2)

		# Then the densities
		for i in range(h.n_part[0]):
			p_list[i].rho = np.fromfile(snapshot, 'float32', 1)
		read_dummy(snapshot, 2)

		# And the smoothing lengths
		for i in range(h.n_part[0]):
			p_list[i].smoothing = np.fromfile(snapshot, 'float32', 1)
		read_dummy(snapshot, 2)
	
	# There are some blocks to be implemented...
	chunk = snapshot.read()
	if chunk:
		print "There still were things to be read..."
	return p_list

#def main():
#	snapshot = open(sys.argv[1], 'r')
#	h = header(snapshot)
#	print h.n_part
#	print h.mass
#	print h.time
#	print h.redshift
#	print h.flag_sfr
#	print h.flag_feedback
#	print h.n_part_total
#	print h.flag_cooling
#	print h.num_files
#	print h.box_size
#	print h.omega0
#	print h.omega_lambda
#	print h.hubble_param
#	print h.fill
#	p_list = read_data(snapshot, h)
#	snapshot.close()

	# From here on, the data for all the particles is accessible
	# As an example, showing the positions of all the particles
#	for i in p_list:
		#print i.vel[0], i.vel[1], i.vel[2]
#		print i.pos[0], i.pos[1], i.pos[2]
		#print i.ID
		#print i.rho
		#print i.smoothing

#if __name__ == '__main__':
#	main()
