from snapread import *
import numpy as np
import matplotlib
matplotlib.use('Agg') # to be able to plot under an SSH session
import matplotlib.pyplot as plt
#from matplotlib.ticker import FuncFormatter
from bisect import *

a = 200
Mh = 100000
G = 43007.1
time = 0.0

def hernquist_density(r):
    if(r == 0):
        return 0
    else:
        return (Mh * a) / (2 * np.pi * r * (r + a)**3)

def radial_velocity(r):
    cte = (G * Mh) / (12 * a)
    t1 = ((12 * r * (r + a)**3) / a**4) * np.log((r + a) / r)
    t2 = -r/(r + a) * (25 + 52 * r / a + 42 * (r / a)**2 + 12 * (r / a)**3)
    return (cte * (t1 + t2))**0.5

# Given a data vector, in which each element represents a different
# particle by a list of the form [radius, radial_velocity^2], ordered
# according to the radii; and a multiplication factor, returns the right
# indexes of a log partition of the vector. Also returns an auxiliary
# vector, which will be useful in the functions that calculate the
# distribution functions.
def log_partition(data, factor):
    limits = []
    auxiliary = []
    radii = [i[0] for i in data]
    left_limit = 0
    right_limit = 0.01
    left_index = 0
    while(right_limit < 200 * a):

        # before right_index, everybody is smaller than right_limit
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
    cte = (10**10 * 3 * Mh) / (4 * np.pi * N)
    for j in np.arange(len(partition)):
        right = partition[j]
        count = aux[j][0]
        middle_radius = aux[j][1]
        if(count > 0):
            density = (cte * count) / (data[right][0]**3 - data[left][0]**3)
            distribution.append([middle_radius, density])
        else:
            distribution.append([middle_radius, 0])
        left = right
    return distribution

# Returns a list containing elements of the form [radius, vr]
def radial_velocity_distribution(data, partition, aux):
    distribution = []
    left = 0
    for j in np.arange(len(partition)):
        right = partition[j]
        count = aux[j][0]
        middle_radius = aux[j][1]
        if(count > 0):
            sum_ = sum(i[1] for i in data[left:right])
            distribution.append([middle_radius, (sum_ / count)**0.5])
        else:
            distribution.append([middle_radius, 0])
        left = right
    return distribution

def density_plot(input_, data, part, aux):
    dist = density_distribution(data, part, aux)
    x_axis = np.logspace(np.log10(dist[0][0]), np.log10(dist[-1][0]), num=500)
    p1, = plt.plot([i[0] for i in dist], [i[1] for i in dist], 'o')
    p2, = plt.plot(x_axis, [10**10 * hernquist_density(i) for i in x_axis])
    plt.legend([p1, p2], ["Simulation", "Theoretical value"], loc=1)
    plt.xlabel("Radius (kpc)")
    plt.ylabel("Density (M$_{\odot}$/kpc$^3$)")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1, 10**4])
    plt.ylim([1, 10**10])
    plt.title("t = %1.2f Gyr" % time)
    plt.savefig(input_ + "-density.png")
    plt.close()
    print "Done with density for " + input_

def radial_velocity_plot(input_, data, part, aux):
    dist = radial_velocity_distribution(data, part, aux)
    x_axis = np.logspace(np.log10(dist[0][0]), np.log10(dist[-1][0]), num=500)
    #formatter = FuncFormatter(lambda x, pos : "%1.2f" % (x / 10**6))
    #ax = plt.subplot(111)
    #ax.yaxis.set_major_formatter(formatter)
    p1, = plt.plot([i[0] for i in dist], [i[1] for i in dist], 'o')
    p2, = plt.plot(x_axis, radial_velocity(x_axis))
    plt.legend([p1, p2], ["Simulation", "Theoretical value"], loc=1)
    plt.xlabel("Radius (kpc)")
    plt.ylabel("$\sqrt{\overline{v_r^2}}$ (km/s)")
    plt.xscale('log')
    plt.title("t = %1.2f Gyr" % time)
    #plt.gcf().subplots_adjust(left=0.17)
    #plt.yscale('log')
    plt.xlim([1, 10**4])
    #plt.ylim([0, 2.5 * 10**6])
    plt.ylim([0, 1650])
    plt.savefig(input_ + "-velocity.png")
    plt.close()
    print "Done with velocity for " + input_

def main():
    global time, N
    input_ = sys.argv[1]
    snapshot = open(input_, 'r')
    h = header(snapshot) # From snapread.py
    time = h.time
    N = sum(h.n_part_total)
    p_list = read_data(snapshot, h) # From snapread.py
    snapshot.close()
    data = []
    for i in p_list:
        r = np.linalg.norm(i.pos)
        vr = np.dot(i.pos, i.vel) / r
        data.append([r, vr**2])
    del(p_list)
    data = sorted(data)
    part, aux = log_partition(data, 1.3)
    radial_velocity_plot(input_, data, part, aux)
    density_plot(input_, data, part, aux)

if __name__ == '__main__':
    main()
