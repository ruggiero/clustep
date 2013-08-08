import numpy as np

Mh = 100000
G = 43007.1

# calculates the numerical potential at point (which is given as a
# NumPy array), given a vector containing the radii of the considered
# particles, and a dictionary where values of potential that have been
# previously calculated are stored
def potential(point, radii, recorded_phis):
    id_ = hash(tuple(point)) # arrays are not hashable, but tuples are
    if(id_ in recorded_phis):
        return recorded_phis[id_]
    else:
        distances = np.array([np.linalg.norm(x - point) for x in radii])
        s = np.sum(1. / distances)
        phi = -G * Mh * s
        recorded_phis[id_] = phi
        return phi

# The displacements that will be added or subtracted to the current
# guess for the location of the COD
def possible_sums(step):
    return [np.array((step, 0, 0)), np.array((0, step, 0)),
            np.array((0, 0, step))]

# Previous knowledge on the distance between the COM and the CED
# is necessary for estimating good values for the step
def COD(p_list):
    COM = sum(p.pos for p in p_list) / len(p_list)
    print "the COM is (%f, %f, %f)" % (COM[0], COM[1], COM[2])
    recorded_phis = {}
    step = 64 # kpc
    psums = possible_sums(step)
    guess = (0, 0, 0)
    new_guess = COM
    close_radii = [] # distant particles won't be considered
    for p in p_list:
        if(np.linalg.norm(p.pos - COM) < 1.25 * 64): # this is arbitrary
            close_radii.append(p.pos)
    while(step >= 4):
        print "finding COD... step = %d kpc" % step

        # if the best guess for a given step has been found, reduce
        # the size of the step by half
        if(np.linalg.norm(new_guess - guess) == 0):
            step /= 2
            psums = possible_sums(step)
        guess = np.array(new_guess)
        g = [guess]
        for s in psums:
            g.append(guess + s)
            g.append(guess - s)
        best = min([[potential(i, close_radii, recorded_phis), i]
                    for i in g])[1]
        new_guess = np.array(best)
        if(step < 4):
            print "the COD is (%f, %f, %f)" % (new_guess[0], new_guess[1],
                                               new_guess[2])
    return new_guess
