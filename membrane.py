import numpy as np
from dataclasses import dataclass

@dataclass
class Material:
    pore_charge = 1
    base_charge = 1

    pore_density = 1  # number density of particles in pores (related to hydration)
    base_density = 1  # number density of particles in backbone
    cf = 1  # continuity factor (see pore_locations)

class Membrane:
    def __init__(self, porosity, psd: dict, thickness, x_centre, height, material, n_slices):
        self.porosity = porosity
        self.x_l = x_centre -(thickness/2)
        self.x_r = x_centre + (thickness/2)

        dx, slices, pores = pore_locations(psd, porosity, x_centre, height, thickness, n_slices, material.cf)
        self.base_particles, self.pore_particles = particle_locations(slices, dx, pores, height, material.pore_density, material.base_density)


def merge_intervals(intervals):
    # Sort the array on the basis of start values of intervals.
    intervals = intervals.sorted(key=lambda tup: tup[0])
    stack = []
    # insert first interval into stack
    stack.append(intervals[0])
    for i in intervals[1:]:
        # Check for overlapping interval,
        # if interval overlap
        if stack[-1][0] <= i[0] <= stack[-1][-1]:
            stack[-1][-1] = max(stack[-1][-1], i[-1])
        else:
            stack.append(i)
    
    return stack


def pore_locations(psd, porosity, x_centre, height, thickness, n_slices, cf):
    # setup membrane pores using the porosity and pore
    # size distribution
    x_l = x_centre -(thickness/2)
    x_r = x_centre + (thickness/2)

    # number of pores is the porosity * membrane length (total amount of desired empty space)
    # divided by the mean pore size
    n_pores = (porosity*height) / psd["loc"]
    
    dx = thickness / n_slices

    # make array of x value ranges
    # array format is [(min 1, max 1), (min 2, max2), ...]
    slices = np.array(range(n_slices)) * dx
    slice_ranges = np.array([(slices[n], slices[n+1]) for n in range(n_slices)])
    
    # list of sets of tuples (LMAO), each set represents one slice in the 
    # x direction.  Each tuple in each set is an interval representing the
    # boundaries of one pore in the form (pore bottom y, pore top y)
    pores = [set()] * n_slices
    
    # iterate through each pore, and within each pore iterate through 
    # thickness, placing new pores offset depending on the continuity factor (cf)
    for pi in range(n_pores):
        # place first pore centre randomly
        m = np.random.uniform(low=0, high=height)
        # get pore size from normal pore size distribution
        d = np.random.normal(psd)

        # place first top and bottom tuples in corresponding set
        pores[0].add((m-(d/2), m+(d/2)))
        for si in range(1, n_slices):
            # move middle of pore depending on continuity factor
            m += cf*np.random.uniform(low=-1, high=1)
            # place top and bottom in appropriate slice
            pores[si].add((m-(d/2), m+(d/2)))
    
    # merge overlapping pores in each slice.  This is important for the next step
    # note that what was previously a set is now a list
    for si in range(0, slices):
        pores[si] = merge_intervals(pores[si])
    
    return dx, slices, pores

def binary_search_intervals(intervals, point):
    # find the adjecent interval minimums that the point is in between
    low = 0
    high = len(intervals)
    while (abs(high-low)) > 1:
        s = (low-high) // 2
        # check which half the point is in
        if intervals[s][0] <= point:
            low = s
        else:
            high = s
    
    # found the interval it could be in, check if it is actually in that interval
    inter = intervals[low]
    return (point < inter[1])

def particle_locations(slices, dx, pores, ly, pore_density, base_density):
    ''' Places membrane material particles in base and pore regions
    '''
    # lists to hold locations of each particle.  Again, each list entry is a tuple
    # in the form (x, y)
    base_particles = []
    pore_particles = []


    # find pore areas to determine number of pore and base particles to place
    p_area = 0
    for si in range(1, len(slices)+1):
        # get slice range
        slice = (slices(si-1), slices(si))
        for pore in pores:
            # calculate area of pore
            area = (slice[1] - slice[0]) * (pore[1] - pore[0])
            p_area += area
    # determine overall area
    t_area = (slices[-1] - slices[0]) * ly
    # determine base area
    b_area = t_area - p_area

    # determine total number of each kind of particle to place
    total_bp = int(np.round(base_density * b_area))
    total_pp = int(np.round(pore_density * p_area))

    n_bp = 0
    n_pp = 0
    while (n_bp < total_bp) or (n_pp < total_pp):
        # generate random point in the membrane area
        point = (np.random.uniform(low=slices[0], high=slices[-1]), np.random.uniform(low=0, high=ly))
        # determine the slice that this point is in
        si = point[0] // dx
        # check if the y value places it in a pore with a very simple binary search. (could do this faster)
        if binary_search_intervals(pores[si], point[1]):
            # place a pore particle
            if n_pp < total_pp:
                pore_particles.append(point)
                n_pp += 1
        else:
            # place a base particle
            if n_bp < total_bp:
                base_particles.append(point)
                n_bp += 1
    
    return base_particles, pore_particles
