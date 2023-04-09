import numpy as np

from matplotlib import pyplot as plt 
from dataclasses import dataclass

from forces import LJF_attractive, LJF_repulsive, coulomb
from types import ParticleType

@dataclass
class Material:
    porosity = 0.4
    psd = {'loc':1, 'scale':0.01}
    pore_charge = 1
    # lennard jones parameters for each ion type interacting with pore
    pore_eps = {ParticleType.PROTON: 1
                ParticleType.CATHODE_ION: 1
                ParticleType.ANODE_ION: 1}
    pore_sigma = {ParticleType.PROTON: 1
                  ParticleType.CATHODE_ION: 1
                  ParticleType.ANODE_ION: 1}

    # lennard jones parameters for each ion type interacting with base
    # NOTE: base is assumed to be neutral and has no charge
    base_eps = {ParticleType.PROTON: 1
                ParticleType.CATHODE_ION: 1
                ParticleType.ANODE_ION: 1}
    base_sigma = {ParticleType.PROTON: 1
                  ParticleType.CATHODE_ION: 1
                  ParticleType.ANODE_ION: 1}

    pore_density = 100  # number density of particles in pores (related to hydration)
    base_density = 10000 # number density of particles in backbone
    cf = 1  # continuity factor (see pore_locations)

class Membrane:
    def __init__(self, thickness, x_centre, height, material, n_slices):
        self.porosity = porosity
        self.x_l = x_centre -(thickness/2)
        self.x_r = x_centre + (thickness/2)

        self.material = material
        psd = self.material.psd
        porosity = self.material.porosity
        self.base_particles, self.pore_particles = place_particles(psd, porosity, x_centre, height, thickness, n_slices, material.cf, material.base_density, material.pore_density)

    def calculate_interactions(self, particle):
        '''
        Calculates total force from membrane particles on a particle
        particle -> particle object
        '''
        force_total = np.zeros(2)
        for p in self.base_particles:
            # use lennard jones attractive and repulsive potentials bc
            # base particles are uncharged??
            # NOTE: not sure if this is right?
            eps = self.material.base_eps[p.type]
            sigma = self.material.base_sigma[p.type]
            force = LJF_attractive(p, particle.q, eps, sigma) + LJF_repulsive(p, particle.q, eps, sigma)
            force_total += force
        for p in self.pore_particles:
            # use coulomb attractive/repulsive and lennard jones repulsive
            # lennard jones repulsive will have almost no impact if the 
            # coulomb force is already repulsive, however if it is attractive
            # it will prevent particles from coming too close to each other
            eps = self.material.pore_eps[p.type]
            sigma = self.material.pore_sigma[p.type]
            force = coulomb(p, particle.q, self.material.pore_charge, p.charge) + LJF_repulsive(p, particle.q, eps, sigma)
            force_total += force
        
        return force_total

def merge_intervals(intervals):
    # Sort the array on the basis of start values of intervals.
    intervals.sort(key=lambda tup: tup[0])
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


def place_particles(psd, porosity, x_centre, height, thickness, n_slices, cf, base_density, pore_density):
    # setup membrane pores using the porosity and pore
    # size distribution
    x_l = x_centre -(thickness/2)
    x_r = x_centre + (thickness/2)

    # number of pores is the porosity * membrane length (total amount of desired empty space)
    # divided by the mean pore size
    n_pores = int(round((porosity*height) / psd["loc"]))

    dx = thickness / n_slices

    # make array of x value ranges
    slices = np.array(range(n_slices)) * dx
    
    # determine mean distance between pores
    dist_mean = (height - n_pores*psd["loc"])/n_pores
    # distance between pores distribution
    dd = {"loc":dist_mean, "scale": psd["scale"]}

    # iterate through the slices, generating pores
    head = 0
    base_particles = []
    pore_particles = []

    # first generate random pore locations in the first slice
    x_range = (slices[0], slices[1])
    pore_centres = []
    pore_sizes = []
    while True:
        # generate a non-porous region
        head, points = generate_region(head, dd, height, base_density, x_range)
        base_particles = base_particles + points
        if head >= height:
            break
        # generate a porous region
        btm = head
        head, points = generate_region(head, psd, height, pore_density, x_range)
        pore_particles = pore_particles + points
        pore_centres.append((head+btm)/2)
        pore_sizes.append((head-btm))
        if head >= height:
            break
    plt.scatter([p[0] for p in base_particles], [p[1] for p in base_particles], label="Base Particles")
    plt.scatter([p[0] for p in pore_particles], [p[1] for p in pore_particles], label="Pore Particles")
    plt.show()
    # now generate pore and base locations in the remaining slices
    for si in range(1, n_slices-1):
        x_range = (slices[si], slices[si+1])
        # update pore centres
        pore_centres = [c + np.random.uniform(low=-(cf/2), high=(cf/2)) for c in pore_centres]
        # generate intervals based on pore locations
        intervals = []
        for pi, c in enumerate(pore_centres):
            intervals.append([c - pore_sizes[pi]/2, c + pore_sizes[pi]/2])
        
        # merge overlapping intervals
        intervals = merge_intervals(intervals)
        head = 0
        for ind, invl in enumerate(intervals):
            # non- porous region
            y_range = (head, invl[0])
            # determin number of points
            b_area = (y_range[1] - y_range[0]) * (x_range[1]-x_range[0])
            n_bpts = int(np.round(base_density * b_area))
            points_b = generate_points(x_range, y_range, n_bpts)
            base_particles += points_b

            # porous region
            y_range = invl
            # determin number of points
            p_area = (y_range[1] - y_range[0]) * (x_range[1]-x_range[0])
            n_ppts = int(np.round(pore_density * p_area))
            points_p = generate_points(x_range, y_range, n_ppts)
            pore_particles += points_p
            head = invl[1]
        
        if head < height:
            # still one more non-porous region to generate
            y_range = (head, height)
            # determin number of points
            b_area = (y_range[1] - y_range[0]) * (x_range[1]-x_range[0])
            n_bpts = int(np.round(base_density * b_area))
            points_b = generate_points(x_range, y_range, n_bpts)
            base_particles += points_b
    
    return base_particles, pore_particles

def generate_region(head, dist, height, density, x_range):
    # generate a non-porous region
    b_n = head
    t_n = head + np.random.normal(**dist)

    # check if t_np is out of range
    if t_n > height:
        t_n = height

    head = t_n
    # figure out area of non-porous region and number of 
    # particles in the area
    area = (t_n - b_n) * (x_range[1]-x_range[0])
    n_p = int(np.round(density * area))

    # generate n_p random points in the non porous region
    y_range = (b_n, t_n)
    particles = (generate_points(x_range, y_range, n_p))

    return head, particles

def generate_points(x_range, y_range, npoints):
    points = []
    for i in range(npoints):
        point_x = np.random.uniform(low=x_range[0], high=x_range[1])
        point_y = np.random.uniform(low=y_range[0], high=y_range[1])
        points.append((point_x, point_y))

    return points

if __name__ == '__main__':
    material = Material()
    psd = {'loc':1, 'scale':0.01}

    membrane = Membrane(0.4, psd, 0.1, 1, 10, material, 30)
    plt.scatter([p[0] for p in membrane.base_particles], [p[1] for p in membrane.base_particles], label="Base Particles")
    plt.scatter([p[0] for p in membrane.pore_particles], [p[1] for p in membrane.pore_particles], label="Pore Particles")
    plt.show()
