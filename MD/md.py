import numpy as np
from matplotlib import pyplot as plt

from MD.forces import pair_pot
from MD.types import ParticleType

# globals
KBT = 1


class Particle:
    def __init__(self, mass, q_init, v_init, charge, x_range, y_range, type: ParticleType, k=None):
        self.k = k
        self.mass = mass
        self.q = q_init
        self.v = v_init

        self.q_arr = []
        self.force_arr = []
        self.v_arr = []

        self.x_range = x_range
        self.y_range = y_range

        # NVT constants
        self.c1 = None
        self.c2 = None

        # for forces
        self.charge = charge
        self.type = type
        self._force = np.zeros(2)

    @property
    def force(self):
        return self._force
    
    def update_force(self, particles, membrane):
        # set the force variable using membrance
        # force field, and pair potential force field
        membrane_force = membrane.calculate_interactions(self)
        pair_pot_force = pair_pot(self, particles)
        self._force = membrane_force + pair_pot_force

    def update_v(self, dt):
        self.v = self.v+(self.force/self.mass)*dt
        
    def update_q(self, dt):
        # check elastic boundaries and update v if necessary
        q = self.q + dt*self.v
        v_mod = elastic_boundaries(q, self.x_range, self.y_range)
        if np.any(v_mod < 1):
            self.v = np.multiply(self.v, v_mod)
            q = self.q + dt*self.v 
        self.q = q
    
    def update_histories(self):
        self.q_arr.append(self.q)
        self.v_arr.append(self.v)
        self.force_arr.append(self.force_arr)


def mag_direction(vec):
    mag = np.sqrt(vec[0]**2 + vec[1]**2)
    dir = (vec[0]/mag, vec[1]/mag)
    return mag, dir


def NVT_constants(gamma, dt, mass):
    c1 = np.exp(-gamma*dt)
    c2 = np.sqrt((1-np.exp(-2*gamma*dt))*KBT/mass)
    return c1, c2


def verlet_integrator(particles: list, pes, membrane, n_steps, dt, rand_max=1, system='NVE', gamma=None, skip=0):
    t_arr = [i*dt for i in range(n_steps+1)]
    n_par = len(particles)
    n_d = len(particles[0].q)

    if system == 'NVT':
        # set NVT constants for each particle
        for p in particles:
            p.c1, p.c2 = NVT_constants(gamma, dt, p.mass)

        # NOTE: I calculate all the random numbers needed ahead of time for speed
        rands = np.random.normal(size=(n_par, n_steps, n_d))

    # iterate over time steps (keep track of step in case steps are skipped)
    # show initial positions
    scatter(particles, membrane)
    i = 0
    while i < n_steps:
        # iterate over particles in particle list
        for pi, p in enumerate(particles):
            p.update_force(particles, membrane)
            p.update_v(dt/2.)
            if system == 'NVE':
                p.update_q(dt)
            elif system == 'NVT':
                p.update_q(dt/2.)

                # apply random motion
                rand = rands[pi][i][:]
                # update x and y velocities (both random)
                p.v = (p.c1*p.v + rand*p.c2)

                p.update_q(dt/2.)

            p.update_force(particles, membrane)
            p.update_v(dt/2.)

            if i % (skip + 1) == 0:
                # update arrays if this is not a step to skip
                p.update_histories()

        if i % 10 == 0:
            scatter(particles, membrane, i)
            pes.calculate(particles)
            plt.plot(pes.x_arr, pes.v_func(pes.x_arr))
            plt.title(f"Potential Distribution, Iteration {i}")
            plt.show()
        
        i += 1
    return t_arr

def scatter(particles, membrane, i):
    plt.scatter([p[0] for p in membrane.base_particles], [p[1] for p in membrane.base_particles], label="Membrane - Base")
    plt.scatter([p[0] for p in membrane.pore_particles], [p[1] for p in membrane.pore_particles], label="Membrane - Pores")
    plt.xlim(0, 10)
    plt.ylim(0,10)
    plt.scatter([p.q[0] for p in particles if p.type is ParticleType.CATHODE_ION], [p.q[1] for p in particles if p.type is ParticleType.CATHODE_ION], label="Cathode Ions")
    plt.scatter([p.q[0] for p in particles if p.type is ParticleType.ANODE_ION], [p.q[1] for p in particles if p.type is ParticleType.ANODE_ION], label="Anode Ions")
    plt.scatter([p.q[0] for p in particles if p.type is ParticleType.PROTON], [p.q[1] for p in particles if p.type is ParticleType.PROTON], label="Protons")
    plt.legend()
    plt.title(f"Particle Positions, Iteration {i}")
    plt.show()

def elastic_boundaries(q, x_range, y_range):
    v_mod = np.array([1, 1])
    if (q[0] > x_range[1]) or (q[0] < x_range[0]):
        v_mod[0] = -1
    
    if (q[1] > y_range[1]) or (q[0] < y_range[1]):
        v_mod[1] = -1
    
    return v_mod
