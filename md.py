import numpy as np
from poisson import poisson

# globals
KBT = 1


class Particle:
    def __init__(self, mass, q_init, v_init, charge, force_fun, k=None):
        self.k = k
        self.mass = mass
        self.q = q_init
        self.v = v_init
        self.force_fun = force_fun

        self.q_arr = []
        self.force_arr = []
        self.v_arr = []

        # NVT constants
        self.c1 = None
        self.c2 = None

        # for poisson eq.
        self.charge = charge

    @property
    def force(self):
        return self.force_fun(self)

    def update_v(self, dt):
        self.v = self.v+(self.force/self.mass)*dt
        
    def update_q(self, dt):
        self.q = self.q + dt*self.v
    
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


def verlet_integrator(particles: list, pes, n_steps, dt, rand_max=1, system='NVE', gamma=None, skip=0):
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
    i = 0
    while i < n_steps:
        # before iterating through all particles, run all updates and checks
        check_update(pes, particles)

        # iterate over particles in particle list
        for pi, p in enumerate(particles):
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

            p.update_v(dt/2.)

            if i % (skip + 1) == 0:
                # update arrays if this is not a step to skip
                p.update_histories()
        
        i += 1
        for p in particles:
            pass
    return t_arr

def check_update(pes, particles):
    pes.calculate(particles)
