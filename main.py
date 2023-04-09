import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

from md import verlet_integrator, Particle
from poisson import PoissonPES
from membrane import Membrane, Material
from types import ParticleType

# GLOBALS
e = 1.602*10**(-19)

# NOTE: This is only for testing purposes
def SHO_force(k, q):
    return -k*q

# NOTE: This is a WIP to do animation of particle positions
'''
def anifunc(n):
    pp_x = [p.q_arr[n][0] for p in particles]
    pp_y = [p.q_arr[n][1] for p in particles]
    points.set_data(pp_x, pp_y)
    return points,
'''

def main():
    lx = 10

    ### STEP 1: setup membrane ###
    material = Material() # NOTE see membrane.py for defaults
    thickness = 3
    x_centre = lx/2 # centre the membrane in the simulation area
    height = 10
    n_slices = 10
    membrane = Membrane(thickness, x_centre, height, material, n_slices)

    ### STEP 2: setup ions in solution ###
    NC = 500  # number of charged catholyte ions
    mass_c = 82.94 # mass of catholyte ions
    charge_c = e # charge of catholyte ions
    NA = 500 # number of charged anolyte ions
    mass_a = 50.94 # mass of anolyte ions
    charge_a = 2*e# charge of anolyte ions
    NP_C = 500 # number of protons cathode side
    NP_A = 500 # number of protons anode side
    mass_p = 1 # mass of protons
    charge_p = e # charge of proton
    n_steps = 10**3

    # set up poisson object
    poisson_pes = PoissonPES(nx=100, lx=1, p0=0, pl=1)

    # randomize initial particle positions and velocities
    # NOTE: putting catholyte on the left and anolyte on the right
    particles = [Particle(mass_c, np.array([np.random.rand()*(x_centre-(thickness/2)), np.random.rand()]), 2*np.random.uniform(1,2)-1, charge_c, type=ParticleType.CATHODE_ION) for n in range(NC)]
    particles += [Particle(mass_a, np.array([np.random.rand()*(x_centre-(thickness/2))+((x_centre+(thickness/2))), np.random.rand()]), 2*np.random.uniform(1,2)-1, charge_a, type=ParticleType.ANODE_ION) for n in range(NA)]
    particles += [Particle(mass_p, np.array([np.random.rand()*(x_centre-(thickness/2)), np.random.rand()]), 2*np.random.uniform(1,2)-1, charge_p, type=ParticleType.PROTON) for n in range(NP_C)]
    particles += [Particle(mass_p, np.array([np.random.rand()*(x_centre-(thickness/2))+((x_centre+(thickness/2))), np.random.rand()]), 2*np.random.uniform(1,2)-1, charge_p, type=ParticleType.PROTON) for n in range(NP_A)]
    
    verlet_integrator(particles, poisson_pes, n_steps, 0.1, system='NVT', gamma=10)

    # NOTE: This is a WIP to do animation of particle positions
    '''
    global points
    points, = ax.plot([], [], 'o')
    ani = animation.FuncAnimation(fig, anifunc, n_steps)
    plt.show()
    '''

if __name__ == '__main__':
    main()
