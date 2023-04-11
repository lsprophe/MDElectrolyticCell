import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

from MD.md import verlet_integrator, Particle
from MD.poisson import PoissonPES
from MD.membrane import Membrane, Material
from MD.types import ParticleType

# GLOBALS
e = 1.602*10**(-19)
AVO = 6.022*10**(23)

# NOTE: This is a WIP to do animation of particle positions
'''
def anifunc(n):
    pp_x = [p.q_arr[n][0] for p in particles]
    pp_y = [p.q_arr[n][1] for p in particles]
    points.set_data(pp_x, pp_y)
    return points,
'''

def main():
    lx = 100  # nm
    x_range = np.array([0, lx])

    ### STEP 1: setup membrane ###
    material = Material(pore_density=50) # NOTE see membrane.py for defaults
    thickness = 10
    x_centre = lx/2 # centre the membrane in the simulation area
    height = 30 # nm
    y_range = np.array([0, height])
    n_slices = 50
    membrane = Membrane(thickness, x_centre, height, material, n_slices)

    ### STEP 2: setup ions in solution ###
    NC = 10  # number of charged catholyte ions
    mass_c = 1000 * 82.94/AVO # mass of catholyte ions
    charge_c = e # charge of catholyte ions
    NA = 10 # number of charged anolyte ions
    mass_a = 1000 * 50.94/AVO # mass of anolyte ions
    charge_a = 2*e# charge of anolyte ions
    NP_C = 5 # number of protons cathode side
    NP_A = 5 # number of protons anode side
    mass_p = 1000 * 1/AVO # mass of protons
    charge_p = e # charge of proton
    NS_C = 5 # number of sulfates cathode side
    NS_A = 5 # number of sulfates anode side
    mass_so4 = 1000 * 96.06/AVO # mass of sulfate
    charge_so4 = -2*e # charge of sulfate
    n_steps = 100

    ### STEP 3: Initialize poisson object to get potential results ###
    poisson_pes = PoissonPES(nx=100, lx=lx, p0=0, pl=1)

    ### STEP 4: Initialize ions in solution ###
    # randomize initial particle positions and velocities
    # NOTE: putting catholyte on the left and anolyte on the right
    particles = [Particle(mass_c, np.array([np.random.rand()*(x_centre-(thickness/2)), height*np.random.rand()]), 2*np.random.uniform(1,2)-1, charge_c, x_range, y_range, type=ParticleType.CATHODE_ION) for n in range(NC)]
    particles += [Particle(mass_a, np.array([np.random.rand()*(x_centre-(thickness/2))+((x_centre+(thickness/2))), height*np.random.rand()]), 2*np.random.uniform(1,2)-1, charge_a, x_range, y_range, type=ParticleType.ANODE_ION) for n in range(NA)]
    particles += [Particle(mass_p, np.array([np.random.rand()*(x_centre-(thickness/2)), height*np.random.rand()]), 2*np.random.uniform(1,2)-1, charge_p, x_range, y_range, type=ParticleType.PROTON) for n in range(NP_C)]
    particles += [Particle(mass_p, np.array([np.random.rand()*(x_centre-(thickness/2))+((x_centre+(thickness/2))), height*np.random.rand()]), 2*np.random.uniform(1,2)-1, charge_so4, x_range, y_range, type=ParticleType.PROTON) for n in range(NP_A)]
    particles += [Particle(mass_so4, np.array([np.random.rand()*(x_centre-(thickness/2)), height*np.random.rand()]), 2*np.random.uniform(1,2)-1, charge_p, x_range, y_range, type=ParticleType.SULFATE) for n in range(NS_C)]
    particles += [Particle(mass_so4, np.array([np.random.rand()*(x_centre-(thickness/2))+((x_centre+(thickness/2))), height*np.random.rand()]), 2*np.random.uniform(1,2)-1, charge_so4, x_range, y_range, type=ParticleType.SULFATE) for n in range(NS_A)]
    
    ### STEP 5: initialize membrane geometry ###
    membrane = Membrane(thickness, x_centre, height, material, n_slices)
    verlet_integrator(particles, poisson_pes, membrane, n_steps, 5e-8, system='NVT', gamma=10)

    # NOTE: This is a WIP to do animation of particle positions
    '''
    global points
    points, = ax.plot([], [], 'o')
    ani = animation.FuncAnimation(fig, anifunc, n_steps)
    plt.show()
    '''

if __name__ == '__main__':
    main()
