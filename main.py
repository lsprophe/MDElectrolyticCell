import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

from md import verlet_integrator, Particle
from poisson import PoissonPES

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
    k = 1
    q = np.array([0, 0])
    v = np.array([0, 0])
    mass = 1
    NP = 500  # number of positive particles
    NN = 500 # number of negative particles
    n_steps = 10**3

    poisson_pes = PoissonPES(nx=100, lx=1, p0=0, pl=1)

    # randomize initial particle positions and velocities
    particles = [Particle(mass, np.random.rand(2), 2*np.random.uniform(1,2)-1, -e, force_fun=poisson_pes.force) for n in range(NN)]
    particles += [Particle(mass, np.random.rand(2), 2*np.random.uniform(1,2)-1, e, force_fun=poisson_pes.force) for n in range(NP)]

    verlet_integrator(particles, poisson_pes, n_steps, 0.1, system='NVT', gamma=10)

    fig, (ax1, ax2) = plt.subplots(1,2)
    
    # For testing, plots x and y positions on their own axes for each particle
    for n in range(NN + NP):
        X = [point[0] for point in particles[n].q_arr]
        Y = [point[0] for point in particles[n].q_arr]
        print(X)
        ax1.plot(X)
        ax2.plot(Y)
    
    ax1.title.set_text(" Y positions")
    ax2.title.set_text("X positions")
    plt.show()

    # NOTE: This is a WIP to do animation of particle positions
    '''
    global points
    points, = ax.plot([], [], 'o')
    ani = animation.FuncAnimation(fig, anifunc, n_steps)
    plt.show()
    '''

if __name__ == '__main__':
    main()
