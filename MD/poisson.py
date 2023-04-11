import numpy as np
from scipy.interpolate import CubicSpline

# GLOBALS
E_R = 1
E_0 = 1

class PoissonPES:
    def __init__(self, nx, lx, p0, pl):
        self.nx = nx
        self.lx = lx
        self.p0 = p0
        self.pl = pl

        dx = lx/nx
        self.x_arr = dx * np.array(range(nx))
        self.v_func = None
        self.E_func = None

    def calculate(self, particles):
        self.E_func, self.v_func = poisson(self.nx, self.lx, particles, self.p0, self.pl)


def poisson(nx, lx, particles, p0, pl):
    ''' Arguments:
    nx -> number of slices along x axis
    lx -> length of x axis
    particles -> list of particle objects
    p0 -> potential at x=0
    p1 -> potential at x=l
    '''
    dx = lx/nx
    x_arr = dx * np.array(range(nx))

    # array to hold the charge densities of each slice
    rho_e = np.zeros(nx)

    # set 0 charge density at anode and cathode
    rho_e[0] = 0
    rho_e[nx-1] = 0

    for n in range(1, nx-1):
        # find particles that are in the slice
        xmin, xmax = (n-1)*dx, n*dx
        for p in particles:
            if xmin <= p.q[0] < xmax:
                # add charge of particle to rho e for slice
                rho_e[n] += p.charge

    dEdx = rho_e / (E_R*E_0)
    # boundary condition: E(0) = 0
    E = np.cumsum(dEdx)

    # Apply fixed boundary conditions: step 1: find mean field value
    #e_mean = np.sum(E) * (1/lx)
    # step 2: add the following value to simulate eigenvalue problem
    #E += (pl - p0)/(lx - e_mean)

    # This is the electric field value, so it needs to be returned too
    # to calculate the force on a given particle
    E_func = CubicSpline(x_arr, E)

    # integrate again to get potential?
    v = np.cumsum(-1*E)

    # interpolate potential
    v_func = CubicSpline(x_arr, v)

    return E_func, v_func
