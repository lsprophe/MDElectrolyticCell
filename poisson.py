import numpy as np
from scipy.interpolate import CubicSpline

# GLOBALS
E_R = 1
E_0 = 1

def integrate_simple(pts, dx):
    # literally simplest possible, most stupid integration
    int_arr = np.zeros_like(pts)
    for i in range(1, len(pts)):
        int_arr[i] = pts[i-1]*dx + 0.5*pts[i]*dx
    
    return int_arr

def poisson_force(E_x_func, pos, charge):
    # NOTE: because we are not looking at potential distribution in
    # y that there is no y component to this force
    E_xp = E_x_func(pos[0])
    F_xp = E_xp * charge
    return F_xp

def poisson(nx, lx, particles, p0, pl):
    ''' Arguments:
    nx -> number of slices along x axis
    lx -> length of x axis
    particles -> list of particle objects
    p0 -> potential at x=0
    p1 -> potential at x=l
    '''
    dx = lx/nx
    x_arr = dx * range(nx)

    # array to hold the charge densities of each slice
    rho_e = np.zeros(nx)

    # set 0 charge density at anode and cathode
    rho_e[0] = 0
    rho_e[nx] = 0

    for n in range(1, nx-1):
        # find particles that are in the slice
        xmin, xmax = (n-1)*dx, n*dx
        for p in particles:
            if xmin <= p.q[0] < xmax:
                # add charge of particle to rho e for slice
                rho_e[n] += p.charge

    d2vdx2 = -rho_e / (E_R*E_0)

    # i guess integrate twice now?? applying the boundary conditions
    # here is kinda confusing to me...
    dvdx = integrate_simple(d2vdx2)

    # boundary conditions: dvdx is 0 at cathode and anode?
    dvdx[0] = 0
    dvdx[nx] = 0

    # This is the electric field value, so it needs to be returned too
    # to calculate the force on a given particle
    E_x = CubicSpline(x_arr, -1*dvdx)

    # integrate again to get potential?
    v = integrate_simple(dvdx)

    # set boundary potentials
    v[0] = p0
    v[nx] = pl

    # interpolate potential
    v_x = CubicSpline(x_arr, v)

    return E_x, v_x
