import numpy as np

# GLOBALS
E_R = 1
E_0 = 1

def integrate_simple(pts, dx):
    # literally simplest possible, most stupid integration
    int_arr = np.zeros_like(pts)
    for i in range(1, len(pts)-1):
        int_arr[i] = pts[i-1]*dx + 0.5*pts[i]*dx
    
    return int_arr

def poisson(nx, lx, particles, p0, pl):
    ''' Arguments:
    nx -> number of slices along x axis
    lx -> length of x axis
    particles -> list of particle objects
    p0 -> potential at x=0
    p1 -> potential at x=l
    '''
    dx = lx/nx

    # array to hold the charge densities of each slice
    rho_e = np.empty(nx)

    # set 0 charge density at anode and cathode
    rho_e[0] = 0
    rho_e[nx] = 0
    for n in range(1, nx):
        # find particles that are in the slice
        xmin, xmax = (n-1)*dx, n*dx
        for p in particles:
            if xmin <= p.q[0] < xmax:
                # add charge of particle to rho e for slice
                rho_e[n] += p.charge

    d2vdx2 = -rho_e / (E_R*E_0)

    # i guess integrate twice now?? applying the boundary conditions
    # here is kinda confusing to me...
