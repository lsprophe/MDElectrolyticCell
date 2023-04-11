import numpy as np
from MD.types import LJ_DICT, ParticleType

# globals
K_E = 8.988e27  #coulomb constant (N*nm^2/C^2)

# simple two atom potential functions
def LJF_attractive(l1, l2, eps, sigma):
    ''' Arguments:
    l1 -> location of particle 1
    l2 -> location of particle 2
    eps -> epsilon for lennard jones
    sigma -> sigma for lennard jones
    '''
    # Calculate interparticle distances
    dx = l1[0] - l2[0]
    dy = l1[1] - l2[1]
    r_ij = np.sqrt(dx**2 + dy**2)

    # Calculate magnitude of attractive field according to Lennard-Jones Potential
    field_mag = -24e-9*eps*(sigma**6)/(r_ij**7) # Includes negative to account for attraction

    # Update electric field vector
    f_x = field_mag * dx/r_ij
    f_y = field_mag * dy/r_ij
    return np.array([f_x, f_y])

def coulomb(l1, l2, charge1, charge2):
    # Calculate interparticle distances
    dx = l1[0] - l2[0]
    dy = l1[1] - l2[1]
    r_ij = np.sqrt(dx**2 + dy**2)

    f_mag = K_E*(charge1*charge2 / (r_ij**2)) # Includes sign to indicate attractive (-) or repulsive (+)

    f_x = f_mag * (dx/r_ij)
    f_y = f_mag * (dy/r_ij)
    return np.array([f_x, f_y])

def LJF_repulsive(l1, l2, eps, sigma):
    # Calculate interparticle distances
    dx = l1[0] - l2[0]
    dy = l1[1] - l2[1]
    r_ij = np.sqrt(dx**2 + dy**2)

    # Calculate magnitude of repulsive field according to Lennard-Jones Potential
    field_mag = 48e-9*eps*(sigma**12)/(r_ij**13)

    # Update electric field vector
    f_x = field_mag * dx/r_ij
    f_y = field_mag * dy/r_ij
    return np.array([f_x, f_y])

def ext_field(lx, particles, p0, pl):
    ''' Arguments:
    lx -> length of x axis
    particles -> list of particle objects
    p0 -> potential at x=0
    p1 -> potential at x=l
    '''
    # Define array to hold force due to external field for each particle
    # ef : "External Field"
    ef_force = np.zeros([len(particles),2])

    # Loop for each particle
    for idx, p in enumerate(particles):
        # Calculate applied electric field
        # This assumes that applied electric field is constant!
        ef_field = -(pl-p0)/lx
        # Calculate force due to electric field
        # Note that force is only in the x-direction
        ef_force[idx][0] = ef_field*p.charge
    
    return ef_force

def pair_pot(particle, particles): 
    ''' Arguments:
    particle -> particular particle that is being investigated
    particles -> list of particle objects
    '''
    force = np.zeros(2)

    for p_j in particles:
        # Check for self-interaction
        if p_j is not particle:
            # Calculate interparticle distances
            l_j = p_j.q # Interacting particle position
            charge_j = p_j.charge # Interacting particle charge

            # Calculate electric field according to Coulomb interactions
            (f_x,f_y) = coulomb(particle.q, l_j, particle.charge, charge_j)
            # Update pair-potential force
            force += np.array([f_x, f_y])

            # if one of the particles is a sulfate group, also use the repulsive
            # lennard-jones force
            t1 = particle.type
            t2 = p_j.type
            if t1 is ParticleType.SULFATE:
                lf_p = LJ_DICT[t1][t2]
                ljf = LJF_repulsive(particle.q, l_j, lf_p["eps"], lf_p["sigma"])
                force += ljf
            elif t2 is ParticleType.SULFATE:
                lf_p = LJ_DICT[t1][t2]
                ljf = LJF_repulsive(particle.q, l_j, lf_p["eps"], lf_p["sigma"])
                force += ljf

    return force
