import numpy as np

# globals
K_E = 8.988e9  #coulomb constant

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
    field_mag = 24*eps*((sigma/r_ij)**5)

    # Update electric field vector
    f_x = field_mag * dx/r_ij
    f_y = field_mag * dy/r_ij

    return (f_x, f_y)

def coulomb_attractive(l1, l2, charge1, charge2):
    # Calculate interparticle distances
    dx = l1[0] - l2[0]
    dy = l1[1] - l2[1]
    r_ij = np.sqrt(dx**2 + dy**2)

    f_mag = K_E*(abs(charge1)*abs(charge2)) / (r_ij**2)

    f_x = f_mag * (dx/r_ij)
    f_y = f_mag * (dy/r_ij)

    return (f_x, f_y)

def LJF_repulsive(l1, l2, eps, sigma):
    # Calculate interparticle distances
    dx = l1[0] - l2[0]
    dy = l1[1] - l2[1]
    r_ij = np.sqrt(dx**2 + dy**2)

    # Calculate magnitude of repulsive field according to Lennard-Jones Potential
    field_mag = 24*eps*(2*(sigma/r_ij)**11)

    # Update electric field vector
    f_x = field_mag * dx/r_ij
    f_y = field_mag * dy/r_ij

    return (f_x, f_y)

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

def pair_pot(particles): 
    ''' Arguments:
    particles -> list of particle objects
    '''
    ### WIP - pair potential should likely use coulomb interaction instead of Lennard-Jones!!

    # Define array to hold electric field for each particle
    # pp : "Pair Potential"
    pp_force = np.zeros([len(particles),2])

    # Loop for each particle
    for idx_i, p_i in range(len(pp_force)):
        x_i, y_i = p_i.q # Set initial position
        pp_field = [0,0] # Initialize field

        # Loop over each interacting particle
        idx_j = 0
        for p_j in particles:
            # Check for self-interaction
            if idx_j != idx_i:
                # Calculate interparticle distances
                dx = p_j.q[0] - x_i
                dy = p_j.q[1] - y_i
                r_ij = np.sqrt(dx**2 + dy**2)

                # Calculate electric field according to Lennard-Jones Potential
                field_mag = -p_j.charge*24*EPS*(2*(SIGMA/r_ij)**11-(SIGMA/r_ij)**5)

                # Update electric field vector
                pp_field[0] += field_mag * dx/r_ij
                pp_field[1] += field_mag * dy/r_ij

            idx_j += 1

        # Update force based on electric field
        pp_force[idx_i] = pp_field

    return pp_force

def combined_forcefield(lx, particles, p0, pl)