import numpy as np

# globals
EPS = 40*4184/
SIGMA = 2.4e-10 # m

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