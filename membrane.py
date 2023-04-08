import numpy as np
from dataclasses import dataclass

@dataclass
class Material:
    pore_cat_sigma = 1
    pore_cat_eps = 1

    pore_an_sigma = 1
    pore_an_eps = 1

    base_cat_sigma = 1
    base_cat_eps = 1

    base_an_sigma = 1
    base_an_eps = 1

    hydration = 5 # mols h2o/mol pore molecule

class Membrane:
    def __init__(self, porosity, psd: tuple, thickness, x_centre, height, material, n_slices):
        self.porosity = porosity
        self.x_l = x_centre -(thickness/2)
        self.x_r = x_centre + (thickness/2)


def pore_locations(psd, porosity, x_centre, height, thickness, n_slices, cf):
    # setup membrane pores using the porosity and pore
    # size distribution
    n_pores = (0.4*height) / psd(0)
    

