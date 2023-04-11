from enum import Enum

class ParticleType(Enum):
    PROTON = 1
    CATHODE_ION = 2
    ANODE_ION = 3
    SULFATE = 4

AVO = 6.022*10**(23)
eps = 1000 * (1/AVO) * (10**9)
# lennard jones information for solution ions
LJ_DICT = {
    ParticleType.SULFATE: {
    ParticleType.PROTON: {
    "eps": eps,
    "sigma": 0.1215
    },
    ParticleType.SULFATE: {
    "eps": eps,
    "sigma": 0.242
    },
    ParticleType.ANODE_ION: {
    "eps": eps,
    "sigma": 0.71
    },
    ParticleType.CATHODE_ION: {
    "eps": eps,
    "sigma": 0.4983
    }
    }
}
