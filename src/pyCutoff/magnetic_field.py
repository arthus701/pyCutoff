import numpy as np

from pymagglobal.utils import i2lm_l
from paleokalmag.utils import dsh_basis

from .constants import REARTH


class MagneticField(object):
    def __init__(self, coeffs):
        self.coeffs = coeffs
        self.l_max = i2lm_l(int(coeffs.shape[0]-1))

    def __call__(self, position):
        z_at = np.atleast_2d(
            [
                np.rad2deg(position[1]),
                np.rad2deg(position[2]),
                position[0] * REARTH,
            ],
        ).T
        # N,E,C in nT
        b = self.coeffs @ dsh_basis(self.l_max, z_at)
        # B_r = -B_C = -B[2]
        # B_t = -B_N = -B[0]
        # B_p = B_E = B[1]

        # R,T,P in Gauss (1 Gauss = 1e-4 T = 1e5 nT)
        return 1e-5 * np.array([-b[2], -b[0], b[1]])
