import unittest
from pathlib import Path

import numpy as np

import pyproj
from scipy.special import factorial, factorial2

from pymagglobal.utils import lm2i, lmax2N

from pyCutoff import constants
# Monkey patch older Ellipsoid
constants.transformer = pyproj.Transformer.from_crs(
    pyproj.CRS.from_proj4(
        "+proj=latlon +a=6378160.001128852 +b=6356774.732519629"
    ),
    pyproj.CRS.from_proj4(
        "+proj=geocent +a=6378160.001128852 +b=6356774.732519629"
    ),
)

from pyCutoff.magnetic_field import MagneticField
from pyCutoff.singletj import singletj

path = str(Path(__file__).parent) + '/'


l_max = 10

coeffs = np.zeros(lmax2N(l_max))
with open(path + 'dat/coeffs.txt') as fh:
    for line in fh.readlines():
        ell, emm, val = np.fromstring(line, dtype=float, sep=' ')
        ell = int(ell)
        emm = int(emm)
        # Conversion factor for Gauss normalization to Schmidt-semi
        # see for example https://ntrs.nasa.gov/api/citations/19900004113/
        # downloads/19900004113.pdf Page 272
        fac = (
            np.sqrt(
                (2 - (emm == 0)) * factorial(ell - emm)
                / (factorial(ell + emm))
            ) * factorial2(2 * ell - 1) / factorial(ell - emm)
        )
        coeffs[lm2i(ell, emm)] = val / fac


magField = MagneticField(coeffs)


class TestLines(unittest.TestCase):
    def test_lines(self):
        with open(path + 'dat/TAPE8.TXT', 'r') as fh:
            it = 0
            for line in fh.readlines():
                values = line.split()
                gdlat = float(values[0])
                gdlon = float(values[2])
                rigidity = float(values[5])

                res = singletj(
                    gdlat,
                    gdlon,
                    rigidity,
                    magField,
                )

                self.assertTrue(round(res[1], 2) == float(values[1]))
                self.assertTrue(1 - res[7] == int(values[9]))
                if res[7] == 1:
                    self.assertTrue(round(res[4], 2) == float(values[6]))
                    self.assertTrue(round(res[5], 2) == float(values[7]))

                it += 1
                if 4 <= it:
                    break


if __name__ == '__main__':
    unittest.main()
