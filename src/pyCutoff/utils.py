import numpy as np

import pyproj

from pathlib import Path

from scipy.special import factorial, factorial2

from pymagglobal.utils import REARTH, lm2i, lmax2N
from paleokalmag.utils import dsh_basis


path = str(Path(__file__).parent) + '/'

# XXX Replace by WGS84 once checks are complete
transformer = pyproj.Transformer.from_crs(
    pyproj.CRS.from_proj4(
        "+proj=latlon +a=6378160.001128852 +b=6356774.732519629"
    ),
    pyproj.CRS.from_proj4(
        "+proj=geocent +a=6378160.001128852 +b=6356774.732519629"
    ),
)
l_max = 10

coeffs = np.zeros(lmax2N(l_max))
with open(path + 'coeffs.txt') as fh:
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


def geodetic_to_geocentric(gdlat, gdlon, alt=0):
    # Transform locations
    ecef_x, ecef_y, ecef_z = transformer.transform(
        gdlon,
        gdlat,
        alt,
        radians=False,
    )

    rad = np.sqrt(
        ecef_x**2 + ecef_y**2 + ecef_z**2
    )
    gccolat = np.rad2deg(np.arccos(ecef_z / rad))
    gclon = np.rad2deg(np.arctan2(ecef_y, ecef_x))
    gclat = 90 - gccolat
    rad /= 1e3

    return gclat, gclon, rad


def rotate_direction_geodetic_to_geocentric(direction, gdlat, gdlon, alt=0):
    gclat, _, _ = geodetic_to_geocentric(gdlat, gdlon, alt=alt)
    angle = gdlat - gclat
    _sin = np.sin(np.deg2rad(angle))
    _cos = np.cos(np.deg2rad(angle))

    return np.array(
        [
            direction[0]*_cos + direction[1] * _sin,
            -direction[0]*_sin + direction[1] * _cos,
            direction[2]
        ]
    )


def get_magnetic_field(position):
    z_at = np.atleast_2d(
        [
            np.rad2deg(position[1]),
            np.rad2deg(position[2]),
            position[0] * REARTH,
        ],
    ).T
    # N,E,C in nT
    b = coeffs @ dsh_basis(l_max, z_at)
    # B_r = -B_C = -B[2]
    # B_t = -B_N = -B[0]
    # B_p = B_E = B[1]

    # R,T,P in Gauss (1 Gauss = 1e-4 T = 1e5 nT)
    return 1e-5 * np.array([-b[2], -b[0], b[1]])
