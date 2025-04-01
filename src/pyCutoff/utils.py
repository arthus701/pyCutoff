from pathlib import Path

import numpy as np

from .constants import transformer

path = str(Path(__file__).parent) + '/'


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
            direction[0] * _cos + direction[1] * _sin,
            -direction[0] * _sin + direction[1] * _cos,
            direction[2] * np.ones_like(gdlat)
        ]
    )
