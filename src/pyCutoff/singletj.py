import numpy as np

from tqdm import tqdm

from pymagglobal.utils import REARTH

from .utils import (
    geodetic_to_geocentric,
    rotate_direction_geodetic_to_geocentric,
    get_magnetic_field,
)

# Oxgen-16
ANUC = 16
ZCHARGE = 8
# Mass equivalent energy for one atomic mass unit
EMCSQ = 0.931141            # GeV
C = 2.99792458E5 / REARTH   # Speed of light per Earth's radius
DISOUT = 25.0               # Earth's radii; outer boundary
RHT = 20.0                  # km; inner boundary

# Geoid definition; XXX unify with proj transformer
ERPLSQ = 40408585.0
EREQSQ = 40680925.0
ERADPL = np.sqrt(ERPLSQ)
ERECSQ = EREQSQ / ERPLSQ - 1.0


def rhs(y, EOMC):
    # B_r, B_t, B_p
    b = get_magnetic_field(y[:3])

    res = np.array([
        C * y[3],
        C * y[4] / y[0],
        C * y[5] / y[0] / np.sin(y[1]),
        EOMC * (y[4] * b[2] - y[5] * b[1])
        + C * (y[4]**2 + y[5]**2) / y[0],
        EOMC * (y[5] * b[0] - y[3] * b[2])
        - C * (
            (y[3] * y[4]) / y[0]
            - y[5]**2 / y[0] / np.tan(y[1])
        ),
        EOMC * (y[3] * b[1] - y[4] * b[0])
        - C * (
            (y[3] * y[5]) / y[0]
            + y[4] * y[5] / y[0] / np.tan(y[1])
        ),
    ])

    return res


def singletj(gdlat, gdlon, rigidity, direction_gd=np.array([1, 0, 0])):
    TENG = np.sqrt(
        (rigidity * ZCHARGE)**2 + (ANUC * EMCSQ)**2
    )
    EOMC = -8987.566297 * ZCHARGE / TENG
    BETAST = 2

    zed = 0
    gclat, gclon, rad = geodetic_to_geocentric(gdlat, gdlon, 20e3)
    gccolat = 90 - gclat
    # This is the initial point of the trajectory
    initial_position = np.array(
        [rad / 6371.2, np.deg2rad(gccolat), np.deg2rad(gclon)]
    )

    # Transform directions
    direction_gc = rotate_direction_geodetic_to_geocentric(
        direction_gd,
        gdlat,
        gdlon,
        alt=20e3,
    )

    # Lorentz factor
    # Calculate relativistic gamma factor from rigidity = velocity / charge
    # XXX Careful with the units!
    gamma_factor = np.sqrt(
        1.0 + ((rigidity * ZCHARGE) / (EMCSQ * ANUC))**2
    )
    beta_factor = np.sqrt(1. - 1. / gamma_factor**2)

    edif = beta_factor * 1.0e-4
    if (edif < 1.0-5):
        edif = 1.0e-5
    if (beta_factor < 0.1):
        edif = 1.0e-4

    inital_velocity = beta_factor * direction_gc

    y = np.concatenate(
        (
            initial_position, inital_velocity
        )
    )

    tcy2 = np.cos(initial_position[1])
    ptcy2 = abs(tcy2)

    ahlt = (1.0 + ptcy2)**2
    h_start = 6.0e-6 * rigidity / (beta_factor * ahlt + zed * ptcy2)
    if (h_start <= 1.0e-6):
        h_start = 1.0e-6

    # hold = h_start
    h_ck = h_start
    h_cng = h_start
    h_max = 1 / C / beta_factor
    # disck = DISOUT - 1.1 * h_max * C * beta_factor
    kbf = 0

    # One RK4 step
    def RK4_step(y, h):
        k1 = rhs(y, EOMC)
        k2 = rhs(y + h / 2 * k1, EOMC)
        k3 = rhs(y + h / 2 * k2, EOMC)
        k4 = rhs(y + h * k3, EOMC)

        y_dot = (k1 + 2 * k2 + 2 * k3 + k4) / 6
        return y + h * y_dot

    # Checks go here
    def check_RK4_step():
        nonlocal kbf, BETAST

        # CHECK FOR  UNACCEPTABLE CHANGES IN BETA
        rck_beta = np.sqrt(np.sum(y[3:]**2))
        tbeta = beta_factor - rck_beta
        if abs(tbeta) > edif:
            kbf += 1
            BETAST = BETAST + ahlt
            BETAST = BETAST + kbf * (1.0 + ahlt)
        if kbf > 5:
            print(
                'Irrecoverable beta error for trajectory '
                f'{gdlat:.2f}, {gdlon:.2f}, {rigidity:.2f}.'
            )
            return False

        return True

    tau = 0.
    n_steps = 0
    irslt = np.nan
    acc_old = np.nan * np.ones(3)
    acc_abs_old = np.nan
    for it in tqdm(range(5000), disable=True):
        # XXX this should be retrieved from the rhs!
        B = get_magnetic_field(y[:3])
        B_mag = np.sqrt(np.sum(B**2))
        hb = 1.6e-5 * rigidity / (B_mag * beta_factor)
        h = hb / BETAST

        if (h_ck < 1.0e-6):
            h_ck = 1.0e-6
        if (kbf > 0):
            h = h / (2 * kbf)
        if (h > h_max):
            h = h_max
        if (h > h_cng):
            h = h_cng
        if (h > h_ck):
            h = h_ck

        # perform RK4 step
        y_new = RK4_step(y, h)
        acc = rhs(y_new, EOMC)[3:]
        B = get_magnetic_field(y)

        acc_abs = np.sqrt(np.sum(acc**2))
        h_cng = h * 1.2
        h_ck = h_cng

        if not check_RK4_step():
            irslt = 0
            break

        if 0 < it:
            if (acc_abs - acc_abs_old > 5):
                h_ck = h_ck / (1.0 + ahlt)
                if (acc_abs / acc_abs_old > 2):
                    h_ck = h_ck / (1.0 + ahlt)

            if np.any(np.abs(acc_old) > 3):
                # XXX dirty hack; is a for loop faster?
                p = np.sum(np.abs(acc / acc_old) > 3)
                h_ck = h_ck / (1.0 + ahlt)**p

        if DISOUT < y_new[0]:
            if (
                h < 1.0e-5
                or h_ck < 1.0e-5
                or h_cng < 1.0e-5
            ):
                irslt = 1
                break
            else:
                h_ck /= 2.0
                h_cng /= 2.0
                continue
        if 4 < it:
            # XXX Unify with proj transformer
            grnd_km = (
                ERADPL / np.sqrt(
                    1.0 - ERECSQ * np.sin(y_new[1])**2
                )
            )

            if (y_new[0] < (RHT + grnd_km) / REARTH):
                irslt = -1
                break

        # increment
        tau += h
        n_steps += 1
        y = y_new
        acc_abs_old = acc_abs
        acc_old = acc

    # If irslt still 0, no conclusion has been reached within the maximum
    # number of steps
    if np.isnan(irslt):
        irslt = 0
    # GDLATD,GLOND,PC,ZED,AZD,ISALT,NSTEP,IFATE,CNAME

    TCY2 = np.cos(y_new[1])
    TSY2 = np.sin(y_new[1])
    YDA5 = y_new[4] * TCY2 + y_new[3] * TSY2
    ATRG1 = y_new[3] * TCY2 - y_new[4] * TSY2
    ATRG2 = np.sqrt(y[5]**2 + YDA5**2)

    as_lat = np.rad2deg(np.arctan2(ATRG1, ATRG2))
    as_lon = np.rad2deg(y_new[2])
    if (y_new[5] != 0 and YDA5 != 0.0):
        as_lon += np.rad2deg(np.arctan2(y[5], YDA5))

    if (as_lon < 0.0):
        as_lon += 360
    if (as_lon > 360.0):
        as_lon -= 360

    return gdlat, gclat, gdlon, rigidity, as_lat, as_lon, n_steps, irslt
