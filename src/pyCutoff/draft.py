from .singletj import singletj

# from scipy.integrate import RK45

# Input
# -----------------------------------------------------------------------------
rigidity = 25
gdlat = 89.00
gdlon = 0.00
# TODO: Calculate direction from azimuth and zenith
# direction_gd = np.array([1, 0, 0])

rigidities = [25, 35, 5, 50]
gdlats = [89.00, 40, 0, -63]
gdlons = [0.00, -75, 169, 23]

for gdlat, gdlon, rigidity in zip(gdlats, gdlons, rigidities):
    singletj(gdlat, gdlon, rigidity)
