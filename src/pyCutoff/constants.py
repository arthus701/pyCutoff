import pyproj

REARTH = 6371.2             # Earth's radius in km

# Oxgen-16
ANUC = 16
ZCHARGE = 8
# Mass equivalent energy for one atomic mass unit
EMCSQ = 0.931141            # GeV
C = 2.99792458E5 / REARTH   # Speed of light per Earth's radius

transformer = pyproj.Transformer.from_crs(
    pyproj.CRS.from_proj4(
        "+proj=latlon +ellps=WGS84 +datum=WGS84"
    ),
    pyproj.CRS.from_proj4(
        "+proj=geocent +ellps=WGS84 +datum=WGS84"
    ),
)
