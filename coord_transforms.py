from _imports import *
from _functions import *

# Example problem from pg 88
# Given position vector of a satellite, rvec = 2Svec - Evec + 0.5 Zvec (DU). Radar site located at 169º W Long, 30º N Lat. Angle to Greenwich is 304º. Find rvec in IJK coordinates. Assume site is at sea level on a spherical Earth.
rvec_sez = np.array([[2, -1, 1.5]]).T # DU

# Longitude will be theta = theta_g + lambda_E
latitude = 30 # degrees
longitude = 304 + (-169) # degrees

ctm = coordTransformMatrix(latitude, longitude) # coordinate transform matrix
rvec_ijk = np.matmul(ctm, rvec_sez)
print(rvec_ijk)