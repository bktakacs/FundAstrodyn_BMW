# Solve Kepler's Problem
# Given r0 and v0 at time t0 = 0, solve for r and v at time t.

# Imports
import numpy as np
from _functions import magnitude, sign

# Constants
mu = 1 # DU^3/TU^2
mur = np.sqrt(mu) # sqrt(mu), since it pops up so much

# Given values (r0 and v0 in IJK system)
r0 = np.array([1, 0, 0]) # initial position vector, DU
v0 = np.array([0, 0, 1.1]) # intitial velocity vector, DU/TU
t0 = 0 # initial epoch, TU
t = 2 # epoch at which we want to find r and v

# Derivatives
r0m = magnitude(r0) # magnitude of r0, DU
v0m = magnitude(v0) # magnitude of v0, DU/TU

# Step 1: From r0, v0, determine r0m and a ------------------------------------

# Orbital constants
h = np.dot(r0, v0.T) # angular momentum, DU^2/TU
p = h**2 / mu # orbital parameter, DU
sme =v0m**2 / 2 - mu / r0m # specific mechanical energy, DU^2/TU^2
a = - mu / (2 * sme) # semi-major axis, DU

# Step 2: Given (t - t0), solve universal time-of-flight (TOF) equation for x using Newton iteration scheme -------------------------------------------------

# Define C(z), S(z)
C = lambda z: (1 - np.cos(np.sqrt(z))) / 2
S = lambda z: (np.sqrt(z) - np.sin(np.sqrt(z))) / np.sqrt(z**3)

# Newton Iteration Scheme -----------------------------------------------------

# First guess for x:
xn = mur * (t - t0) / a

# Define variable z
z = lambda x: x**2 / a

# Define differential equaiton dt/dx
dtdx = lambda x: x**2 * C(z(x)) / mur + h / mur * (1 - z(x) * S(z(x)) + r0m / mur * (1 - z(x) * C(z(x))))

# Define t_n equation
tn = lambda xn: h / mu * xn**2 * C(z(xn)) + (1 - r0m/a) * xn**3 * S(z(xn)) / mur + r0m * xn / mur

# Iterate until convergence (t ~ tn)
for n in range(10):
    tn1 = tn(xn)
    xn1 = xn + (t - tn1) / dtdx(xn)

    tn2 = tn(xn1)

    if np.abs(t - tn2) < 1e-4:
        break
    else:
        tn1 = tn2
        xn = xn1

# Universal variable, x
x = xn1

# Step 3: Evaluate f & g then compute r ---------------------------------------
f = 1 - (a / r0m) * (1 - np.cos(x / np.sqrt(a)))
g = t - x**3 / mur * S(z(x))

# Compute r using f & g
r = f * r0 + g * v0
rm = magnitude(r)

# Step 4: Evaluate fdot & gdot then compute v ---------------------------------
fdot = - np.sqrt(mu * a) / r0m / rm * np.sin(x / np.sqrt(a))
gdot = 1 - a / rm + a / rm * np.cos(x / np.sqrt(a))

# Compute v using fdot & gdot
v = fdot * r0 + gdot * v0
vm = magnitude(v)

# Output results
w = 2 # number of places to round to
print('\n'
      'Final Position Vector r(t={}) = {}{} I {} {} J {} {} K\n'
      'Final Velocity Vector v(t={}) = {}{} I {} {} J {} {} K\n'
      ''.format(
          t, sign(r[0]), np.round(np.abs(r[0]), w), sign(r[1]), np.round(np.abs(r[1]), w), sign(r[2]), np.round(np.abs(r[2]), w),
          t, sign(v[0]), np.round(np.abs(v[0]), w), sign(v[1]), np.round(np.abs(v[1]), w), sign(v[2]), np.round(np.abs(v[2]), w),
      ))