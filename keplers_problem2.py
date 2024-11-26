# Solve Kepler's Problem
# Given r0, v0 and ∆nu, solve for r and v.

# Imports
import numpy as np
from _functions import magnitude, sign
import math as m

# Constants
mu = 1 # DU^3/TU^2
mur = np.sqrt(mu) # sqrt(mu), since it pops up so much

# Given values (r0 and v0 in IJK system)
r0 = np.array([1, 1, 0]) # initial position vector, DU
v0 = np.array([0, 0, 2]) # intitial velocity vector, DU/TU
dnu = 60 * np.pi / 180 # ∆nu in radians

# Derivatives
r0m = magnitude(r0) # magnitude of r0, DU
v0m = magnitude(v0) # magnitude of v0, DU/TU

# Step 1: From r0, v0, determine r0m, a, e, p, and nu0-------------------------

# Orbital constants
hvec = np.cross(r0, v0) # angular momentum, DU^2/TU
h = magnitude(hvec) # magnitude of hvec, DU^2/TU
p = h**2 / mu # orbital parameter, DU
sme = v0m**2 / 2 - mu / r0m # specific mechanical energy, DU^2/TU^2
a = - mu / (2 * sme) # semi-major axis, DU
if a > 0:
    # e = np.sqrt(1 - (2 * sme * h**2 / mu**2)) # eccentricity, unitless
    e = np.sqrt(1 - p / a)
    orbit = 'e' if 0 < e < 1 else 'p'
else:
    e = np.sqrt(1 - p / a)
    orbit = 'h'

# Calculate true anomaly at t0
nu0 = np.arccos(np.round(1 / e * (p / r0m - 1), 10)) # true anomaly at t0, radians
nu  = nu0 + dnu # true anomaly at t, radians
if nu > 2 * np.pi:
    nu -= 2* np.pi
    print('nu reduction')

# Calculate eccentric anomalies and TOF for specific orbit type
if orbit == 'e':    # elliptical orbit
    ea0 = np.arccos((e + np.cos(nu0)) / (1 + e * np.cos(nu0)))
    ea  = np.arccos((e + np.cos(nu))  / (1 + e * np.cos(nu)))
    k = 1 if nu0 > nu else 0
    deltat = np.sqrt(a**3 / mu) * \
            (2 * np.pi * k + (ea - e * np.sin(ea)) - (ea0 - e * np.sin(ea0)))
    
elif orbit == 'h':   # hyperbolic orbit
    fa0 = np.arccosh((e + np.cos(nu0)) / (1 + e * np.cos(nu0)))
    fa0 *= -1 if np.pi < nu0 < 2 * np.pi else fa0
    fa  = np.arccosh((e + np.cos(nu))  / (1 + e * np.cos(nu)))
    fa *= -1 if np.pi < nu0 < 2 * np.pi else fa
    deltat = np.sqrt((-a)**3 / mu) * \
            ((e * np.sinh(fa) - fa) - (e * np.sinh(fa0) - fa0))
    
elif orbit == 'p':   # parabolic orbit
    da0 = np.sqrt(p) * np.tan(nu0 / 2)
    da  = np.sqrt(p) * np.tan(nu / 2)
    deltat = (2 * mur)**-1 * ((p * da + 1/3 * da**3) - (p * da0 + 1/3 * da0**3))


# Step 2: Given (t - t0), solve universal time-of-flight (TOF) equation for x using Newton iteration scheme -------------------------------------------------

# Define C(z), S(z)
if orbit == 'e':
    C = lambda z: (1 - np.cos(np.sqrt(z))) / 2
    S = lambda z: (np.sqrt(z) - np.sin(np.sqrt(z))) / np.sqrt(z**3)
elif orbit == 'h':
    C = lambda z: (1 - np.cosh(np.sqrt(-z))) / z
    S = lambda z: (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / np.sqrt((-z)**3)
elif orbit == 'p':
    C = lambda z: np.sum((-z)**k / m.factorial(2 * k + 2) for k in range(10))
    S = lambda z: np.sum((-z)**k / m.factorial(2 * k + 3) for k in range(10))

# x = np.sqrt(a) * (ea - ea0) if orbit == 'e' else np.sqrt(-a) * (fa - fa0)
# z = lambda x: x**2 / a

# t = x**3 * S(z(x)) + np.dot(r0, v0) / mur * x**2 * C(z(x)) + r0m * x * (1 - z(x) * S(z(x)))



# # Newton Iteration Scheme -----------------------------------------------------

# First guess for x:
# if orbit == 'e':
#     xn = mur * deltat / a
# elif orbit == 'h':
#     xn = np.sin(deltat) * np.sqrt(-a) * np.log(
#         (-2 * mu * deltat) /
#         (a * (np.dot(r0, v0) + np.sign(deltat) * np.sqrt(-a * mu) * (1 - r0m/a)))
#     )

# Define variable z
# z = lambda x: x**2 / a

# # Define differential equaiton dt/dx
# dtdx = lambda x: x**2 * C(z(x)) / mur + h / mur * (1 - z(x) * S(z(x)) + r0m / mur * (1 - z(x) * C(z(x))))
# dmde = lambda ea: 1 - e * np.cos(ea) # need to go back and convert f to e

# # Define t_n equation
# tn = lambda xn: h / mu * xn**2 * C(z(xn)) + (1 - r0m/a) * xn**3 * S(z(xn)) / mur + r0m * xn / mur
# mn = lambda en: en - e * np.sin(en)

# # Solve for mean anomaly, M
# ma = np.sqrt(mu / )

# # Iterate until convergence (t ~ tn)
# for n in range(10):
#     # tn1 = tn(xn)
#     # xn1 = xn + (t - tn1) / dtdx(xn)
#     # tn2 = tn(xn1)
#     mn1 = mn(xn)
#     xn1 = xn + (mn - mn1) / dmde(xn1)
#     mn2 = mn(xn1)

#     # if np.abs(t - tn2) < 1e-4:
#     if np.abs(mn - mn2) < 1e-4:
#         break
#     else:
#         # tn1 = tn2
#         xn = xn1
#         mn1 = mn2

# # Universal variable, x
# x = xn1


# Step 3: Evaluate f & g then compute r ---------------------------------------
# f = 1 - (a / r0m) * (1 - np.cos(x / np.sqrt(a)))
# g = t - x**3 / mur * S(z(x))
# If ∆nu known:
dfa = fa - fa0
print(dfa)
rm = a * (1 - e * np.cosh(fa))
f = 1 - (a / r0m) * (1 - np.cosh(dfa))
fdot = - np.sqrt(- mu * a) * np.sinh(dfa) / rm / r0m

g = deltat - np.sqrt((-a)**3 / mu) * (np.sinh(dfa) - dfa)
gdot = 1 - (a / rm) * (1 - np.cosh(dfa))
# f = 1 - (rm / p) * (1 - np.cos(dnu))
# g = (rm * r0m * np.sin(dnu)) / np.sqrt(mu * p)
# galt = deltat - np.sqrt((-a)**3/mu) * ((-fa + fa0) - np.sin(-fa + fa0))
# print(g, galt)

# Compute r using f & g
r = f * r0 + g * v0
# rm = magnitude(r)

# Step 4: Evaluate fdot & gdot then compute v ---------------------------------
# fdot = - np.sqrt(mu * a) / r0m / rm * np.sin(x / np.sqrt(a))
# gdot = 1 - a / rm + a / rm * np.cos(x / np.sqrt(a))
# fdot = np.sqrt(mu / p) * np.tan(dnu / 2) * ((1 - np.cos(dnu)) / p - rm**-1 - r0m**-1)
# gdot = 1 - (r0m / p) * (1 - np.cos(dnu))

# Compute v using fdot & gdot
v = fdot * r0 + gdot * v0
vm = magnitude(v)

t = deltat - np.sqrt((-a)**3 / mu) * (e * np.sinh(fa0) - fa0)
# np.sqrt((-a)**3 / mu) * \
#             ((e * np.sinh(fa) - fa) - (e * np.sinh(fa0) - fa0))

# Output results
w = 3 # number of places to round to
print('\n'
      'Final Position Vector r(t={:.2f}) = {}{} I {} {} J {} {} K\n'
      'Final Velocity Vector v(t={:.2f}) = {}{} I {} {} J {} {} K\n'
      ''.format(
          t, sign(r[0]), np.round(np.abs(r[0]), w), sign(r[1]), np.round(np.abs(r[1]), w), sign(r[2]), np.round(np.abs(r[2]), w),
          t, sign(v[0]), np.round(np.abs(v[0]), w), sign(v[1]), np.round(np.abs(v[1]), w), sign(v[2]), np.round(np.abs(v[2]), w),
      ))