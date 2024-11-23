# Chapter 2 Exercises

# Imports
from _imports import *
from _functions import *

# 2.13 Determine orbital elements of following objects using Gibbs method. Units are earth canonical units.

# a.
r1a = np.array([[1.41422511,          0,  1.414202]])
r2a = np.array([[1.81065659, 1.06066883,  0.3106515]])
r3a = np.array([[1.35353995, 1.41422511, -0.6464495]])
print('\n\tPart a:')
orbitalElemGibbs(r1a, r2a, r3a)

# b.
r1b = [[0.70711255,            0,  0.70710101]]
r2b = [[-0.89497879,  0.56568081, -0.09496418]]
r3b = [[-0.09497879, -0.56568081, -0.89497724]]
print('\n\tPart b:')
orbitalElemGibbs(r1b, r2b, r3b)

# c.
r1c = [[   1,    0, 0]]
r2c = [[-0.8,  0.6, 0]]
r3c = [[ 0.8, -0.6, 0]]
print('\n\tPart c:')
orbitalElemGibbs(r1c, r2c, r3c)

# d. 
r1d = [[0.20709623, 3.53552813, 1.2071255]]
r2d = [[0.91420062, 4.9497417,  1.91423467]]
r3d = [[1.62130501, 6.36395526, 2.62134384]]
print('\n\tPart d:')
orbitalElemGibbs(r1d, r2d, r3d)

# e. 
r1e = [[1, 0, 0]]
r2e = [[0, 1, 0]]
r3e = [[-1, 0, 0]]
print('\n\tPart e:')
orbitalElemGibbs(r1e, r2e, r3e)

# f.
r1f = [[7, 2, 0]]
r2f = [[1, 1, 0]]
r3f = [[2, 7, 0]]
print('\n\tPart f:')
orbitalElemGibbs(r1f, r2f, r3f)

# g.
r1g = [[    0, 2.7, 0]]
r2g = [[ 2.97,   0, 0]]
r3g = [[-2.97,   0, 0]]
print('\n\tPart g:')
orbitalElemGibbs(r1g, r2g, r3g)