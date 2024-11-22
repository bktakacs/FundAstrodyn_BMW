from _imports import *
from _functions import *

# Here we implement the algorithm of computing residuals to improve a preliminary orbit found on pg 125.

# n   number of observations
# p   numer of elements (parameters in the equation)
# ~W  nxn diagonal matrix representing square of confidence
# alpha, beta     elements (only 2 here, as p = 2 in our example)
# ~A  nxp matrix of partial derivatives of each quantity with respenct to each 
#     element measured
# ~b  nx1 matrix fo residuals based on previous estimate of elements
# ~∆z px1 matrix of computed corrections to estimates of elements
# x_i independent variable measurement
# y_i dependent variable measurement corresponding to x_i
# y_i_bar     computed value of dependent variable using previous values for 
#             elements

# if p > n, not enough information to solve problem
# if p = n, there's an exact solution
# if p < n, more equations than unknowns, infinite solutions, no unique one. Find best solution in the "least-squares" sense.

# ~∆z = (~A.T ~W ~A)^-1 ~A.T ~W ~b

# Algorithm:
# 1. Solve eq. 1 for changes to elements.
# 2. Correct elements (alpha_new = alpha_old + ∆alpha)
# 3. Compute new residuals using same data with new elements
# 4. Repeat steps 1-3 until residuals are:
#    a. zero, for the exactly determined case (p = n)
#    b. minimum, as described below. 

# Example problem: Let there be two elements, a & b, with the relationship y = a + bx, and two observations (n=p=2). Assume equal confidence in data. Choose a=2 and b=3 as first estimates.
elem = np.array([[2, 3]]).T
# Given:
# Observation      xi       yi      yi_bar      Residual
#      1           2        1         8          -7
#      2           3        2        11          -9
x = np.array([[2, 3]]).T
y = np.array([[1, 2]]).T
yeq = lambda x: elem[0][0] + elem[1][0]*x
ybar = np.array(yeq(x))

# yi_bar was predicted as follows:
# x1 = 2, y1_bar = 2 + 3x1 = 8
# x2 = 3, y2_bar = 2 + 3x2 = 11

# Calculate residuals:
# ~b = [y1 - y1_bar] = [-7] , ~W = [1 0]
#      [y2 - y2_bar] = [-9] ,      [0 1]
btild = y - ybar


# From fitting equation y = a + bx, the partial derivatives are:
# dyi/da = 1, dyi/db = xi

# Then ~A is:
# ~A = [dy1/da  dy1/db] = [1 2]
#      [dy2/da  dy2/db] = [1 3]
atild = np.array([[1, 2], [1, 3]])

# Now solve for ~∆z:
# ~A.T ~A = 
ata = np.matmul(atild.T, atild)
# (~A.T ~A)^-1 = 
ata_inv = np.linalg.inv(ata)
# ~A.T ~b = 
atb = np.matmul(atild.T, btild)
# ~∆z = 
deltaztild = np.matmul(ata_inv, atb)

# Update fitting equation: y = 
elem = elem + deltaztild
ybar_new = np.array(yeq(x))
# Calculate residuals (off by 1e-14 for some reason, maybe from calculating inverse)
btild_new = y - ybar_new
print(btild_new)

# Alternatively, for the example problem on page 128:
data = np.array([[1, 2.5], [2, 8], [3, 19], [4, 50]])
estimate = np.array([0.474, 3.360])

def par_eq_1(x, estimate):
    return estimate[0] * x**(estimate[1])

def par_dif_1(x, estimate):
    return x**(estimate[1])

def par_diff_2(x, estimate):
    return estimate[0] * x**(estimate[1]) * np.log(x)

partial_diff_eqs = [par_dif_1, par_diff_2]

differentialCorrectionAlgorithm(data, estimate, par_eq_1, partial_diff_eqs)