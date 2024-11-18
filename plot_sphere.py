from _imports import *

# Spherical coordinates
theta = np.linspace(0, 2 * np.pi, 100)
phi = np.linspace(0, np.pi, 100)

theta, phi = np.meshgrid(theta, phi)

radius = 3444


# Plot ellipse in ecliptic
theta_ell = np.linspace(0, 2 * np.pi, 100)
phi_ell = np.linspace(0, np.pi, 100)
# a = 3
# b = 2
# e = np.sqrt(1 - (b / a)**2)
# c = a * e
# Given e and p:
p = 700 + radius
e = 0.4
a = p / (1 - e**2)
c = a * e
b = np.sqrt(a**2 - c**2)
xe = a * np.cos(theta_ell)
ye = b * np.sin(theta_ell)
ze = 0 * theta_ell
print(a)

# Rotate ellipse
alpha = 0#-30 * np.pi / 360 # degrees
# conversion matrix for rotation about y axis
conv_matrix = [[np.cos(alpha),  0, np.sin(alpha)],
               [0,              1,             0],
               [-np.sin(alpha), 0, np.cos(alpha)]]
orb_coord = [[xe, ye, ze]]
orb_coord_trans = np.matmul(conv_matrix, orb_coord)
xe_c = orb_coord_trans[:,0].T
ye_c = orb_coord_trans[:,1].T
ze_c = orb_coord_trans[:,2].T


# Convert to cartesian
x = radius * np.sin(phi) * np.cos(theta)
y = radius * np.sin(phi) * np.sin(theta)
z = radius * np.cos(phi)

# Mark points of interest
rp = (1 + e) * (p / (1 + e))
ra = (1 - e) * (p / (1 - e))
print(rp, np.max(xe_c), np.max(xe_c) - rp)
print(ra)

# Plot
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(x+c, y, z, cmap='viridis', alpha=0.1)
ax.plot(xe_c, ye_c, ze_c)
ax.scatter(rp, 0, 0, c='k')
ax.scatter(-ra, 0, 0, c='k')
ax.scatter(np.max(xe_c), 0, 0, c='orange')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('3D Sphere')
ax.set_aspect('equal', adjustable='box')
plt.show()

