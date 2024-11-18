from _imports import *

# Spherical coordinates
theta = np.linspace(0, 2 * np.pi, 100)
phi = np.linspace(0, np.pi, 100)

theta, phi = np.meshgrid(theta, phi)

radius = 3443.922786 # earth radius nmi
nmi_ft = 2.092567257e7 / 3443.922786 # ft / nmi
ft_nmi = nmi_ft**-1 # nmi / ft


# Plot ellipse in ecliptic
theta_ell = np.linspace(0, 2 * np.pi, 100)
phi_ell = np.linspace(0, np.pi, 100)

# Determine orbital elements given r and v
mu = 1.406e16 * (ft_nmi**3) #earth gravitational parameter, ft^3/s^2
rvec = np.array([radius + 200, 0, 0]).T # nmi
vvec = np.array([1e4*ft_nmi,          3e4*ft_nmi, 0]).T # nmi/s
r = np.sqrt(rvec.dot(rvec)) # another way is np.linalg.norm(rvec) 
v = np.sqrt(vvec.dot(vvec))

hvec = np.cross(rvec, vvec) # angular momentum ft^2/s
h = np.sqrt(hvec.dot(hvec))

p = h**2 / mu # semi latus rectum, ft

e_spec = (v**2 / 2) - (mu/r) # specific mechanical energy, ft^2/s^2
evec = (1/mu) * (e_spec*rvec - (rvec.dot(vvec)*vvec))
e = np.sqrt(evec.dot(evec))

a = p / (1 - e**2)
c = a * e
b = np.sqrt(a**2 - c**2)
xe = a * np.cos(theta_ell)
ye = b * np.sin(theta_ell)
ze = 0 * theta_ell

print('Orbital Elements:\np = \t{:.2f} n.mi.\ne = \t{:.2f}\na = \t{:.2f} n.mi.\nc = \t{:.2f} n.mi.'.format(p, e, a, c))

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
rp = a * (1 - e)

# Plot
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(x, y, z, cmap='viridis', alpha=0.1)
ax.plot(xe_c-a, ye_c, ze_c)
ax.scatter(rp, 0, 0, c='k')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('3D Sphere')
ax.set_aspect('equal', adjustable='box')
plt.show()

