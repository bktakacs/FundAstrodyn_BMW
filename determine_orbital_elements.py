from _imports import *

# Spherical coordinates
theta = np.linspace(0, 2 * np.pi, 100)
phi = np.linspace(0, np.pi, 100)

theta, phi = np.meshgrid(theta, phi)

radius = 3443.922786 # earth radius nmi
nmi_ft = 2.092567257e7 / 3443.922786 # ft / nmi
ft_nmi = nmi_ft**-1 # nmi / ft
sec_day = (60*60*24)**-1
sec_hour = (60*60)**-1


# Plot ellipse in ecliptic
theta_ell = np.linspace(0, 2 * np.pi, 100)
phi_ell = np.linspace(0, np.pi, 100)

# Determine orbital elements given r and v in IJK geocentric equatorial system
mu = 1.406e16 * (ft_nmi**3) #earth gravitational parameter, nmi^3/s^2
rvec = np.array([radius + 100, 0, 0]).T # nmi
vvec = np.array([0,          3e4*ft_nmi, 0]).T # nmi/s
r = np.sqrt(rvec.dot(rvec)) # another way is np.linalg.norm(rvec) 
v = np.sqrt(vvec.dot(vvec))

hvec = np.cross(rvec, vvec) # angular momentum nmi^2/s
h = np.sqrt(hvec.dot(hvec))

p = h**2 / mu # semi latus rectum, nmi

e_spec = (v**2 / 2) - (mu/r) # specific mechanical energy, nmi^2/s^2

evec = (1/mu) * (e_spec*rvec - (rvec.dot(vvec)*vvec)) # eccentricity
e = np.sqrt(evec.dot(evec))

a = p / (1 - e**2) # semi major axis, nmi
cvec = a * evec # c, distance from origin to focus, nmi
c = np.sqrt(cvec.dot(cvec))
b = np.sqrt(a**2 - c**2) # semi minor axis, nmi

t = 2 * np.pi / np.sqrt(mu) * a**(3/2) * sec_hour # period, hours

xe = cvec[0] + (a * np.cos(theta_ell)) # coordinates of ellipse
ye = cvec[1] + (b * np.sin(theta_ell))
ze = cvec[2] + (0 * theta_ell)

print('Orbital Elements:\np = \t{:.2f} n.mi.\ne = \t{:.2f}\na = \t{:.2f} n.mi.\nT = \t{:.2f} hours\nc = \t{:.2f} n.mi.'.format(p, e, a, t, c))

# Rotate ellipse
alpha = -30 * np.pi / 360 # degrees
# conversion matrix for rotation about y axis
conv_matrix = [[np.cos(alpha),  0, np.sin(alpha)],
               [0,              1,             0],
               [-np.sin(alpha), 0, np.cos(alpha)]]
orb_coord = [[xe, ye, ze]]
orb_coord_trans = np.matmul(conv_matrix, orb_coord)
xe_c = orb_coord_trans[:,0].T
ye_c = orb_coord_trans[:,1].T
ze_c = orb_coord_trans[:,2].T


# Convert sphere to cartesian
x = radius * np.sin(phi) * np.cos(theta)
y = radius * np.sin(phi) * np.sin(theta)
z = radius * np.cos(phi)

# Mark points of interest
rp = a * (1 - e)
ra = a * (1 + e)

# Plot
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(x, y, z, cmap='viridis', alpha=1, zorder=0)
ax.plot(xe_c, ye_c, ze_c, zorder=10)
ax.scatter(rp, 0, 0, c='k', zorder=10)
ax.scatter(-ra, 0, 0, c='k', zorder=10)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('3D Sphere')
ax.set_aspect('equal', adjustable='box')
plt.show()

