from _imports import *

def coordTransformMatrix(
        lat: float = 0.0, lon: float = 0.0, units: str = 'd',to_ijk: bool = True
    ):
    """
    Calculate a coordinate transformation matrix from either IJK -> SEZ or viceversa, given the latittue and longitude of the site (topos).
    """
    # Convert degrees to radians for np.sin/cos function. 
    if units == 'd':
        lat *= np.pi / 180
        lon *= np.pi / 180

    m = np.array([
        [np.sin(lat) * np.cos(lon),  -np.sin(lon),  np.cos(lat) * np.cos(lon)],
        [np.sin(lat) * np.sin(lon),   np.cos(lon),  np.cos(lat) * np.sin(lon)],
        [             -np.cos(lat),            0,   np.sin(lat)              ]
    ])

    m = m.T if to_ijk == False else m

    return m

###############################################################################

def differentialCorrectionAlgorithm(
        data: np.ndarray,
        estimate: np.ndarray,
        par_eq,
        partial_diff_eqs: list,
        ):
    """
    This function takes in orbital data (observations) and an intitial estimate of coefficients related to how x and y are related. We then use partial differentitation of the x-y relationship as well as the residuals of y to calculate an updated y. 

    Parameters:
    -----------
    data : ndarray
    estimate : ndarray
    par_eq : callable
    partial_diff_eqs : list of callables

    Returns:
    --------

    """

    n = len(data) # number of observations
    p = len(estimate) # number of elements


    # Separate data into x and y
    x = data[:, 0]
    y = data[:, 1]
    step = 0

    # while eval(condition) and step < 1:
    while step < 20:

        step += 1

        ybar = par_eq(x, estimate)
        b = y - ybar
        rms = np.sqrt(np.sum([x**2 for x in b]) / n)

        # Matrix Algebra
        # Build ~A matrix
        a = np.zeros((n, p))
        for index, val in enumerate(x):
            for j in range(p):
                a[index][j] = partial_diff_eqs[j](val, estimate)

        ata = np.matmul(a.T, a)
        ata_inv = np.linalg.inv(ata)
        atb = np.matmul(a.T, b)

        deltaz = np.matmul(ata_inv, atb)
        estimate = estimate + deltaz


        ybar_n = par_eq(x, estimate)
        b_n = y - ybar_n
        rms_n = np.sqrt(np.sum([x**2 for x in b_n]) / n)

        if np.abs(rms_n - rms) < 0.001:
            break
        
    print('Final Parameter Estimate:', estimate, '\nAfter ', step, ' steps.')

###############################################################################

def magnitude(r: list):
    return np.linalg.norm(r)

###############################################################################

def orbitalElemGibbs(
        r1: list, r2: list, r3: list, units: str = 'CANONUNITS',
        angle_units: str = 'DEG'
):
    """
    Function which determines orbital parameters given 3 position vectors.

    Parameters
    ----------
    r1 : list
        Position vector of the 1st observation
    r2 : list
        Position vector of the 2rd observation
    r3 : list
        3rd observation
    units : str
        Units of position vectors. Options are 'CANONUNITS', 'NM', or 'KM'.
    angle_units : str
        Units of inclination angle. Options are 'RAD', 'DEG'.

    Returns
    -------
    Nothing, prints a statement describing the orbit.
    """

    # Check position vectors
    if (
        (type(r1) != np.ndarray) or (type(r2) != np.ndarray) or              (type(r3) != np.ndarray)
    ):
        r1 = np.array(r1); r2 = np.array(r2); r3 = np.array(r3)
    if r1.shape != (1, 3) or r2.shape != (1, 3) or r3.shape != (1, 3): 
        raise ValueError('Position vectors must be of shape (1, 3).')
    
    
    # Test that position vectors are coplanar:
    coplan_test = np.dot(r1, (np.cross(r2, r3)).T)
    if coplan_test > 1e-6:  # add a little bit of tolerance
        raise ValueError('Position vectors must be coplanar.')
    
    # Calculate in rads or degs
    angle = 1 if angle_units == 'RAD' else (180 / np.pi)

    # Calculate vector magnitudes
    r1m = magnitude(r1)
    r2m = magnitude(r2)
    r3m = magnitude(r3)

    # Find D, N, S vectors
    d = np.cross(r1, r2) + np.cross(r2, r3) + np.cross(r3, r1)
    n = r3m * np.cross(r1, r2) + r1m * np.cross(r2, r3) + r2m * np.cross(r3, r1)
    s = (r2m - r3m) * r1 + (r3m - r1m) * r2 + (r1m - r2m) * r3

    # Test D != 0, N != 0, and D dot N > 0 to assure vectors decsribe possible two-body orbit
    print_statement = True
    dtest = magnitude(d); ntest = magnitude(n); dntest = np.dot(d, n.T)
    if dtest == 0 or ntest == 0 or dntest <= 0:
        # raise ValueError('Vectors do not describe possible two-body orbit.')
        print_statement = False
        print('WARNING: Vectors do not describe possible two-body orbit.')
        
    
    # Calculate B = D cross r
    b1 = np.cross(d, r1); b2 = np.cross(d, r2); b3 = np.cross(d, r3)
    # Calculate L based on units
    if units == 'CANONUNITS':
        mu = 1 # canonical units DU^3 / TU^2
    elif units == 'NM':
        mu = 1.407646882e16 * (ft_nmi**3) # nmi^3 / s^2
    elif units == 'KM':
        mu = 3.986012e5 # km^3 / s^2

    l = np.sqrt(mu / dtest / ntest)

    # Calculate velocity vectors and their magnitudes
    v1 = (l / r1m) * b1 + l * s # velocity L/T
    v2 = (l / r2m) * b2 + l * s
    v3 = (l / r3m) * b3 + l * s
    v1m = magnitude(v1); v2m = magnitude(v2); v3m = magnitude(v3)

    # Calculate orbital parameters
    e = magnitude(s) / dtest # eccentricity
    p = ntest / dtest# semi latus rectum, L
    h = np.cross(r1, v1)
    hm = np.sqrt(ntest * mu / dtest) # angular momentum magnitude L^2/T^2
    inc = np.arccos(h[0][0] / hm) * angle # inclination, degrees
    a = np.inf if e == 1 else p / (1 - e**2) # semi major axis, L
    tperiod = 2 * np.pi / np.sqrt(mu) * a**(3/2) # period, S

    # Calculate perifocal basis vectors in IJK
    # Test for zero vectors
    zero_vec = np.zeros((1, 3))
    qvec = zero_vec if (s == zero_vec).all() else s / magnitude(s)
    wvec = zero_vec if (n == zero_vec).all() else n / magnitude(n)
    pvec = np.cross(qvec, wvec)

    # Convert r, v, pvec, qvec, wvec to strings and round to 2 places for printing
    r1s = np.round(r1, 2).astype(str); r2s = np.round(r2, 2).astype(str)
    r3s = np.round(r3, 2).astype(str)
    v1s = np.round(v1, 2).astype(str); v2s = np.round(v2, 2).astype(str)
    v3s = np.round(v3, 2).astype(str)
    ps = np.round(pvec, 2).astype(str); qs = np.round(qvec, 2).astype(str)
    ws = np.round(wvec, 2).astype(str)

    # Units
    if units == 'CANONUNITS':
        runit = 'DU'; vunit = 'DU/TU'; punit = 'DU'; aunit = 'DU'; tunit = 'TU'
    elif units == 'NM':
        runit = 'n.mi.'; vunit = 'n.mi./s'; punit = 'n.mi.'; aunit = 'n.mi.'; tunit = 'hours'; tperiod *= sec_hour

    w = 5 if units == 'CANONUNITS' else 9
    if print_statement:
        print(
            '\nOrbital Parameters given the position vectors:\n'
            '  r1 = [{}]  ;  r2 = [{}]  ;  r3 = [{}] {}\n'
            '       |{}|          |{}|          |{}|\n'
            '       [{}]          [{}]          [{}]\n'
            '|r1| =  {}    |r2| =  {}    |r3| = {}\n\n'
            'Perifocal basis vectors in IJK system:\n'
            '   P = [{}]  ;   Q = [{}]  ;   W = [{}]\n'
            '       |{}|          |{}|          |{}|\n'
            '       [{}]          [{}]          [{}]\n\n'
            'Velocity vectors at each position:\n'
            '  v1 = [{}]  ;  v2 = [{}]  ;  v3 = [{}] {}\n'
            '       |{}|          |{}|          |{}|\n'
            '       [{}]          [{}]          [{}]\n'
            '|v1| =  {}    |v2| =  {}    |v3| = {}\n\n'
            'Orbital elements:\n'
            'semi-latus pectum p = {:.2f} {} ; eccentricity e = {:.2f}\n'
            'semi-major axis   a = {:.2f} {} ; period       T = {:.2f} {}\n'
            'inclination       i = {:.2f} {}'
            ''.format(
                # r1, r2, and r3 and their units and magnitudes
                str.center(r1s[0][0], w), str.center(r2s[0][0], w),
                str.center(r3s[0][0], w), runit,
                str.center(r1s[0][1], w),
                str.center(r2s[0][1], w), str.center(r3s[0][1], w),
                str.center(r1s[0][2], w), str.center(r2s[0][2], w),
                str.center(r3s[0][2], w),
                str.center(str(np.round(r1m, 2)), w),
                str.center(str(np.round(r2m, 2)), w),
                str.center(str(np.round(r3m, 2)), w),
                # perifocal basis vectors p, q, w
                str.center(ps[0][0], w), str.center(qs[0][0], w),
                str.center(ws[0][0], w), str.center(ps[0][1], w),
                str.center(qs[0][1], w), str.center(ws[0][1], w),
                str.center(ps[0][2], w), str.center(qs[0][2], w),
                str.center(ws[0][2], w),
                # v1, v2, v3
                str.center(v1s[0][0], w), str.center(v2s[0][0], w),
                str.center(v3s[0][0], w), vunit,
                str.center(v1s[0][1], w),
                str.center(v2s[0][1], w), str.center(v3s[0][1], w),
                str.center(v1s[0][2], w), str.center(v2s[0][2], w),
                str.center(v3s[0][2], w),
                str.center(str(np.round(v1m, 2)), w),
                str.center(str(np.round(v2m, 2)), w),
                str.center(str(np.round(v3m, 2)), w),
                # orbital elements and their units
                p, punit, e,
                a, aunit, tperiod, tunit,
                inc, angle_units
            )
        )

###############################################################################

def sign(x):
    return '+' if x >= 0 else '-'

###############################################################################

def kepler_problem(
        r0: list, v0: list, t: float, t0: float = 0.0, units: str = 'CANON', tolerance: float = 1e-4, w: int = 3
):
    """
    Solves Kepler's problem for a given initial position vector, velocity vector, time at which to find r and v.

    Parameters:
    r0 (list): Initial position vector [x, y, z] in the IJK system.
    v0 (list): Initial velocity vector [vx, vy, vz] in the IJK system.
    t (float): Time at which to find r and v.
    t0 (float): Initial time. Default is 0.0.
    units (str): Unit system.
    tolerance (float): Tolerance for convergence in Newton iteration scheme
    w (int): Number of decimal places to round results to

    Returns:
    r (list): Position vector [x, y, z] in the IJK system at time t.
    v (list): Velocity vector [vx, vy, vz] in the IJK system at time t.
    """

    # Make r0 and v0 into np.array
    if type(r0) != np.ndarray or type(v0) != np.ndarray:
        r0 = np.array(r0)
        v0 = np.array(v0)

    # Get correct units for m, and create sqrt(m) var
    mu = 1 if units == 'CANON' else mu_nmi_s if units == 'NMI' else mu_ft_s
    mur = np.sqrt(mu)

    # Derive r0, v0 magnitudes
    r0m = magnitude(r0)
    v0m = magnitude(v0)

    # Step 1: Fro`m r0, v0, determine r0m and a -------------------------------

    # Orbital constants
    h = np.dot(r0, v0.T) # not angular momentum, just saving space
    sme =v0m**2 / 2 - mu / r0m # specific mechanical energy, DU^2/TU^2
    a = - mu / (2 * sme) # semi-major axis, DU

    # Step 2: Given (t - t0), solve universal time-of-flight (TOF) equation for x using Newton iteration scheme -------------------------------------------

    # Define C(z), S(z)
    C = lambda z: (1 - np.cos(np.sqrt(z))) / 2
    S = lambda z: (np.sqrt(z) - np.sin(np.sqrt(z))) / np.sqrt(z**3)

    # Newton Iteration Scheme -------------------------------------------------

    # First guess for x:
    xn = mur * (t - t0) / a

    # Define variable z
    z = lambda x: x**2 / a

    # Define differential equaiton dt/dx
    dtdx = lambda x: x**2 * C(z(x)) / mur + h / mur * (1 - z(x) * S(z(x)) + r0m / mur * (1 - z(x) * C(z(x))))

    # Define t_n equation
    tn = lambda xn: h / mu * xn**2 * C(z(xn)) + (1 - r0m/a) * xn**3 * S(z(xn)) / mur + r0m * xn / mur

    # Iterate until convergence (t ~ tn)
    step = 1
    while True:
        tn1 = tn(xn)
        xn1 = xn + (t - tn1) / dtdx(xn)

        tn2 = tn(xn1)

        if np.abs(t - tn2) < 1e-4:
            break
        else:
            tn1 = tn2
            xn = xn1
            step += 1
        
        if step > 100:
            print('Newton iteration did not converge after 100 steps.')
            break

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
    print(v)
    vm = magnitude(v)

    # Output results
    print('\n'
        'Final Position Vector r(t={}) = {}{} I {} {} J {} {} K\n'
        'Final Velocity Vector v(t={}) = {}{} I {} {} J {} {} K\n'
        ''.format(
            t, sign(r[0]), np.round(np.abs(r[0]), w), sign(r[1]), np.round(np.abs(r[1]), w), sign(r[2]), np.round(np.abs(r[2]), w),
            t, sign(v[0]), np.round(np.abs(v[0]), w), sign(v[1]), np.round(np.abs(v[1]), w), sign(v[2]), np.round(np.abs(v[2]), w),
        ))

###############################################################################