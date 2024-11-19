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
        data: list,
        initial_estimate: list,
        parametric_equation,
        partial_differential_matrix,
        ):
    """
    Doesn't work for now
    """
    # Correct formats:
    data = np.array(data).T
    initial_estimate = np.array(initial_estimate).T

    # Separate data into x and y
    x = np.array(data[:, 0]).T
    y = np.array(data[:, 1]).T

    # Calculate initial residuals
    ybar = np.array(parametric_equation(x))
    btild = y - ybar

    # Matrix Algebra
    ata = np.matmul(partial_differential_matrix.T, partial_differential_matrix)
    ata_inv = np.linalg.inv(ata)
    atb = np.matmul(partial_differential_matrix.T, btild)
    deltaz = np.matmul(ata_inv, atb)

    new_estimate = initial_estimate + deltaz

###############################################################################

def magnitude(r: list):
    return np.linalg.norm(r)

###############################################################################

def orbitDetermThreePosVec(
        r1: list, r2: list, r3: list, units: str = 'CANONUNITS'
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

    Returns
    -------
    Nothing, prints a statement describing the orbit.
    """

    # Check position vectors
    if r1.shape != (1, 3) or r2.shape != (1, 3) or r3.shape != (1, 3): 
        raise ValueError('Position vectors must be of shape (1, 3).')
    
    if (
        (type(r1) != np.ndarray) or (type(r2) != np.ndarray) or              (type(r3) != np.ndarray)
    ):
        r1 = np.array(r1); r2 = np.array(r2); r3 = np.array(r3)
    
    # Test that position vectors are coplanar:
    coplan_test = np.dot(r1, (np.cross(r2, r3)).T)
    if coplan_test != 0:
        raise ValueError('Position vectors must be coplanar.')

    # Calculate vector magnitudes
    r1m = magnitude(r1)
    r2m = magnitude(r2)
    r3m = magnitude(r3)

    # Find D, N, S vectors
    d = np.cross(r1, r2) + np.cross(r2, r3) + np.cross(r3, r1)
    n = r3m * np.cross(r1, r2) + r1m * np.cross(r2, r3) + r2m * np.cross(r3, r1)
    s = (r2m - r3m) * r1 + (r3m - r1m) * r2 + (r1m - r2m) * r3

    # Test D != 0, N != 0, and D dot N > 0 to assure vectors decsribe possible two-body orbit
    dtest = magnitude(d); ntest = magnitude(n); dntest = np.dot(d, n.T)
    if dtest == 0 or ntest == 0 or dntest <= 0:
        raise ValueError('Vectors do not describe possible two-body orbit.')
    
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
    p = ntest / dtest# semi latus rectum, L
    e = magnitude(s) / dtest # eccentricity
    hm = np.sqrt(ntest * mu / dtest) # angular momentum magnitude L^2/T^2
    a = p / (1 - e**2) # semi major axis, L
    tperiod = 2 * np.pi / np.sqrt(mu) * a**(3/2) # period, S

    # Calculate perifocal basis vectors in IJK
    qvec = s / magnitude(s)
    wvec = n / magnitude(n)
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
        'semi-major axis   a = {:.2f} {} ; period       T = {:.2f} {}'
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
            a, aunit, tperiod, tunit
        )
    )

###############################################################################
