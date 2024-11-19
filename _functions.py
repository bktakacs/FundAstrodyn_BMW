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