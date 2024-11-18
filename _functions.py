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
    else:
        pass

    m = np.array([
        [np.sin(lat) * np.cos(lon),  -np.sin(lon),  np.cos(lat) * np.cos(lon)],
        [np.sin(lat) * np.sin(lon),   np.cos(lon),  np.cos(lat) * np.sin(lon)],
        [             -np.cos(lat),            0,   np.sin(lat)              ]
    ])

    m = m.T if to_ijk == False else m

    return m
