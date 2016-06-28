import numpy as np
from math import pi, sin, cos, sqrt
from obspy.geodetics import gps2dist_azimuth

def wgs84():

    # semi-major axis, in m
    a = 6378137.0

    # semi-minor axis, in m
    b = 6356752.314245

    # inverse flattening f
    f = a/(a-b)

    # squared eccentricity e
    e_2 = (a**2-b**2)/a**2
    
    return(a,b,e_2,f)

def len_deg_lon(lat):
    
    (a,b,e_2,f) = wgs84()
    # This is the length of one degree of longitude 
    # approx. after WGS84, at latitude lat
    # in m
    lat = pi/180*lat
    dlon = (pi*a*cos(lat))/180*sqrt((1-e_2*sin(lat)**2))
    return round(dlon,2)

def len_deg_lat(lat):
    # This is the length of one degree of latitude 
    # approx. after WGS84, between lat-0.5deg and lat+0.5 deg
    # in m
    lat = pi/180*lat
    dlat = 111132.954 - 559.822 * cos(2*lat) + 1.175*cos(4*lat)
    return round(dlat,2)


#ToDo: Tests
def points_on_sphere(dx,xmin=-180.,xmax=180.,ymin=-90.,ymax=90.,c_centr=None,\
radius=None):
    """
    Calculate a more or less equally spaced grid on spherical Earth's surface.
    :param dx: spacing in latitudinal and longitudinal direction in meter
    :type c_centr: Tuple
    :param c_centr: Specify a central location
    :type radius: float
    :param radius: Radius around central location in m; no sources beyond this will be included
    :returns: np.array(latitude, longitude) of grid points, where -180<=lon<180     and -90 <= lat < 90
    """
    
    if xmax <= xmin or ymax <= ymin:
        msg = 'Lower bounds must be lower than upper bounds.'
        raise ValueError(msg)

    
    gridx = []
    gridy = []
    
    lat = ymin
    
    while lat <= ymax:
        d_lat = dx / len_deg_lat(lat)
        lon = xmin
        while lon <= xmax:
            
            if c_centr and radius:
                if gps2dist_azimuth(lat,lon,c_centr[0],c_centr[1])[0] > radius:
                    if abs(lat) != 90.:
                        d_lon = dx / len_deg_lon(lat)
                        lon += d_lon
                        continue
                    else:
                        break
                    
            gridx.append(lon)
            gridy.append(lat)
            
            if abs(lat) == 90:
                # length of a degree longitude will be 0.
                break
            else:
                d_lon = dx / len_deg_lon(lat)
                lon += d_lon
        lat += d_lat # do not start at pole or zero division will raise...
        
            
    # return values sorted by longitude, because basemap complains otherwise.
    grid = list(zip(*sorted(zip(gridx, gridy), key=lambda it: it[0])))
    return grid
    