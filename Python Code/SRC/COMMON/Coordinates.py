import math
import numpy as np

# Ref.: ESA_GNSS-Book_TM-23_Vol_I.pdf Section B.1.2 (Appendix B)
def xyz2llh(x,y,z):
    # --- WGS84 constants
    a = 6378137.0
    f = 1.0 / 298.257223563
    # --- derived constants
    b = a - f*a
    e = math.sqrt(math.pow(a,2.0)-math.pow(b,2.0))/a
    clambda = math.atan2(y,x)
    p = math.sqrt(pow(x,2.0)+pow(y,2))
    h_old = 0.0
    # first guess with h=0 meters
    theta = math.atan2(z,p*(1.0-math.pow(e,2.0)))
    cs = math.cos(theta)
    sn = math.sin(theta)
    N = math.pow(a,2.0)/math.sqrt(math.pow(a*cs,2.0)+math.pow(b*sn,2.0))
    h = p/cs - N
    while abs(h-h_old) > 1.0e-6:
        h_old = h
        theta = math.atan2(z,p*(1.0-math.pow(e,2.0)*N/(N+h)))
        cs = math.cos(theta)
        sn = math.sin(theta)
        N = math.pow(a,2.0)/math.sqrt(math.pow(a*cs,2.0)+math.pow(b*sn,2.0))
        h = p/cs - N
    Rad2Deg = 180.0 / math.pi
    return clambda * Rad2Deg, theta * Rad2Deg, h

# Ref.: ESA_GNSS-Book_TM-23_Vol_I.pdf Section B.1.1 (Appendix B)
def llh2xyz(lon,lat,h):
    N = (6378137.0 / math.sqrt(1 - 0.0066943799901*(math.sin(math.radians(lat))**2)))

    X = (N+h)*(math.cos(math.radians(lat))*math.cos(math.radians(lon)))
    Y = (N+h)*(math.cos(math.radians(lat))*math.sin(math.radians(lon)))
    Z = ((1-0.0066943799901)*N + h)*(math.sin(math.radians(lat))) 

    return X,Y,Z

# From ECEF to ENU
def ecef2enu(lon, latit, Coordinates, refcoord):
    Long = math.radians(lon)
    lat = math.radians(latit)
    
    Rot = np.array([[-np.sin(Long), np.cos(Long), 0], \
            [-np.cos(Long)*np.sin(lat), -np.sin(Long)*np.sin(lat), np.cos(lat)], \
                [np.cos(Long)*np.cos(lat), np.sin(Long)*np.cos(lat), np.sin(lat)]])

    CoordinatesENU = np.dot(Rot, Coordinates-refcoord)

    return CoordinatesENU
