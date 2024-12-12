
import sys, os
import numpy as np
from Dates import convertYearDoy2JulianDay
from math import atan2


def crossProd(a, b):
    res = np.zeros(3)
    res[0] = a[1]*b[2] - a[2]*b[1]
    res[1] = a[2]*b[0] - a[0]*b[2]
    res[2] = a[0]*b[1] - a[1]*b[0]

    return res


def rotate(v, angle, axis):
    vRes = np.ones(3)
    if(axis==3):
        vRes[0] =  v[0]*np.cos(angle) + v[1]*np.sin(angle)
        vRes[1] = -v[0]*np.sin(angle) + v[1]*np.cos(angle)
        vRes[2] =  v[2]
    elif(axis==2):
        vRes[0] =  v[0]*np.cos(angle) - v[2]*np.sin(angle)
        vRes[1] =  v[1]
        vRes[2] =  v[0]*np.sin(angle) + v[2]*np.cos(angle)
    else:
        vRes[0] =  v[0]
        vRes[1] =  v[1]*np.cos(angle) + v[2]*np.sin(angle)
        vRes[2] = -v[1]*np.sin(angle) + v[2]*np.cos(angle)

    return vRes


def modulo(a,b):
    return a%b


def findSun(Year, Doy, Sod):
    
    d2r = np.pi/180
    AU = 1.49597870e8

    fday = Sod/86400
    JDN = convertYearDoy2JulianDay(Year, Doy, Sod) - 2415020

    vl = modulo(279.696678 + 0.9856473354*JDN,360)
    gstr = modulo(279.690983 + 0.9856473354*JDN + 360*fday + 180,360)
    g = modulo(358.475845 + 0.985600267*JDN,360)*d2r
    
    slong = vl + (1.91946-0.004789*JDN/36525)*np.sin(g) + 0.020094*np.sin(2*g)
    obliq = (23.45229-0.0130125*JDN/36525)*d2r
    
    slp = (slong-0.005686)*d2r
    sind = np.sin(obliq)* np.sin(slp)
    cosd = np.sqrt(1-sind*sind)
    sdec = atan2(sind,cosd)/d2r
    
    sra = 180 - atan2(sind/cosd/np.tan(obliq),-np.cos(slp)/cosd)/d2r
    
    sunPosition = np.ones(3)
    sunPosition[0] = np.cos(sdec*d2r) * np.cos((sra)*d2r) * AU
    sunPosition[1] = np.cos(sdec*d2r) * np.sin((sra)*d2r) * AU
    sunPosition[2] = np.sin(sdec*d2r) * AU
    
    # Rotate from inertial to non inertial system (ECI to ECEF)
    sunPosition = rotate(sunPosition, gstr*d2r, 3)
    
    # print("SUNPOS", JDN, Sod, sdec, sra, slp, sind, cosd, sunPosition)

    return sunPosition

