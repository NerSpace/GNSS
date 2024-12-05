import numpy as np
from COMMON.Coordinates import llh2xyz
from COMMON import GnssConstants as Const

def buildGmatrix(G, SatCoord, i, Constel, RcvrPosClk):
    
    '''
    Purpose: Construct the geometry matrix (G) for satellite positioning

    Parameters
    ==========
    G: ndarray
        Geometry matrix to be filled with satellite position data.
    SatCoord: list or ndarray
        Coordinates of the satellite in the format [East, North, Up].
    i: int
        Index of the current satellite to fill in the geometry matrix.
    Constel: str
        Satellite constellation type ('G' for GPS, 'E' for Galileo).
    RcvrPosClk: list or ndarray
        Receiver position and clock bias.

    Returns
    =======
    G: ndarray
        Updated geometry matrix with the unitary vectors and clock bias included.

    '''

    # Satellite position
    SatX = SatCoord[0]
    SatY = SatCoord[1]
    SatZ = SatCoord[2]

    # Rcvr position
    LeoX = RcvrPosClk[0]
    LeoY = RcvrPosClk[1]
    LeoZ = RcvrPosClk[2]

    # Range Vector and geom range
    r = np.array([SatX-LeoX, SatY-LeoY, SatZ-LeoZ])
    rho = np.linalg.norm(r)

    # Fill G matrix with unitary vectors
    G[i, 0] = - r[0] / rho
    G[i, 1] = - r[1] / rho  
    G[i, 2] = - r[2] / rho

    # Clock bias for multi-constellation    
    if Constel == "G":
        G[i, 3] = 1 # GPS clock bias column
    elif Constel == "E":
        G[i, 4] = 1 # GAL clock bias column
    return G

def buildGmatrixEl(G, SatCorr, i, Constel):

    '''
    Purpose: Construct the geometry matrix (G) using satellite azimuth and elevation

    Parameters
    ==========
    G: ndarray
        Geometry matrix to be filled with satellite direction components.
    SatCorr: dict
        Dictionary containing satellite correction data, including azimuth and elevation.
    i: int
        Index of the current satellite to fill in the geometry matrix.
    Constel: str
        Satellite constellation type ('G' for GPS, 'E' for Galileo).

    Returns
    =======
    G: ndarray
        Updated geometry matrix with the unitary vectors and clock bias included.

    '''

    # Extract azimuth and elevation
    Azim = np.deg2rad(SatCorr["Azimuth"])
    Elev = np.deg2rad(SatCorr["Elevation"])

    # Compute east, north and up components
    East = -np.cos(Elev) * np.sin(Azim)
    North = - np.cos(Elev) * np.cos(Azim)
    Up = -np.sin(Elev)

    # Fill G matrix with unitary vectors
    G[i, 0] = East
    G[i, 1] = North
    G[i, 2] = Up

    # Clock bias for multi-constellation    
    if Constel == "G":
        G[i, 3] = 1 # GPS clock bias column
    elif Constel == "E":
        G[i, 4] = 1 # GAL clock bias column
    return G

def buildRes(SatCorr, i, Constel, RcvrPosClk, D_rho):

    '''
    Purpose: Compute the residual for satellite positioning

    Parameters
    ==========
    SatCorr: dict
        Dictionary containing satellite correction data, including satellite positions and code corrections.
    i: int
        Index of the current satellite for which the residual is being calculated.
    Constel: str
        Satellite constellation type ('G' for GPS, 'E' for Galileo).
    RcvrPosClk: list or ndarray
        Receiver position and clock bias in the format [X, Y, Z, ClockBiasG, ClockBiasE].
    D_rho: ndarray
        Array to store the calculated residuals for each satellite.

    Returns
    =======
    D_rho: ndarray
        Updated array containing the residuals for the specified satellite.

    '''

    # Satellite position
    SatX = SatCorr["SatX"]
    SatY = SatCorr["SatY"]
    SatZ = SatCorr["SatZ"]

    # Rcvr position
    LeoX = RcvrPosClk[0]
    LeoY = RcvrPosClk[1]
    LeoZ = RcvrPosClk[2]

    # Range Vector and geom range
    r = np.array([SatX-LeoX, SatY-LeoY, SatZ-LeoZ])
    rho = np.linalg.norm(r)

    # Clock bias for multi-constellation    
    if Constel == "G":
        D_rho[i] = SatCorr["CorrCode"] - rho - RcvrPosClk[3]
    elif Constel == "E":
        D_rho[i] = SatCorr["CorrCode"] - rho - RcvrPosClk[4]
    
    return D_rho

def buildWmatrix(W, sigma, i):

    '''
    Purpose: Construct the weight matrix for satellite measurements

    Parameters
    ==========
    W: ndarray
        Weight matrix to be filled with the inverse of the variances of the satellite measurements.
    sigma: float
        Standard deviation of the measurement for the i-th satellite.
    i: int
        Index of the current satellite for which the weight is being set.

    Returns
    =======
    W: ndarray
        Updated weight matrix with the weight for the specified satellite.

    '''

    W[i, i] = 1 / sigma**2  # Set the weight for the i-th satellite

    return W

def computeDOPs(G, RcvrRefPosLlh):

    '''
    Purpose: Compute the Dilution of Precision (DOP) values for the given geometry matrix.

    Parameters
    ==========
    G: ndarray
        Geometry matrix containing satellite position information in relation to the receiver.
    RcvrRefPosLlh: list
        Reference position of the receiver in Latitude, Longitude, and Height (LLH) coordinates.

    Returns
    =======
    PDOP: float
        Position Dilution of Precision, indicating the geometric quality of satellite positions.
    HDOP: float
        Horizontal Dilution of Precision, representing the horizontal positioning accuracy.
    VDOP: float
        Vertical Dilution of Precision, representing the vertical positioning accuracy.

    '''

    # Compute Q matrix
    Q = np.linalg.inv(np.dot(G.T, G))

    q_ee = Q[0, 0]
    q_nn = Q[1, 1]
    q_uu = Q[2, 2]

    # Calculate HDOP and VDOP
    HDOP = np.sqrt(q_ee + q_nn)
    VDOP = np.sqrt(q_uu)

    # Calculate PDOP
    PDOP = np.sqrt(q_ee+q_nn+q_uu)

    return PDOP, HDOP, VDOP


def buildSmatrix(G, W):

    '''
    Purpose: Construct the S matrix used in the Weighted Least Squares (WLS) adjustment.

    Parameters
    ==========
    G: ndarray
        Geometry matrix that represents the relationship between satellite measurements and receiver position.
    W: ndarray
        Weight matrix that contains weights for the satellite measurements based on their accuracies.

    Returns
    =======
    S: ndarray
        The constructed S matrix, which is used in the adjustment process.

    '''

    G_T = G.T
    GT_W = np.dot(G_T, W)
    GTW_G = np.dot(GT_W, G)
    GTWG_inv = np.linalg.inv(GTW_G)
    GTWG_GT = np.dot(GTWG_inv, G_T)

    S = np.dot(GTWG_GT, W)

    return S

def runWlsqIteration(CorrInfo, S, RcvrPosClk, NumIter):
    
    '''
    Purpose: Execute one iteration of the Weighted Least Squares (WLS) adjustment process to refine the receiver position.

    Parameters
    ==========
    CorrInfo: dict
        Dictionary containing the valid corrected satellite data.
    S: ndarray
        S matrix used for the adjustment calculations.
    RcvrPosClk: list
        Current estimated position and clock bias of the receiver in the format [Longitude, Latitude, Altitude, ClkGPS, ClkGAL].
    NumIter: int
        Current iteration number for the WLS process.

    Returns
    =======
    NumIter: int
        Updated iteration count after this adjustment.
    RcvrPosClkDelta: ndarray
        Change in the receiver position and clock biases after the iteration.
    RcvrPosClk: list
        Updated receiver position and clock biases.

    '''

    # Initialize variables
    idx=0
    NSats = len(CorrInfo)
    D_Rho = np.zeros(NSats)

    # LEO Coordinates in ECEF
    LEO_ECEF = llh2xyz(RcvrPosClk[0], RcvrPosClk[1], RcvrPosClk[2])
    
    RcvrPosClkECEF = np.array([LEO_ECEF[0], LEO_ECEF[1], LEO_ECEF[2], RcvrPosClk[3], RcvrPosClk[4]])

    # Compute Residuals
    for SatLabel, SatCorr in CorrInfo.items():
        # Get Constel
        Constel = SatLabel[0]
        D_Rho = buildRes(SatCorr, idx, Constel, RcvrPosClkECEF, D_Rho)
        idx += 1

    # Compute RcvrPosClkDelta in ENU Coord
    RcvrPosClkDelta = np.dot(S, D_Rho)

    DeltaLong = np.degrees(RcvrPosClkDelta[0]/(Const.EARTH_RADIUS*np.cos(np.radians(RcvrPosClk[1]))))
    DeltaLat= np.degrees(RcvrPosClkDelta[1]/(Const.EARTH_RADIUS))

    # New RcvrPosClk
    RcvrPosClk[0] = RcvrPosClk[0] + DeltaLong
    RcvrPosClk[1] = RcvrPosClk[1] + DeltaLat
    RcvrPosClk[2] = RcvrPosClk[2] + RcvrPosClkDelta[2]
    RcvrPosClk[3] = RcvrPosClk[3] + RcvrPosClkDelta[3]
    RcvrPosClk[4] = RcvrPosClk[4] + RcvrPosClkDelta[4]

    NumIter += 1

    return NumIter, RcvrPosClkDelta, RcvrPosClk
