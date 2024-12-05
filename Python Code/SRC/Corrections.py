#!/usr/bin/env python

########################################################################
# SENTUS-V3/SRC/Corrections.py:
# This is the Corrections Module of SENTUS tool
#
#  Project:        SENTUS
#  File:           Corrections.py
#  Date(YY/MM/DD): 14/07/24
#
#   Author: Nerea Sánchez
#   Copyright 2024 GNSS Academy / Nerea Sánchez
#
########################################################################

# Import External and Internal functions and Libraries
#----------------------------------------------------------------------
import sys, os
# Add path to find all modules
Common = os.path.dirname(os.path.dirname(
    os.path.abspath(sys.argv[0]))) + '/COMMON'
sys.path.insert(0, Common)
from collections import OrderedDict
from COMMON import GnssConstants as Const
from COMMON.Misc import findSun, crossProd
from COMMON.Coordinates import xyz2llh
from COMMON.Dates import convertYearDoy2JulianDay
import numpy as np
from InputOutput import ObsIdxC, ObsIdxP, LeoPosIdx, LeoQuatIdx
from InputOutput import SatPosIdx, SatApoIdx, SatClkIdx, SatBiaIdx


def runCorrectMeas(Year,
                   Doy,
                   ObsInfo,
                   Conf, 
                   PreproObsInfo, 
                   LeoPosInfo,
                   LeoQuatInfo,
                   SatPosInfo, 
                   SatApoInfo,
                   SatClkInfo,
                   SatBiaInfo,
                   SatComPos_1,
                   Sod_1
                   ):

    '''
    
    Purpose: correct GNSS preprocessed measurements and compute the first
             pseudo range residuals

             More in detail, this function handles the following:
             tasks:

                *  Compute the Satellite Antenna Phase Center position at the transmission time and corrected from the Sagnac
                   effect interpolating the SP3 file positions
                *  Compute the Satellite Clock Bias interpolating the biases coming from the RINEX CLK file and
                   applying the Relativistic Correction (DTR)
                *  Correct the Pre-processed measurements from Geometrical Range, Satellite clock and Troposphere. 
                *  Build the Corrected Measurements and Measurement Residuals
                *  Build the Sigma UERE


    Parameters
    ==========
    Conf: dict
            Configuration dictionary
    Rcvr: list
            Receiver information: position, masking angle...
    ObsInfo: list
            OBS info for current epoch
            ObsInfo[1][1] is the second field of the 
            second satellite
    PreproObsInfo: dict
            Preprocessed observations for current epoch per sat
            PreproObsInfo["G01"]["C1"]
    LeoPosInfo: dict
            containing the LEO reference positions
    LeoQuatInfo: dict
            containing the LEO quaternions
    SatPosInfo: dict
            containing the SP3 file info
    SatApoInfo: dict
            containing the ANTEX file info
    SatClkInfo: dict
            containing the RINEX CLK file info
    SatBiaInfo: dict
            containing the BIA file info
    SatComPos_1: dict
            containing the previous satellite positions
    Sod_1: dict
            containing the time stamp of previous satellite positions

    Returns
    =======
    CorrInfo: dict
            Corrected measurements for current epoch per sat
            CorrInfo["G01"]["CorrectedPsr"]

        '''

    # Initialize output
    CorrInfo = OrderedDict({})
    RcvrRefPosXyz = np.zeros(3)
    RcvrRefPosLlh = np.zeros(3)

    # Initialize accumulations for the weighted average
    Res_accum_GPS = 0
    UERE_accum_GPS = 0
    Res_accum_GAL = 0
    UERE_accum_GAL = 0

    # Get Sod and Receiver Reference Positio CoM
    Sod = int(ObsInfo[1][0][ObsIdxP["SOD"]])
    RcvrRefPosXyzCom = computeLeoComPos(Sod, LeoPosInfo)

    # Remove information from constelations not used
    Constel = Conf["NAV_SOLUTION"]
    SatPosInfo = filterByConstel(SatPosInfo, Constel)
    SatApoInfo = filterByConstel(SatApoInfo, Constel)
    SatClkInfo = filterByConstel(SatClkInfo, Constel)
    SatBiaInfo = filterByConstel(SatBiaInfo, Constel)
    PreproObsInfo = filterPreproObsByConstel(PreproObsInfo, Constel)

    # Loop over satellites
    for SatLabel, SatPrepro in PreproObsInfo.items():

        # Get constellation
        Constel = SatLabel[0]
        
        Wave = {}
        # Get wavelengths
        if Constel == 'G':

            # L1 wavelength
            Wave["F1"] = Const.GPS_L1_WAVE

            # L2 wavelength
            Wave["F2"] = Const.GPS_L2_WAVE

            # Gamma GPS
            GammaF1F2 = Const.GPS_GAMMA_L1L2

        elif Constel == 'E':

            # E1 wavelength
            Wave["F1"] = Const.GAL_E1_WAVE

            # E5a wavelength
            Wave["F2"] = Const.GAL_E5A_WAVE

            # Gamma Galileo
            GammaF1F2 = Const.GAL_GAMMA_E1E5A

        # Initialize output info
        SatCorrInfo = {
            "Sod": Sod,             # Second of day
            "Doy": Doy,               # Day of year
            "Elevation": SatPrepro["Elevation"],       # Elevation
            "Azimuth": SatPrepro["Azimuth"],         # Azimuth
            "Flag": 1,              # 0: Not Used 1: Used
            "LeoX": 0.0,            # X-Component of the Receiver CoP Position 
            "LeoY": 0.0,            # Y-Component of the Receiver CoP Position  
            "LeoZ": 0.0,            # Z-Component of the Receiver CoP Position  
            "LeoApoX": 0.0,         # X-Component of the Receiver APO in ECEF
            "LeoApoY": 0.0,         # Y-Component of the Receiver APO in ECEF
            "LeoApoZ": 0.0,         # Z-Component of the Receiver APO in ECEF
            "SatX": 0.0,            # X-Component of the Satellite CoP Position 
                                    # at transmission time and corrected from Sagnac
            "SatY": 0.0,            # Y-Component of the Satellite CoP Position  
                                    # at transmission time and corrected from Sagnac
            "SatZ": 0.0,            # Z-Component of the Satellite CoP Position  
                                    # at transmission time and corrected from Sagnac
            "SatApoX": 0.0,         # X-Component of the Satellite APO in ECEF
            "SatApoY": 0.0,         # Y-Component of the Satellite APO in ECEF
            "SatApoZ": 0.0,         # Z-Component of the Satellite APO in ECEF
            "ApoProj": 0.0,         # Projection of the Satellite APO
            "SatClk": 0.0,          # Satellite Clock Bias
            "SatCodeBia": 0.0,      # Satellite Code Bias
            "SatPhaseBia": 0.0,     # Satellite Phase Bias
            "FlightTime": 0.0,      # Signal Flight Time
            "Dtr": 0.0,             # Relativistic correction
            "CorrCode": 0.0,        # Code corrected from delays
            "CorrPhase": 0.0,       # Phase corrected from delays
            "GeomRange": 0.0,       # Geometrical Range (distance between Satellite 
                                    # Position and Receiver Reference Position)
            "CodeResidual": 0.0,    # Code Residual
            "PhaseResidual": 0.0,   # Phase Residual
            "RcvrClk": 0.0,         # Receiver Clock estimation
            "SigmaUere": 0.0,       # Sigma User Equivalent Range Error (Sigma of 
                                    # the total residual error associated to the 
                                    # satellite)

        } # End of SatCorrInfo

        CorrInfo[SatLabel] = SatCorrInfo

        # Only for those satellites with Status OK and there is previous information
        if SatPrepro["Status"] == 1 and Sod_1[SatLabel] != []:
            # Compute Satellite Clock Bias (linear interpolation between closer inputs)
            SatClkBias = computeSatClkBias(Sod, SatClkInfo, SatLabel)
            if SatClkBias == Const.NAN:
                SatCorrInfo["Flag"] = 0
                break
            SatClkBiam = SatClkBias*Const.SPEED_OF_LIGHT

            # Compute Delta t
            DeltaT = SatPrepro["C1"] / Const.SPEED_OF_LIGHT

            # Compute Transmission Time
            TransmissionTime = Sod - DeltaT - SatClkBias

            # Compute Receiver Position at Reception Time
            # Get the Receiver APO in ECEF coordinates
            RcvrApoXyz = computeRcvrApo(Conf, Year, Doy, Sod, SatLabel, LeoQuatInfo)

            # Apply the APO
            RcvrRefPosXyz = RcvrRefPosXyzCom + RcvrApoXyz
            RcvrLos = RcvrRefPosXyz / np.linalg.norm(RcvrRefPosXyz)

            # Compute Long-Lat-H
            RcvrRefPosLlh = xyz2llh(RcvrRefPosXyz[0], RcvrRefPosXyz[1], RcvrRefPosXyz[2])

            # Compute Satellite CoM Position at Transmission Time
            # 10-point Lagrange interpolation between closer inputs (SP3 positions)
            SatComPos1 = computeSatComPos(TransmissionTime, SatPosInfo, SatLabel)

            # Compute Flight Time
            distance = np.linalg.norm(SatComPos1 - RcvrRefPosXyz)
            FlightTime = distance / Const.SPEED_OF_LIGHT 

            # Apply Sagnac correction
            SatComPos = applySagnac(SatComPos1, FlightTime)

            # Compute APO in ECEF from ANTEX APOs in SFR
            SunPos = findSun(Year, Doy, Sod)
            SatApo = computeSatApo(SatLabel, SatComPos, SunPos, SatApoInfo, GammaF1F2)

            # Apply APOs to the Satellite Position
            SatCopPos = SatComPos + SatApo

            # Get Satellite Biases in meters
            [CodeSatBias, PhaseSatBias] = getSatBias(GammaF1F2, SatLabel, SatBiaInfo)

            # Compute DTR (Relativistic correction)
            Dtr = computeDtr(SatComPos_1, SatComPos, Sod, Sod_1, SatLabel)
            
            if Dtr != Const.NAN:
                # Apply Dtr to Clock Bias if Dtr is computed
                SatClkBiam = SatClkBiam + Dtr 
                # Corrected Measurements from previous information
                CorrCode = SatPrepro["SmoothIF"] + SatClkBiam + CodeSatBias
                CorrPhase = SatPrepro["IF_P"] + SatClkBiam + PhaseSatBias
                # Compute the Geometrical Range
                GeomRange = ComputeGeomRange(SatCopPos, RcvrRefPosXyz)
                # Get Sigma UERE from Conf
                UERE = getUERE(Conf, SatLabel)
            else:
                SatCorrInfo["Flag"] = 0
                CorrCode = 0
                CorrPhase = 0
                GeomRange = 0
                UERE = 0
            
            
            # Compute the first Residual removing the geometrical range 
            # They include Receiver clock estimation
            CodeResidual = CorrCode - GeomRange
            PhaseResidual = CorrPhase - GeomRange
            

            # Update CorrInfo dict
            SatCorrInfo["LeoApoX"] = RcvrApoXyz[0]
            SatCorrInfo["LeoApoY"] = RcvrApoXyz[1]
            SatCorrInfo["LeoApoZ"] = RcvrApoXyz[2]
            SatCorrInfo["LeoX"] = RcvrRefPosXyz[0]
            SatCorrInfo["LeoY"] = RcvrRefPosXyz[1]
            SatCorrInfo["LeoZ"] = RcvrRefPosXyz[2]
            SatCorrInfo["FlightTime"] = FlightTime * 1000
            SatCorrInfo["SatComPos"] = SatComPos
            SatCorrInfo["SatApoX"] = SatApo[0]
            SatCorrInfo["SatApoY"] = SatApo[1]
            SatCorrInfo["SatApoZ"] = SatApo[2]
            SatCorrInfo["ApoProj"] = SatApo*RcvrLos
            SatCorrInfo["SatX"] = SatCopPos[0]
            SatCorrInfo["SatY"] = SatCopPos[1]
            SatCorrInfo["SatZ"] = SatCopPos[2]
            SatCorrInfo["SatCodeBia"] = CodeSatBias
            SatCorrInfo["SatPhaseBia"] = PhaseSatBias
            SatCorrInfo["Dtr"] = Dtr
            SatCorrInfo["CorrCode"] = CorrCode
            SatCorrInfo["CorrPhase"] = CorrPhase
            SatCorrInfo["GeomRange"] = GeomRange
            SatCorrInfo["SatClk"] = SatClkBiam
            SatCorrInfo["SigmaUere"] = UERE
            SatCorrInfo["CodeResidual"] = CodeResidual
            SatCorrInfo["PhaseResidual"] = PhaseResidual

            # Sum of CodeResiduals
            if Constel == "G":
                Res_accum_GPS +=  CodeResidual*UERE
                UERE_accum_GPS += UERE
            elif Constel == "E":
                Res_accum_GAL += CodeResidual*UERE
                UERE_accum_GAL += UERE

            # Save data for using it in next Epoch
            SatComPos_1[SatLabel] = [SatComPos[0], SatComPos[1], SatComPos[2]]
            Sod_1[SatLabel] = Sod

        else:
            # Set sat to "don't use"
            SatCorrInfo["Flag"] = 0

        # End of if SatPrepro["Status"] == 1

        # Update info for next epochs
        Sod_1[SatLabel] = Sod
        
        try:
            SatComPos_1[SatLabel] = [SatComPos[0], SatComPos[1], SatComPos[2]]
        except:
            SatComPos_1[SatLabel] = [0, 0, 0]

    # End of for SatLabel, SatPrepro in PreproObsInfo.items()

        try:
            del SatComPos
        except:
            pass   
    

    # Estimate the Receiver Clock first guess as a weighted average of the Residuals
    if UERE_accum_GPS != 0:
        RcvrClk_GPS = Res_accum_GPS / UERE_accum_GPS
    if UERE_accum_GAL != 0:
        RcvrClk_GAL = Res_accum_GAL / UERE_accum_GAL

    for SatLabel, SatPrepro in PreproObsInfo.items():
        
        # Get constellation
        Constel = SatLabel[0]

        if SatCorrInfo["Flag"] == 1:
            if CorrInfo[SatLabel]["Dtr"] == Const.NAN:
                RcvrClk = 0
                CodeResidual = 0
                PhaseResidual = 0
            else:
                if Constel == "G":
                    RcvrClk = RcvrClk_GPS
                elif Constel == "E":
                    RcvrClk = RcvrClk_GAL
                CodeResidual = CorrInfo[SatLabel]["CodeResidual"] - RcvrClk
                PhaseResidual = CorrInfo[SatLabel]["PhaseResidual"] - RcvrClk
            CorrInfo[SatLabel]["RcvrClk"] = RcvrClk
            # Remove the receiver clock from the Residuals
            CorrInfo[SatLabel]["CodeResidual"] = CodeResidual
            CorrInfo[SatLabel]["PhaseResidual"] = PhaseResidual
    
    return CorrInfo, RcvrRefPosXyz, RcvrRefPosLlh, Sod_1, SatComPos_1

#######################################################
# OTHER INTERNAL FUNCTIONS
#######################################################

def computeLeoComPos(Sod, LeoPosInfo):
    '''
    Purpose: compute the Leo COM position for the specified SOD

    Parameters
    ==========
    Sod: int
        Current epoch
    LeoPosInfo: dict
        Dictionary containing the data from SP3: LEO_POS.dat file

    Returns
    =======
    xyzLeoComPos: vect
        Vector containing the three coordinates of the position (meters)
    '''
    # Find the index of the row corresponding to the specified SOD
    rowIdx = LeoPosInfo[LeoPosIdx["SOD"]].index(Sod)

    # Extract the xCoM, yCoM, and zCoM columns from the LeoPosInfo dictionary
    xCoM = LeoPosInfo[LeoPosIdx["xCM"]][rowIdx] * 1000
    yCoM = LeoPosInfo[LeoPosIdx["yCM"]][rowIdx] * 1000
    zCoM = LeoPosInfo[LeoPosIdx["zCM"]][rowIdx] * 1000

    # Save the values as a vector
    xyzLeoComPos = [xCoM, yCoM, zCoM]

    return xyzLeoComPos

# End of computeLeoComPos()

def computeQuat(Sod, LeoQuatInfo):
    '''
    Purpose: get the quaternions from LeoQuatInfo

    Parameters
    ==========
    Sod: int
        Current epoch
    LeoQuatInfo: dict
        Dictionary containing the data from ATT: LEO_QUAT.dat file

    Returns
    =======
    LeoQuat: vect
        Vector containing the quaternions

    '''

    # Find the index of the row corresponding to the specified SOD
    rowIdx = LeoQuatInfo[LeoPosIdx["SOD"]].index(Sod)

    # Extract the xCoM, yCoM, and zCoM columns from the LeoQuatInfo dictionary
    q0 = LeoQuatInfo[LeoQuatIdx["q0"]][rowIdx]
    q1 = LeoQuatInfo[LeoQuatIdx["q1"]][rowIdx]
    q2 = LeoQuatInfo[LeoQuatIdx["q2"]][rowIdx]
    q3 = LeoQuatInfo[LeoQuatIdx["q3"]][rowIdx]

    # Save the values as a vector
    LeoQuat = [q0, q1, q2, q3]

    return LeoQuat

# End of computeLeoComPos()


def computeSatClkBias(Sod, SatClkInfo, SatLabel):
    '''
    Purpose: compute satellite clock bias using linear interpolation

    Parameters
    ==========
    Sod: int
        Current epoch
    SatClkInfo: OrderedDict
        Dictionary containing satellite clock information
    SatLabel: str
        Satellite label (e.g. 'G01', 'E02')

    Returns
    =======
    SatClkBias: float
        Satellite clock bias (seconds)
    '''

    # Filter SatClkInfo based on SatLabel
    filtered_idx = [i for i, label in enumerate(SatClkInfo[SatClkIdx["SatLabel"]]) if label == SatLabel]
    filtered_SatClkSod = np.array([SatClkInfo[SatClkIdx["SOD"]][i] for i in filtered_idx])
    filtered_SatClkBias = np.array([SatClkInfo[SatClkIdx["CLK-BIAS"]][i] for i in filtered_idx])

    # Find the indices of the 2 closest points to the actual Sod
    idx = find_n_closest_points(filtered_SatClkSod, Sod, 2)
    selected_times = filtered_SatClkSod[idx]

    # Extract the position data for the 2 closest points
    selected_SatClkBias = filtered_SatClkBias[idx]

    # Check if there are values greater than the current SOD
    if( ((selected_times[1] - Sod) <= 300) and  ((Sod - selected_times[0]) <= 300)):
        # If Sod is directly in the information, no need to interpolate
        if Sod in filtered_SatClkSod:
            idx = np.argwhere(filtered_SatClkSod == Sod)
            SatClkBias = filtered_SatClkBias[idx][0][0]
        else:
            # Compute the linear interpolation
            SatClkBias = np.interp(Sod, selected_times, selected_SatClkBias)
    else:
        SatClkBias = Const.NAN
        
    return SatClkBias


# End of computeSatClkBias()

def computeRcvrApo(Conf, Year, Doy, Sod, SatLabel, LeoQuatInfo):
    '''
    Purpose: compute the receiver APO in ECEF coordinates

    Parameters
    ==========
    Conf: dict
        Configuration dictionary
    Year: int
        Year
    Doy: int
        Day of year
    Sod: int
        Current epoch
    SatLabel: str
        Satellite label (e.g. 'G01', 'E02')
    LeoQuatInfo: dict
        Dictionary containing the data from ATT: LEO_QUAT.dat file

    Returns
    =======
    RcvrApoXyz: numpy array
        Receiver APO in ECEF coordinates (meters)
    '''

    # Get the configured PCOs in SFR
    if SatLabel[0] == "G":
        Pco = np.array(Conf['LEO_PCO_GPS'])
    elif SatLabel[0] == "E":
        Pco = np.array(Conf['LEO_PCO_GAL'])

    # Define the matrix A_GP1 to transfrom ARF to SRF
    A_GP1 = np.array([[1.0, -0.000707, -0.000236],
                      [-0.000707, -1.0, 0.0],
                      [-0.000236, 0.0, -1.0]])
    # APO in SRF
    ApoXyzSfr = np.array(Conf["LEO_ARP_POS"]) - np.array(Conf["LEO_COM_POS"]) + A_GP1 @ Pco

    # Get the satellite quaternions for the given SatLabel and time
    [q0, q1, q2, q3] = computeQuat(Sod, LeoQuatInfo)

    # Build the rotation matrix from the quaternions
    R = np.array([[1-2*q2**2-2*q3**2, 2*(q1*q2-q0*q3), 2*(q0*q2+q1*q3)],
                  [2*(q1*q2+q0*q3), 1-2*q1**2-2*q3**2, 2*(q2*q3-q0*q1)],
                  [2*(q1*q3-q0*q2), 2*(q0*q1+q2*q3), 1-2*q1**2-2*q2**2]])

    # Rotate the APO from SFR to ECI
    ApoXyzEci = R.dot(ApoXyzSfr)

    # Convert ECI coordinates to ECEF coordinates
    jday = convertYearDoy2JulianDay(Year, Doy, Sod) - 2415020
    fday = Sod / 86400
    gst = np.mod((279.690983 + 0.9856473354*jday + 360*fday + 180), 360)
    gstr = np.deg2rad(gst)
    R_eci_to_ecef = np.array([[np.cos(gstr), np.sin(gstr), 0],
                               [-np.sin(gstr), np.cos(gstr), 0],
                               [0, 0, 1]])
    
    RcvrApoXyz = np.dot(R_eci_to_ecef, ApoXyzEci)

    return RcvrApoXyz

# End of computeRcvrApo()

def computeSatComPos(TransmissionTime, SatPosInfo, SatLabel):
    '''
    Purpose: compute satellite CoM position at transmission time 
    using 10-point Lagrange interpolation

    Parameters
    ==========
    TransmissionTime: float
        Transmission time (seconds)
    SatPosInfo: dict
        Dictionary containing satellite position information

    Returns
    =======
    SatComPos: numpy array
        Satellite CoM position at transmission time (meters)
    '''

    # Filter SatPosInfo based on SatLabel
    filtered_idx = [i for i, label in enumerate(SatPosInfo[SatPosIdx["SatLabel"]]) if label == SatLabel]
    filtered_SatPosSod = np.array([SatPosInfo[SatPosIdx["SOD"]][i] for i in filtered_idx])
    filtered_SatPosX = np.array([SatPosInfo[SatPosIdx["xCM"]][i] for i in filtered_idx])
    filtered_SatPosY = np.array([SatPosInfo[SatPosIdx["yCM"]][i] for i in filtered_idx])
    filtered_SatPosZ = np.array([SatPosInfo[SatPosIdx["zCM"]][i] for i in filtered_idx])


    # Find the indices of the 10 closest points to the transmission time
    idx = find_n_closest_points(filtered_SatPosSod, TransmissionTime, 10)
    selected_times = filtered_SatPosSod[idx]

    # Extract the position data for the 10 closest points (meters)
    selected_positionsX = filtered_SatPosX[idx]* 1000
    selected_positionsY = filtered_SatPosY[idx]* 1000
    selected_positionsZ = filtered_SatPosZ[idx] * 1000

    satPosXYZ = np.array([
        selected_positionsX,
        selected_positionsY,
        selected_positionsZ
    ])

    Lx, Ly, Lz = lagrange_int(TransmissionTime, selected_times, satPosXYZ)
    SatComPos = np.array([Lx, Ly, Lz])

    return SatComPos

# End of computeSatComPos()

def lagrange_int(x, x_points, y_points):

    '''
    Purpose: build Lagrange interpolation at a given time

    Parameters
    ==========
    x: float
        Point where the interpolation must be done
    x_points: list
        List of x-points for the interpolation
    y_points: list
        List of y-points for the interpolation

    Returns
    =======
    Lx, Ly, Lz: float
        Results from the interpolation at each coordinate
    '''
    assert len(x_points) == len(y_points[0]), "values must have the same length."

    # Build Lagrange basis 
    Lx, Ly, Lz = 0, 0, 0
    
    for i in range(len(x_points)):
        # Calculate Lagrange basis polynomial L_i(x)
        L_i = 1
        for j in range(len(x_points)):
            if i != j:
                L_i *= (x - x_points[j]) / (x_points[i] - x_points[j])
        
        # Add the term y_i * L_i(x) to the interpolating polynomial
        Lx += y_points[0][i] * L_i
        Ly += y_points[1][i] * L_i
        Lz += y_points[2][i] * L_i

    return Lx, Ly, Lz

# End of lagrange_int()

def applySagnac(SatComPos, FlightTime):
    '''
    Purpose: Apply Sagnac correction to the satellite position

    Parameters
    ==========
    SatComPos: numpy array
        Satellite CoM position at transmission time (meters)
    FlightTime: float
        Time of flight (s)

    Returns
    =======
    SatComPos: numpy array
        Satellite CoM position at transmission time corrected by Sagnac (meters)
    '''
    theta = Const.OMEGA_EARTH * FlightTime
    R3= np.array([[np.cos(theta), np.sin(theta), 0],
                  [-np.sin(theta), np.cos(theta), 0],
                  [0, 0, 1]])
    
    SatComPos = np.dot(R3, SatComPos)
    
    return SatComPos
# End of applySagnac()

def computeSatApo(SatLabel, SatComPos, SunPos, SatApoInfo, gamma):
    '''
    Purpose: Compute the APO in ECEF from ANTEX APOs in SFR

    Parameters
    ==========
    SatComPos: numpy array
       Satellite Center of Mass (CoM) position in ECEF corrected by Sagnac
    RcvrPos: numpy array
       Receiver (observer) position in ECEF
    SatApoInfo: numpy array
       Antenna Phase Offset (APO) in satellite-body reference frame

    Returns
    ==========
    SatApo: np array
        Satellite antenna phase offset for each coordinate (meters)
    '''

    # Calculate the unit vectors
    k_hat = -SatComPos / np.linalg.norm(SatComPos)
    e_hat = (SunPos-SatComPos) / np.linalg.norm(SunPos-SatComPos)
    j_hat = crossProd(k_hat, e_hat)
    j_hat = j_hat / np.linalg.norm(j_hat)
    i_hat = crossProd(j_hat, k_hat)

    # Construct rotation matrix from i, j, k unit vectors
    Rot = np.column_stack([i_hat, j_hat, k_hat])

    # Rotate the APO (meters)
    idx = SatApoInfo[SatApoIdx["SatLabel"]].index(SatLabel)

    xIF = (SatApoInfo[SatApoIdx["x_f2"]][idx] - gamma * SatApoInfo[SatApoIdx["x_f1"]][idx]) / (1-gamma)
    yIF = (SatApoInfo[SatApoIdx["y_f2"]][idx] - gamma * SatApoInfo[SatApoIdx["y_f1"]][idx]) / (1-gamma)
    zIF = (SatApoInfo[SatApoIdx["z_f2"]][idx] - gamma * SatApoInfo[SatApoIdx["z_f1"]][idx]) / (1-gamma)
    
    SatApo = [xIF/1000, yIF/1000, zIF/1000]
    SatApo = np.dot(Rot, SatApo)

    return SatApo
# End of computeSatApo

def getSatBias(gamma, SatLabel, SatBiaInfo):
    '''
    Purpose: get Satellite Biases from SatBiaInfo 
    split into Code and Phase biases

    Parameters
    ==========
    gamma: float
        GammaF1F2 for GPS or GAL
    SatBiaInfo: dict
        Dictionary containing information from Satellite Signal Biases

    Returns
    =======
    SatelliteBiases: np array
        SatelliteBiases = [CodeBiases PhaseBiases] (meters)
    '''

    # Filter SatBiaInfo based on SatLabel
    filtered_idx = [i for i, label in enumerate(SatBiaInfo[SatBiaIdx["SatLabel"]]) if label == SatLabel]
    filtered_CLKF1C = np.array([SatBiaInfo[SatBiaIdx["CLK_f1_C"]][i] for i in filtered_idx])
    filtered_CLKF2C = np.array([SatBiaInfo[SatBiaIdx["CLK_f2_C"]][i] for i in filtered_idx])
    filtered_OBSF1C = np.array([SatBiaInfo[SatBiaIdx["OBS_f1_C"]][i] for i in filtered_idx])
    filtered_OBSF2C = np.array([SatBiaInfo[SatBiaIdx["OBS_f2_C"]][i] for i in filtered_idx])
    filtered_CLKF1P = np.array([SatBiaInfo[SatBiaIdx["CLK_f1_P"]][i] for i in filtered_idx])
    filtered_CLKF2P = np.array([SatBiaInfo[SatBiaIdx["CLK_f2_P"]][i] for i in filtered_idx])
    filtered_OBSF1P = np.array([SatBiaInfo[SatBiaIdx["OBS_f1_P"]][i] for i in filtered_idx])
    filtered_OBSF2P = np.array([SatBiaInfo[SatBiaIdx["OBS_f2_P"]][i] for i in filtered_idx])

    # Compute code biases (meters)
    BIA_CLK_C = ((filtered_CLKF2C-gamma*filtered_CLKF1C) / (1-gamma)) * Const.SPEED_OF_LIGHT / 10**9 
    BIA_OBS_C = ((filtered_OBSF2C-gamma*filtered_OBSF1C) / (1-gamma)) * Const.SPEED_OF_LIGHT / 10**9 
    
    CodeSatBias = BIA_CLK_C - BIA_OBS_C

    # Compute phase biases (meters)
    BIA_CLK_P = ((filtered_CLKF2P-gamma*filtered_CLKF1P) / (1-gamma)) * Const.SPEED_OF_LIGHT / 10**9 
    BIA_OBS_P = ((filtered_OBSF2P-gamma*filtered_OBSF1P) / (1-gamma)) * Const.SPEED_OF_LIGHT / 10**9 
    
    PhaseSatBias =  BIA_CLK_P - BIA_OBS_P

    SatelliteBiases = [CodeSatBias, PhaseSatBias]

    return SatelliteBiases
#End of getSatBias()

def computeDtr(SatComPos_1, SatComPos, Sod, Sod_1, SatLabel):
    '''
    Purpose: compute the Relativistic Correction (DTR)

    Parameters
    ==========
    SatComPos_1: numpy array
       Satellite Center of Mass (CoM) position in ECEF corrected by Sagnac from previous Epoch
    SatComPos: numpy array
       Satellite Center of Mass (CoM) position in ECEF corrected by Sagnac
    Sod: int
        Current epoch
    Sod_1: int
        Previous epoch
            

    Returns
    =======
    DTR: float
      Relativistic correction (meters)
    '''
    Sod_1 = Sod_1[SatLabel]
    SatComPos_1 = np.array(SatComPos_1[SatLabel])

    # If SatComPos has been computed previously Dtr is computed
    if SatComPos_1[0] == 0:
        Dtr = Const.NAN
    else:
        # Estimate satellite speed
        vsat = (SatComPos - SatComPos_1) / (Sod - Sod_1)

        # Compute Dtr
        Dtr = -2*np.dot(SatComPos,vsat) / Const.SPEED_OF_LIGHT

    return Dtr
# End of computeDtr()

def getUERE(Conf, SatLabel):
    '''
    Purpose: get UERE from configuration for each satellite

    Parameters
    ==========
    Conf: dict
        Configuration dictionary
    SatLabel: str
        Satellite label (e.g. 'G01', 'E02')

    Returns
    =======
    UERE: float
        Sigma UERE depending on the constellation 
    '''
    if SatLabel.startswith("G"):
        UERE = Conf["GPS_UERE"]
    else:
        UERE = Conf["GAL_UERE"]
    return UERE
# End of getUERE()

def ComputeGeomRange(SatCopPos, RcvrRefPosXyz):
    '''
    Purpose: compute the Geometrical Range

    Parameters
    ==========
    SatComPos: numpy array
       Satellite Position at Center of Phases (m)
    RcvrRefPosXyz: numpy array
        Leo position at Center of Phases (m)

    Returns
    =======
    GeomRange: float
      Geometrical Range (m)
    '''

    GeomRange = np.linalg.norm(SatCopPos-RcvrRefPosXyz)
    return GeomRange

### Find closest points ###

def find_n_closest_points(points, reference, n):
    '''
    Purpose: find the points for the interpolation

    Parameters
    ==========
    points: list
        The list of points to check
        
    reference: float
        The reference value to compute distances
    
    n: int
        The number of closest points to find

    Returns
    =======
    points: list
        The closest points in the original list
    
    idx: list
        The sorted indices (from smallest to largest) of the closest points
    '''
    # Calculate the absolute differences between each point and the reference value
    differences = [(abs(point - reference), index) for index, point in enumerate(points)]
    
    # Sort based on the absolute differences
    differences.sort()
    
    # Get the n closest points and their indices
    idx = [diff[1] for diff in differences[:n]]

    # Sort the results by the indices
    sorted_idx = sorted(idx)
    
    return sorted_idx


### Filter by constellation ###

def filterByConstel(data, constel):
    '''
    Purpose: filter the data by constellation

    Parameters
    ==========
    data: dict
        Dictionary containing the data from the data file
    constel: str
        Constellation to filter by (GPS, GAL, or GPSGAL)

    Returns
    =======
    filtered data: dict
        OrderedDict with info from file filtered by constellation
    '''

    filteredData = {}

    # Get the index of the PRN column
    prnIndex = len(data) - 1

    # Filter the rows based on the PRN column and the specified constellation
    for i in range(len(data[prnIndex])):
        if constel == 'GPS' and data[prnIndex][i].startswith('G'):
            for j in range(len(data)):
                if j not in filteredData:
                    filteredData[j] = []
                filteredData[j].append(data[j][i])
        elif constel == 'GAL' and data[prnIndex][i].startswith('E'):
            for j in range(len(data)):
                if j not in filteredData:
                    filteredData[j] = []
                filteredData[j].append(data[j][i])
        elif constel == 'GPSGAL':
            for j in range(len(data)):
                if j not in filteredData:
                    filteredData[j] = []
                filteredData[j].append(data[j][i])

    return filteredData

# End of filterByConstel


def filterPreproObsByConstel(data, constel):
    '''
    Purpose: filter the data from PreproObs by constellation

    Parameters
    ==========
    data: OrderedDict
        Dictionary containing the data from the data file
    constel: str
        Constellation to filter by (GPS, GAL, or GPSGAL)

    Returns
    =======
    filtered data: dict
        OrderedDict with info from file filtered by constellation
    '''

    filteredData = OrderedDict()

    # Filter the rows based on the PRN column and the specified constellation
    for prn, values in data.items():
        if constel == 'GPS' and prn.startswith('G'):
            filteredData[prn] = values
        elif constel == 'GAL' and prn.startswith('E'):
            filteredData[prn] = values
        elif constel == 'GPSGAL':
            filteredData[prn] = values

    return filteredData

def filterObsInfoByConstel(data, constel): # **********************************
    '''
    Purpose: filter the data from ObsInfo by constellation

    Parameters
    ==========
    data: list of lists
            list containing the data from the data file
    constel: str
            constellation to filter by (GPS, GAL, or GPSGAL)

    Returns
    =======
    list of lists with info from file filtered by constellation
    '''

    filteredData = []

    # Filter the rows based on the PRN column and the specified constellation
    for row in data:
        prn = row[2]
        if constel == 'GPS' and prn.startswith('G'):
            filteredData.append(row)
        elif constel == 'GAL' and prn.startswith('E'):
            filteredData.append(row)
        elif constel == 'GPSGAL':
            filteredData.append(row)

    return filteredData