#!/usr/bin/env python

########################################################################
# SENTUS-V3/SRC/Pvts.py:
# This is the PVT Solution Module of SENTUS tool
#
#  Project:        SENTUS
#  File:           Pvts.py
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
from COMMON.Coordinates import ecef2enu
import numpy as np
from COMMON.Wlsq import buildGmatrixEl, buildGmatrix
from COMMON.Wlsq import buildWmatrix
from COMMON.Wlsq import computeDOPs
from COMMON.Wlsq import buildSmatrix
from COMMON.Wlsq import runWlsqIteration
from Perf import updatePerfEpoch

def computeWlsqSolution(CorrInfo, Conf, Sod, RcvrRefPosLlh, PerfInfo):
    '''
    Purpose:
        This function computes the position and clock bias estimation for a GNSS receiver
        using a Weighted Least Squares (WLS) approach based on satellite observations.
        It constructs the observation and weighting matrices from valid satellite corrections
        and iteratively solves for the receiver's position in LLH (Longitude, Latitude, Altitude)
        coordinates, as well as the clock bias relative to GPS time.

    Parameters
    ==========
    CorrInfo : dict
        A dictionary containing satellite correction information, where each key is a satellite label
        and each value is another dictionary with the satellite's position, clock corrections, 
        measurement flags, and other relevant data.
    
    Conf : dict
        A configuration dictionary containing settings for the computation.

    Sod : int
        The second of the day, representing the time of the measurements in seconds.

    RcvrRefPosLlh : list
        The reference position of the receiver in latitude, longitude, and height (LLH) coordinates.

    PerfInfo : dict
        A dictionary to track performance metrics, such as sample counts and solution statuses,
        updated during the computation process.

    Returns
    ==========
    PosInfo : dict
        A dictionary containing the computed position and clock bias estimations, as well as status
        indicators.

    PerfInfo : dict
        An updated performance dictionary reflecting the current state of the receiver's performance.
    '''


    # Initialize output
    PosInfo = OrderedDict({})
    PosInfo = {
        "SOD": Sod,            # Second of day
        "LONEST": 0.0,         # Longitude estimation
        "LATEST": 0.0,         # Latitude estimation
        "ALTEST": 0.0,         # Altitude estimation
        "CLKEST": 0.0,         # Clock bias estimation
        "GGTO": 0.0,           # GNSS-GPS Time Offset
        "SOL": 0,              # Solution status
        "NSVVIS": 0,           # Number of satellites visible
        "NSV": 0,              # Number of satellites used in the solution
        "HPE": 0.0,            # Horizontal Position Error
        "VPE": 0.0,            # Vertical Position Error
        "EPE": 0.0,            # Total Position Error
        "NPE": 0.0,            # North Position Error
        "UPE": 0.0,            # Up Position Error
        "HDOP": 0.0,           # Horizontal Dilution of Precision
        "VDOP": 0.0,           # Vertical Dilution of Precision
        "PDOP": 0.0            # Position Dilution of Precision
    } # End of PosInfo

    # Build G: Observation and W: Weighting Matrices
    # -------------------------------------------------------------------------------------------

    # Filter valid satellites
    Val_CorrInfo = {SatLabel: info for SatLabel, info in CorrInfo.items() if info["Flag"] == 1}

    # Check number of visible sats and valid sats    
    PosInfo["NSVVIS"] = len(CorrInfo)
    NSats = len(Val_CorrInfo)
    PosInfo["NSV"] = NSats

    # Initialize matrices
    GMat = np.zeros((NSats, 5)) # *********** Multi-const 4 --> 5
    Wmat = np.zeros((NSats, NSats))
    RcvrPosClkECEF = np.zeros(5)
    RcvrPosENU = np.zeros(3)
    RcvrPosClkENU = np.zeros(5)
    idx=0

    for SatLabel, SatCorr in Val_CorrInfo.items():
        if SatLabel.startswith("G"): 
            RcvrPosClkECEF[0] = SatCorr["LeoX"]
            RcvrPosClkECEF[1] = SatCorr["LeoY"]
            RcvrPosClkECEF[2] = SatCorr["LeoZ"]
            RcvrPosClkECEF[3] = SatCorr["RcvrClk"]
            break 
    for SatLabel, SatCorr in Val_CorrInfo.items():
        if SatLabel.startswith("E"): 
            RcvrPosClkECEF[4] = SatCorr["RcvrClk"]
            break 

    ## ECEF --> ENU with LEO as reference
    RcvrPosENU = ecef2enu(RcvrRefPosLlh[0], RcvrRefPosLlh[1], RcvrPosClkECEF[:3], RcvrPosClkECEF[:3])
    RcvrPosClkENU = np.array([RcvrPosENU[0], RcvrPosENU[1], RcvrPosENU[2], RcvrPosClkECEF[3], RcvrPosClkECEF[4]])

    # Loop over all valid Satellite Corrected Measurements
    for SatLabel, SatCorr in Val_CorrInfo.items():
        # Get Constel
        Constel = SatLabel[0]

        # Get SatPos
        # Satellite position in ECEF
        SatCoord= np.array([SatCorr["SatX"], SatCorr["SatY"], SatCorr["SatZ"]])

        # ECEF --> ENU
        SatCoordENU = ecef2enu(RcvrRefPosLlh[0], RcvrRefPosLlh[1], SatCoord, RcvrPosClkECEF[:3])

        # Build Geometry matrix in line with SBAS Standard in ENU Coordinates
        GMat = buildGmatrix(GMat, SatCoordENU, idx, Constel, RcvrPosClkENU) # In case of computing G matrix from ENU Coordinates
        # GMat = buildGmatrixEl(GMat, SatCorr, idx, Constel) # In case of computing G matrix from elevation and azimuth
        
        # Build Weighting matrix in line with SBAS Standard
        Wmat = buildWmatrix(Wmat, SatCorr["SigmaUere"], idx)

        idx += 1  # Increment the index


    # Perform PVT Solution using all SBS with Flag == 1
    # -------------------------------------------------------------------------------------------

    # Compute solution if the minimum required satellites are available
    if PosInfo["NSV"] >= Const.MIN_NUM_SATS_PVT:
        # Compute DOPs
        PDOP, HDOP, VDOP = computeDOPs(GMat, RcvrRefPosLlh)

        PosInfo["PDOP"] = PDOP
        PosInfo["HDOP"] = HDOP
        PosInfo["VDOP"] = VDOP

        # Check if the PDOP is below the configured threshold
        if PDOP <= Conf["PDOP_MAX"]:

            # Compute the Projection Matrix in ENU Coordinates
            Smat = buildSmatrix(GMat, Wmat)

            # Get the estimated PVT solution: Position and Clock applying WLSE filter
            # ------------------------------------------------------------------------------------
            # Initial Position guess is the Receiver Reference Position in LLH
            RcvrPosClk = np.array([RcvrRefPosLlh[0], RcvrRefPosLlh[1], \
                                   RcvrRefPosLlh[2], 0, 0])

            NumIter = 0
            RcvrPosClkDelta_accum = 0
            RcvrPosClkDelta = np.array([np.inf, 0, 0, 0, 0])

            while NumIter <= Conf["MAX_LSQ_ITER"] and np.linalg.norm(RcvrPosClkDelta) > Const.LSQ_DELTA_EPS:

                NumIter, RcvrPosClkDelta, RcvrPosClk = runWlsqIteration(Val_CorrInfo, Smat, RcvrPosClk, NumIter)

                # Estimate Errors  in ENU Coord
                RcvrPosClkDelta_accum += RcvrPosClkDelta

            GGTO_meters = RcvrPosClk[4] - RcvrPosClk[3]

            PosInfo["LONEST"] = RcvrPosClk[0]
            PosInfo["LATEST"] = RcvrPosClk[1]
            PosInfo["ALTEST"] = RcvrPosClk[2]
            PosInfo["CLKEST"] = RcvrPosClk[3]
            PosInfo["GGTO"] = GGTO_meters/(Const.SPEED_OF_LIGHT)*10**9
            PosInfo["SOL"] = 1
            PosInfo["NPE"] = RcvrPosClkDelta_accum[1]
            PosInfo["EPE"] = RcvrPosClkDelta_accum[0]
            PosInfo["UPE"] = RcvrPosClkDelta_accum[2]
            PosInfo["HPE"] = np.sqrt(PosInfo["EPE"] **2+PosInfo["NPE"] **2)
            PosInfo["VPE"] = np.abs(PosInfo["UPE"])

            # Update Performance Receiver Information per epoch
            # ----------------------------------------------------------
            PerfInfo["SAMPLES"] += 1
            PerfInfo = updatePerfEpoch(PosInfo, PerfInfo)
        else:
            PosInfo["SOL"] = 0
            PerfInfo["SAMNOSOL"] += 1
    else:
        PosInfo["SOL"] = 0
        PerfInfo["SAMNOSOL"] += 1

    return PosInfo, PerfInfo


