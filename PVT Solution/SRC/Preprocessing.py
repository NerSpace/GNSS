#!/usr/bin/env python

########################################################################
# Preprocessing.py:
# This is the Preprocessing Module of SENTUS tool
#
#  Project:        SENTUS
#  File:           Preprocessing.py
#
#   Author: GNSS Academy
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
from InputOutput import ObsIdxC, ObsIdxP, REJECTION_CAUSE
from InputOutput import FLAG, VALUE, TH, CSNEPOCHS, CSNPOINTS, CSPDEGREE
import numpy as np

# Preprocessing internal functions
#-----------------------------------------------------------------------


def runPreprocessing(Conf, ObsInfo, PrevPreproObsInfo):
    '''
    Purpose: preprocess GNSS raw measurements from OBS file
             and generate PREPRO OBS file with the cleaned,
             smoothed measurements

             More in detail, this function handles:
             
             * Measurements cleaning and validation and exclusion due to different 
             criteria as follows:
                - Minimum Masking angle
                - Maximum Number of channels
                - Minimum Carrier-To-Noise Ratio (CN0)
                - Pseudo-Range Output of Range 
                - Maximum Pseudo-Range Step
                - Maximum Pseudo-Range Rate
                - Maximum Carrier Phase Increase
                - Maximum Carrier Phase Increase Rate
                - Data Gaps checks and handling 
                - Cycle Slips detection

            * Filtering/Smoothing of Code-Phase Measurements with a Hatch filter 

    Parameters
    ==========
    Conf: dict
            Configuration dictionary
    ObsInfo: list
            OBS info for current epoch
    PrevPreproObsInfo: dict
            Preprocessed observations for previous epoch per sat
            PrevPreproObsInfo["G01"]["C1"]

    Returns
    =======
    PreproObsInfo: dict
            Preprocessed observations for current epoch per sat
            PreproObsInfo["G01"]["C1"]
    '''

    # Get Observations
    CodesObs = ObsInfo[0]
    PhaseObs = ObsInfo[1]

    # if check Cycle Slips activated
    if (Conf["CYCLE_SLIPS"][FLAG] == 1):

        # Loop over Phase measurements
        for iObs, SatPhaseObs in enumerate(PhaseObs):

            # Get SOD
            Sod = float(SatPhaseObs[ObsIdxP["SOD"]])

            # Get satellite label
            SatLabel = SatPhaseObs[ObsIdxP["PRN"]]

            # Compute Geometry-Free in cycles
            GF_Lcy = float(SatPhaseObs[ObsIdxP["L1"]]) - float(SatPhaseObs[ObsIdxP["L2"]])

            # Check Data Gaps
            DeltaT = Sod - PrevPreproObsInfo[SatLabel]["GF_Epoch_Prev"][-1]
            if DeltaT > Conf["MAX_DATA_GAP"][TH]:
                # Reinitialize CS detection
                PrevPreproObsInfo[SatLabel]["GF_L_Prev"] = [0.0] * int(Conf["CYCLE_SLIPS"][CSNPOINTS])
                PrevPreproObsInfo[SatLabel]["GF_Epoch_Prev"] = [0.0] * int(Conf["CYCLE_SLIPS"][CSNPOINTS])
                PrevPreproObsInfo[SatLabel]["CycleSlipBuffIdx"] = 0
                PrevPreproObsInfo[SatLabel]["CycleSlipFlags"] = [0.0] * int(Conf["CYCLE_SLIPS"][CSNEPOCHS])
                PrevPreproObsInfo[SatLabel]["CycleSlipFlagIdx"] = 0

            # Cycle slips detection
            # fit a polynomial using previous GF measurements to compare the predicted value
            # with the observed one
            # --------------------------------------------------------------------------------------------------------------------
            # Get N
            N = PrevPreproObsInfo[SatLabel]["CycleSlipBuffIdx"]

            # If the buffer is full, we can detect the cycle slip with a polynom
            if N == (Conf["CYCLE_SLIPS"][CSNPOINTS]):
                # Adjust polynom to the samples in the buffer
                Polynom = np.polynomial.polynomial.polyfit(PrevPreproObsInfo[SatLabel]["GF_Epoch_Prev"],
                PrevPreproObsInfo[SatLabel]["GF_L_Prev"],
                int(Conf["CYCLE_SLIPS"][CSPDEGREE]))

                # Predict value evaluating the polynom
                TargetPred = np.polynomial.polynomial.polyval(Sod,
                Polynom)

                # Compute Residual
                Residual = abs(GF_Lcy - TargetPred)

                # Compute CS flag
                CsFlag = Residual > Conf["CYCLE_SLIPS"][TH]

                # Update CS flag buffer
                PrevPreproObsInfo[SatLabel]["CycleSlipFlagIdx"] = \
                    (PrevPreproObsInfo[SatLabel]["CycleSlipFlagIdx"] + 1) % \
                        int(Conf["CYCLE_SLIPS"][CSNEPOCHS])
                PrevPreproObsInfo[SatLabel]["CycleSlipFlags"][PrevPreproObsInfo[SatLabel]["CycleSlipFlagIdx"]] = \
                    CsFlag

                # Check if threshold was exceeded CSNEPOCHS times
                if (np.sum(PrevPreproObsInfo[SatLabel]["CycleSlipFlags"]) == \
                    int(Conf["CYCLE_SLIPS"][CSNEPOCHS])):
                    # Cycle Slip detected
                    PrevPreproObsInfo[SatLabel]["CycleSlipDetectFlag"] = 1

                    # Reinitialize some variables
                    PrevPreproObsInfo[SatLabel]["PrevC1"] = Const.NAN                                     # Previous Code
                    PrevPreproObsInfo[SatLabel]["PrevC2"] = Const.NAN                                     # Previous Code
                    PrevPreproObsInfo[SatLabel]["PrevL1"] = Const.NAN                                     # Previous Phase
                    PrevPreproObsInfo[SatLabel]["PrevL2"] = Const.NAN                                     # Previous Phase
                    PrevPreproObsInfo[SatLabel]["PrevRangeRateL1"] = Const.NAN                            # Previous Code Rate
                    PrevPreproObsInfo[SatLabel]["PrevRangeRateL2"] = Const.NAN                            # Previous Code Rate
                    PrevPreproObsInfo[SatLabel]["PrevPhaseRateL1"] = Const.NAN                            # Previous Phase Rate
                    PrevPreproObsInfo[SatLabel]["PrevPhaseRateL2"] = Const.NAN                            # Previous Phase Rate

                    # Reinitialize CS detection
                    PrevPreproObsInfo[SatLabel]["GF_L_Prev"] = [0.0] * int(Conf["CYCLE_SLIPS"][CSNPOINTS])
                    PrevPreproObsInfo[SatLabel]["GF_Epoch_Prev"] = [0.0] * int(Conf["CYCLE_SLIPS"][CSNPOINTS])
                    PrevPreproObsInfo[SatLabel]["CycleSlipBuffIdx"] = 0
                    PrevPreproObsInfo[SatLabel]["CycleSlipFlags"] = [0.0] * int(Conf["CYCLE_SLIPS"][CSNEPOCHS])
                    PrevPreproObsInfo[SatLabel]["CycleSlipFlagIdx"] = 0

                    # Raise flag to reset Hatch filter
                    PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1

                # If threshold was exceeded less than CSNEPOCHS times,
                # don't update the buffer and set the measurement to invalid
                elif CsFlag == 1:
                    # PrevPreproObsInfo[SatLabel]["CycleSlipDetectFlag"] = -1
                    pass

                # If threshold was not exceeded
                else:
                    # Leave space for the new sample
                    PrevPreproObsInfo[SatLabel]["GF_L_Prev"][:-1] = PrevPreproObsInfo[SatLabel]["GF_L_Prev"][1:]
                    PrevPreproObsInfo[SatLabel]["GF_Epoch_Prev"][:-1] = PrevPreproObsInfo[SatLabel]["GF_Epoch_Prev"][1:]

                    # Store new sample
                    PrevPreproObsInfo[SatLabel]["GF_L_Prev"][-1] = GF_Lcy
                    PrevPreproObsInfo[SatLabel]["GF_Epoch_Prev"][-1] = Sod

                # End of if (np.sum(PrevPreproObsInfo[SatLabel]["CycleSlipFlags"])

            # Buffer is not full, need to add new GF observable
            else:
                PrevPreproObsInfo[SatLabel]["GF_L_Prev"][N] = GF_Lcy
                PrevPreproObsInfo[SatLabel]["GF_Epoch_Prev"][N] = Sod
                PrevPreproObsInfo[SatLabel]["CycleSlipBuffIdx"] += 1

                # PrevPreproObsInfo[SatLabel]["CycleSlipDetectFlag"] = -1

            # End of if N == (Conf["CYCLE_SLIPS"][CSNPOINTS])::

        # end of for iObs, SatPhaseObs in enumerate(PhaseObs):

    # End of if (Conf["CYCLE_SLIPS"][FLAG] == 1)

    # Initialize output
    PreproObsInfo = OrderedDict({})

    # Loop over Code measurements
    for iObs, SatCodesObs in enumerate(CodesObs):

        # Get satellite label
        SatLabel = SatCodesObs[ObsIdxC["PRN"]]

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

        # Get Phases
        SatPhaseObs = PhaseObs[iObs]

        while SatPhaseObs[ObsIdxP["PRN"]] != SatLabel:
            PhaseObs.pop(iObs)
            SatPhaseObs = PhaseObs[iObs]

        assert(SatLabel == SatPhaseObs[ObsIdxP["PRN"]])

        # Initialize output info
        SatPreproObsInfo = {
            "Sod": 0.0,                   # Second of day
            
            "Elevation": 0.0,             # Elevation
            "Azimuth": 0.0,               # Azimuth
            
            "C1": 0.0,                    # GPS L1C/A pseudorange
            "C2": 0.0,                    # GPS L1P pseudorange
            "L1": 0.0,                    # GPS L1 carrier phase (in cycles)
            "L1Meters": 0.0,              # GPS L1 carrier phase (in m)
            "S1": 0.0,                    # GPS L1C/A C/No
            "L2": 0.0,                    # GPS L2 carrier phase 
            "L2Meters": 0.0,              # GPS L2 carrier phase  (in m)
            "S2": 0.0,                    # GPS L2 C/No
            
            "GeomFree_P": Const.NAN,      # Geometry-free of Phases
            "IF_C": Const.NAN,            # Iono-free of Codes
            "IF_P": Const.NAN,            # Iono-free of Phases (pre-aligned)
            "SmoothIF": Const.NAN,        # Smoothed Iono-free of Codes 
            
            "Valid": 1,                   # Measurement Status
            "RejectionCause": 0,          # Cause of rejection flag
            "Status": 0,                  # Smoothing status
            
            "RangeRateL1": Const.NAN,     # L1 Code Rate
            "RangeRateStepL1": Const.NAN, # L1 Code Rate Step
            "PhaseRateL1": Const.NAN,     # L1 Phase Rate
            "PhaseRateStepL1": Const.NAN, # L1 Phase Rate Step
            
            "RangeRateL2": Const.NAN,     # L2 Code Rate
            "RangeRateStepL2": Const.NAN, # L2 Code Rate Step
            "PhaseRateL2": Const.NAN,     # L2 Phase Rate
            "PhaseRateStepL2": Const.NAN, # L2 Phase Rate Step

        } # End of SatPreproObsInfo

        # Prepare outputs
        # Get SoD
        SatPreproObsInfo["Sod"] = float(SatCodesObs[ObsIdxC["SOD"]])
        # Get Elevation
        SatPreproObsInfo["Elevation"] = float(SatCodesObs[ObsIdxC["ELEV"]])
        # Get Azimuth
        SatPreproObsInfo["Azimuth"] = float(SatCodesObs[ObsIdxC["AZIM"]])
        # Get C1
        SatPreproObsInfo["C1"] = float(SatCodesObs[ObsIdxC["C1"]])
        # Get C1
        SatPreproObsInfo["C2"] = float(SatCodesObs[ObsIdxC["C2"]])
        # Get L1 in cycles and in m
        SatPreproObsInfo["L1"] = float(SatPhaseObs[ObsIdxP["L1"]])
        SatPreproObsInfo["L1Meters"] = float(SatPhaseObs[ObsIdxP["L1"]]) * Wave["F1"]
        # Get S1
        SatPreproObsInfo["S1"] = float(SatCodesObs[ObsIdxC["S1"]])
        # Get L2 in cycles and in m
        SatPreproObsInfo["L2"] = float(SatPhaseObs[ObsIdxP["L2"]])
        SatPreproObsInfo["L2Meters"] = float(SatPhaseObs[ObsIdxP["L2"]]) * Wave["F2"]
        # Get S2
        SatPreproObsInfo["S2"] = float(SatCodesObs[ObsIdxC["S2"]])
        
        # Prepare output for the satellite
        PreproObsInfo[SatLabel] = SatPreproObsInfo

    # Loop over satellites
    for SatLabel, PreproObs in PreproObsInfo.items():

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

            # Gamma
            GammaF1F2 = Const.GAL_GAMMA_E1E5A

        # Get epoch
        Epoch = PreproObs["Sod"]

        # Check data gaps
        # ----------------------------------------------------------
        # Compute gap between previous and current observation
        DeltaT = Epoch - PrevPreproObsInfo[SatLabel]["PrevEpoch"]

        # If there is a gap
        if DeltaT > Conf["MAX_DATA_GAP"][TH]:
            # If check gaps activated
            if Conf["MAX_DATA_GAP"][FLAG] == 1:
                # Reject the measurement and indicate the rejection cause
                # (if it is not a non-visibility period)
                if(DeltaT < 1000):
                    PreproObs["RejectionCause"] = REJECTION_CAUSE["DATA_GAP"]

            # Reinitialize some variables
            PrevPreproObsInfo[SatLabel]["PrevC1"] = Const.NAN                                     # Previous Code
            PrevPreproObsInfo[SatLabel]["PrevC2"] = Const.NAN                                     # Previous Code
            PrevPreproObsInfo[SatLabel]["PrevL1"] = Const.NAN                                     # Previous Phase
            PrevPreproObsInfo[SatLabel]["PrevL2"] = Const.NAN                                     # Previous Phase
            PrevPreproObsInfo[SatLabel]["PrevRangeRateL1"] = Const.NAN                            # Previous Code Rate
            PrevPreproObsInfo[SatLabel]["PrevRangeRateL2"] = Const.NAN                            # Previous Code Rate
            PrevPreproObsInfo[SatLabel]["PrevPhaseRateL1"] = Const.NAN                            # Previous Phase Rate
            PrevPreproObsInfo[SatLabel]["PrevPhaseRateL2"] = Const.NAN                            # Previous Phase Rate

            # Reinitialize CS detection
            PrevPreproObsInfo[SatLabel]["GF_L_Prev"] = [0.0] * int(Conf["CYCLE_SLIPS"][CSNPOINTS])
            PrevPreproObsInfo[SatLabel]["GF_Epoch_Prev"] = [0.0] * int(Conf["CYCLE_SLIPS"][CSNPOINTS])
            PrevPreproObsInfo[SatLabel]["CycleSlipBuffIdx"] = 0
            PrevPreproObsInfo[SatLabel]["CycleSlipFlags"] = [0.0] * int(Conf["CYCLE_SLIPS"][CSNEPOCHS])
            PrevPreproObsInfo[SatLabel]["CycleSlipFlagIdx"] = 0

            # Raise flag to reset Hatch filter
            PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1

        # End of if DeltaT > Conf["MAX_DATA_GAP"][TH]:

        # If satellite shall be rejected due to mask angle
        # ----------------------------------------------------------
        if PreproObs["Elevation"] < Conf["RCVR_MASK"]:
            # Indicate the rejection cause
            PreproObs["RejectionCause"] = REJECTION_CAUSE["MASKANGLE"]
            PreproObs["Valid"] = 0

        # If satellite shall be rejected due to C/N0 (only if activated in conf)
        # --------------------------------------------------------------------------------------------------------------------
        if (Conf["MIN_SNR"][FLAG] == 1):
            # 42000 E05 <20
            if(PreproObs["S1"] < float(Conf["MIN_SNR"][VALUE])):
                # Indicate the rejection cause
                PreproObs["RejectionCause"] = REJECTION_CAUSE["MIN_SNR_F1"]
                PreproObs["Valid"] = 0
            if(PreproObs["S2"] < float(Conf["MIN_SNR"][VALUE])):
                # Indicate the rejection cause
                PreproObs["RejectionCause"] = REJECTION_CAUSE["MIN_SNR_F2"]
                PreproObs["Valid"] = 0

        # End of if (Conf["MIN_SNR"][FLAG] == 1) \

        # If satellite shall be rejected due to Pseudorange Out-of-range (only if activated in conf)
        # --------------------------------------------------------------------------------------------------------------------
        if (Conf["MAX_PSR_OUTRNG"][FLAG] == 1):
            if (PreproObs["C1"] > float(Conf["MAX_PSR_OUTRNG"][VALUE])):
                # Indicate the rejection cause
                PreproObs["RejectionCause"] = REJECTION_CAUSE["MAX_PSR_OUTRNG_F1"]
                PreproObs["Valid"] = 0
            if (PreproObs["C2"] > float(Conf["MAX_PSR_OUTRNG"][VALUE])):
                # Indicate the rejection cause
                PreproObs["RejectionCause"] = REJECTION_CAUSE["MAX_PSR_OUTRNG_F2"]
                PreproObs["Valid"] = 0

        # End of if (Conf["MAX_PSR_OUTRNG"][FLAG] == 1)

        # 54 E27 Cycle Slip
        # Account for Cycle Slips flags
        if PrevPreproObsInfo[SatLabel]["CycleSlipDetectFlag"] == 1:
            # Indicate the rejection cause
            PreproObs["RejectionCause"] = REJECTION_CAUSE["CYCLE_SLIP"]
            PreproObs["Valid"] = 0
            PrevPreproObsInfo[SatLabel]["CycleSlipDetectFlag"] = 0

        # Build combinations
        # ----------------------------------------------------------
        # Build Geometry-free combination of phases
        PreproObs["GeomFree_P"] = (PreproObs["L2"] - PreproObs["L1"]) / (1 - GammaF1F2)
        # Build Iono-free combination of codes and phases
        PreproObs["IF_C"] = (PreproObs["C2"] - GammaF1F2 * PreproObs["C1"]) / (1 - GammaF1F2)
        PreproObs["IF_P"] = (PreproObs["L2Meters"] - GammaF1F2 * PreproObs["L1Meters"]) / (1 - GammaF1F2)

        # Hatch filter (re)initialization
        # ----------------------------------------------------------
        # If Hatch filter shall be reset
        if PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] == 1:
            # Lower Smoothing filter reset flag
            PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 0

            # Ksmooth: Time index -> is equal to 1 at the beginning and is increasing by one
            # each epoch
            PrevPreproObsInfo[SatLabel]["Ksmooth"] = 1

            # Initialize smoothed values
            PreproObs["SmoothIF"] = PreproObs["IF_C"]
            PrevPreproObsInfo[SatLabel]["PrevSmooth"] = \
                PreproObs["SmoothIF"]
            
            # Reset Prealign Offset
            PrevPreproObsInfo[SatLabel]["PrealignOffset"] = \
                                        PreproObs["IF_C"] - PreproObs["IF_P"]

            # Reinitialize some variables
            PrevPreproObsInfo[SatLabel]["PrevC1"] = Const.NAN                               # Previous Code
            PrevPreproObsInfo[SatLabel]["PrevL1"] = Const.NAN                               # Previous Phase
            PrevPreproObsInfo[SatLabel]["PrevRangeRateL1"] = Const.NAN                      # Previous Code Rate
            PrevPreproObsInfo[SatLabel]["PrevPhaseRateL1"] = Const.NAN                      # Previous Phase Rate
            PrevPreproObsInfo[SatLabel]["PrevC2"] = Const.NAN                               # Previous Code
            PrevPreproObsInfo[SatLabel]["PrevL2"] = Const.NAN                               # Previous Phase
            PrevPreproObsInfo[SatLabel]["PrevRangeRateL2"] = Const.NAN                      # Previous Code Rate
            PrevPreproObsInfo[SatLabel]["PrevPhaseRateL2"] = Const.NAN                      # Previous Phase Rate

        else:
            # Code Carrier Smoothing with a Hatch Filter
            # ----------------------------------------------------------
            # Update Smoothing iterator
            PrevPreproObsInfo[SatLabel]["Ksmooth"] = \
                    PrevPreproObsInfo[SatLabel]["Ksmooth"] + DeltaT

            # Smoothing Time computation
            # Smoothing Time is equal to the time index if the time index 
            # is lower than the Hatch filter and equal to the Hatch filter 
            # time constant otherwise
            SmoothingTime = \
            (PrevPreproObsInfo[SatLabel]["Ksmooth"] <= Conf["HATCH_TIME"]) * \
                            PrevPreproObsInfo[SatLabel]["Ksmooth"] + \
            (PrevPreproObsInfo[SatLabel]["Ksmooth"] > Conf["HATCH_TIME"]) * \
                            Conf["HATCH_TIME"]

            # Weighting factor of the Smoothing filter
            Alpha = float(DeltaT) / SmoothingTime

            # Compute Smoothed C1
            PreproObs["SmoothIF"] = \
                Alpha * PreproObs["IF_C"] + \
                (1-Alpha) * \
                    (PrevPreproObsInfo[SatLabel]["PrevSmooth"] + \
                        PreproObs["IF_P"] - PrevPreproObsInfo[SatLabel]["IF_P_Prev"])

        # End of if PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] == 1:

        # Loop over frequencies
        for freq in ['1', '2']:
            # If previous information is available
            if(PrevPreproObsInfo[SatLabel]["PrevL" + freq] != Const.NAN):
                # Check Phase Rate (only if activated in conf)
                # --------------------------------------------------------------------------------------------------------------------
                # Compute Phase Rate in meters/second
                PreproObs["PhaseRateL" + freq] = ((PreproObs["L" + freq] - PrevPreproObsInfo[SatLabel]["PrevL" + freq]) / DeltaT) * Wave["F" + freq]

                # 55560 E07 -1e6 cy in f2
                # Check Phase Rate
                if (Conf["MAX_PHASE_RATE"][FLAG] == 1) \
                        and (abs(PreproObs["PhaseRateL" + freq]) > Conf["MAX_PHASE_RATE"][VALUE]):
                            # Indicate the rejection cause
                            PreproObs["RejectionCause"] = REJECTION_CAUSE["MAX_PHASE_RATE_F" + freq]
                            PreproObs["Valid"] = 0

                            # Raise flag to reset Hatch filter
                            PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1

                # If there are enough samples
                if (PrevPreproObsInfo[SatLabel]["PrevPhaseRateL" + freq] != Const.NAN):
                    # Check Phase Rate Step (only if activated in conf)
                    # ----------------------------------------------------------
                    # Compute Phase Rate Step in meters/second^2
                    PreproObs["PhaseRateStepL" + freq] = (PreproObs["PhaseRateL" + freq] - PrevPreproObsInfo[SatLabel]["PrevPhaseRateL" + freq]) / DeltaT

                    if (Conf["MAX_PHASE_RATE_STEP"][FLAG] == 1) \
                        and (abs(PreproObs["PhaseRateStepL" + freq]) > Conf["MAX_PHASE_RATE_STEP"][VALUE]):
                            # Indicate the rejection cause
                            PreproObs["RejectionCause"] = REJECTION_CAUSE["MAX_PHASE_RATE_STEP_F" + freq]
                            PreproObs["Valid"] = 0

                            # Raise flag to reset Hatch filter
                            PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1

                else: PreproObs["Valid"] = 0

                # End of if (PrevPreproObsInfo[SatLabel]["PrevPhaseRateL" + freq] != Const.NAN):

            else: PreproObs["Valid"] = 0

            # End of if(PrevPreproObsInfo[SatLabel]["PrevL" + freq] != Const.NAN):

            # If previous information is available
            if(PrevPreproObsInfo[SatLabel]["PrevC" + freq] != Const.NAN):
                # Check Code Step (only if activated in conf)
                # --------------------------------------------------------------------------------------------------------------------
                # Compute Code Rate in meters/second
                PreproObs["RangeRateL" + freq] = (PreproObs["C" + freq] - PrevPreproObsInfo[SatLabel]["PrevC" + freq]) / DeltaT

                # Check Code Rate
                if (Conf["MAX_CODE_RATE"][FLAG] == 1) \
                        and (abs(PreproObs["RangeRateL" + freq]) > Conf["MAX_CODE_RATE"][VALUE]):
                            # Indicate the rejection cause
                            PreproObs["RejectionCause"] = REJECTION_CAUSE["MAX_CODE_RATE_F" + freq]
                            PreproObs["Valid"] = 0

                            # Raise flag to reset Hatch filter
                            PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1

                # If there are enough samples
                if (PrevPreproObsInfo[SatLabel]["PrevRangeRateL" + freq] != Const.NAN):
                    # Compute Code Rate Step in meters/second^2
                    PreproObs["RangeRateStepL" + freq] = (PreproObs["RangeRateL" + freq] - \
                        PrevPreproObsInfo[SatLabel]["PrevRangeRateL" + freq]) / DeltaT
                    # Check Code Rate Step (only if activated in conf)
                    # ----------------------------------------------------------
                    if (Conf["MAX_CODE_RATE_STEP"][FLAG] == 1) \
                            and (abs(PreproObs["RangeRateStepL" + freq]) > Conf["MAX_CODE_RATE_STEP"][VALUE]):
                                # Indicate the rejection cause
                                PreproObs["RejectionCause"] = REJECTION_CAUSE["MAX_CODE_RATE_STEP_F" + freq]
                                PreproObs["Valid"] = 0

                                # Raise flag to reset Hatch filter
                                PrevPreproObsInfo[SatLabel]["ResetHatchFilter"] = 1

                else: PreproObs["Valid"] = 0

                # End of if (PrevPreproObsInfo[SatLabel]["PrevRangeRateL" + freq] != Const.NAN):

            else: PreproObs["Valid"] = 0

            # End of if(PrevPreproObsInfo[SatLabel]["PrevC" + freq] != Const.NAN)

        # End of for freq in ['1', '2']:

        # Set Status flag
        # ----------------------------------------------------------
        # 1 if convergence was reached, 0 otherwise
        if(PrevPreproObsInfo[SatLabel]["Ksmooth"] > \
            Conf["HATCH_STATE_F"] * Conf["HATCH_TIME"]) and \
              (PreproObs["Valid"] != 0) :
            PreproObs["Status"] = 1
        else: 
            PreproObs["Status"] = 0

        # Update previous values
        # ----------------------------------------------------------
        PrevPreproObsInfo[SatLabel]["PrevC1"] = PreproObs["C1"]
        PrevPreproObsInfo[SatLabel]["PrevL1"] = PreproObs["L1"]
        PrevPreproObsInfo[SatLabel]["PrevC2"] = PreproObs["C2"]
        PrevPreproObsInfo[SatLabel]["PrevL2"] = PreproObs["L2"]
        PrevPreproObsInfo[SatLabel]["PrevSmooth"] = PreproObs["SmoothIF"]
        PrevPreproObsInfo[SatLabel]["IF_P_Prev"] = PreproObs["IF_P"]
        PrevPreproObsInfo[SatLabel]["PrevRangeRateL1"] = PreproObs["RangeRateL1"]
        PrevPreproObsInfo[SatLabel]["PrevRangeRateL2"] = PreproObs["RangeRateL2"]
        PrevPreproObsInfo[SatLabel]["PrevPhaseRateL1"] = PreproObs["PhaseRateL1"]
        PrevPreproObsInfo[SatLabel]["PrevPhaseRateL2"] = PreproObs["PhaseRateL2"]
        PrevPreproObsInfo[SatLabel]["PrevRej"] = PreproObs["RejectionCause"]
        PrevPreproObsInfo[SatLabel]["PrevEpoch"] = Epoch

        # Pre-align the Phase
        PreproObs["IF_P"] +=  PrevPreproObsInfo[SatLabel]["PrealignOffset"]

    # End of for SatLabel, PreproObs in PreproObsInfo.items():

    return PreproObsInfo

# End of function runPreprocessing()

########################################################################
# END OF PREPROCESSING FUNCTIONS MODULE
########################################################################
