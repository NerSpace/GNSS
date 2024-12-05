#!/usr/bin/env python

########################################################################
# SENTUS-V3/SRC/InputOutput.py:
# This is the Inputs (conf and input files) Module of SENTUS tool
#
#  Project:        SENTUS
#  File:           Preprocessing.py
#
#   Author: Nerea Sánchez
#   Copyright 2024 GNSS Academy / Nerea Sánchez
#
########################################################################


# Import External and Internal functions and Libraries
#----------------------------------------------------------------------
import sys, os
from collections import OrderedDict
from COMMON.Dates import convertYearMonthDay2JulianDay
from COMMON import GnssConstants as Const
from COMMON.Coordinates import llh2xyz
import numpy as np

# Input interfaces
#----------------------------------------------------------------------
# CONF
FLAG = 0
VALUE = 1
TH = 1
CSNEPOCHS = 2
CSNPOINTS = 3
CSPDEGREE = 4

# OBS file columns
ObsIdxP = OrderedDict({})
ObsIdxP["SOD"]=1
ObsIdxP["PRN"]=2
ObsIdxP["L1"]=3
ObsIdxP["L2"]=4

ObsIdxC = OrderedDict({})
ObsIdxC["SOD"]=1
ObsIdxC["PRN"]=2
ObsIdxC["ELEV"]=3
ObsIdxC["AZIM"]=4
ObsIdxC["C1"]=5
ObsIdxC["C2"]=6
ObsIdxC["S1"]=7
ObsIdxC["S2"]=8

# LEO POS file columns
LeoPosIdx = OrderedDict({})
LeoPosIdx["SOD"]=0
LeoPosIdx["DOY"]=1
LeoPosIdx["YEAR"]=2
LeoPosIdx["xCM"]=3
LeoPosIdx["yCM"]=4
LeoPosIdx["zCM"]=5

# LEO QUAT file columns
LeoQuatIdx = OrderedDict({})
LeoQuatIdx["SOD"]=0
LeoQuatIdx["q0"]=1
LeoQuatIdx["q1"]=2
LeoQuatIdx["q2"]=3
LeoQuatIdx["q3"]=4

# SAT POS file columns
SatPosIdx = OrderedDict({})
SatPosIdx["SOD"]=0
SatPosIdx["DOY"]=1
SatPosIdx["YEAR"]=2
SatPosIdx["CONST"]=3
SatPosIdx["PRN"]=4
SatPosIdx["xCM"]=5
SatPosIdx["yCM"]=6
SatPosIdx["zCM"]=7
SatPosIdx["SatLabel"]=8 # This column is added after being read

# SAT APO file columns
SatApoIdx = OrderedDict({})
SatApoIdx["CONST"]=0
SatApoIdx["PRN"]=1
SatApoIdx["x_f1"]=2
SatApoIdx["y_f1"]=3
SatApoIdx["z_f1"]=4
SatApoIdx["x_f2"]=5
SatApoIdx["y_f2"]=6
SatApoIdx["z_f2"]=7
SatApoIdx["SatLabel"]=8 # This column is added after being read

# SAT CLK file columns
SatClkIdx = OrderedDict({})
SatClkIdx["SOD"]=0
SatClkIdx["DOY"]=1
SatClkIdx["YEAR"]=2
SatClkIdx["CONST"]=3
SatClkIdx["PRN"]=4
SatClkIdx["CLK-BIAS"]=5
SatClkIdx["SatLabel"]=6 # This column is added after being read

# SAT BIA file columns
SatBiaIdx = OrderedDict({})
SatBiaIdx["CONST"]=0
SatBiaIdx["PRN"]=1
SatBiaIdx["CLK_f1_C"] = 2
SatBiaIdx["CLK_f2_C"] = 3
SatBiaIdx["OBS_f1_C"] = 4
SatBiaIdx["OBS_f2_C"] = 5
SatBiaIdx["CLK_f1_P"] = 6
SatBiaIdx["CLK_f2_P"] = 7
SatBiaIdx["OBS_f1_P"] = 8
SatBiaIdx["OBS_f2_P"] = 9
SatBiaIdx["SatLabel"]=10 # This column is added after being read

# Output interfaces
#----------------------------------------------------------------------
# PREPRO OBS 
# Header
PreproHdr = "\
# SOD   PRN    ELEV     AZIM  VALID REJ  STATUS             C1             C2               L1              L2       S1       S2   CODERATE   CODEACC   PHASERATE   PHASEACC          CODEIF         PHASEIF        SMOOTHIF\n"

# Line format
PreproFmt = "%6d %6s %8.3f %8.3f %4d %4d %4d "\
    "%15.3f %15.3f %15.3f %15.3f %8.3f %8.3f %10.3f %10.3f %10.3f %10.3f "\
    "%15.3f %15.3f %15.3f".split()

# File columns
PreproIdx = OrderedDict({})
PreproIdx["SOD"]=0
PreproIdx["PRN"]=1
PreproIdx["ELEV"]=2
PreproIdx["AZIM"]=3
PreproIdx["VALID"]=4
PreproIdx["REJECT"]=5
PreproIdx["STATUS"]=6
PreproIdx["C1"]=7
PreproIdx["C2"]=8
PreproIdx["L1"]=9
PreproIdx["L2"]=10
PreproIdx["S1"]=11
PreproIdx["S2"]=12
PreproIdx["CODE_RATE"]=13
PreproIdx["CODE_RATE_STEP"]=14
PreproIdx["PHASE_RATE"]=15
PreproIdx["PHASE_RATE_STEP"]=16
PreproIdx["CODE_IF"]=17
PreproIdx["PHASE_IF"]=18
PreproIdx["SMOOTH_IF"]=19

# Rejection causes flags
REJECTION_CAUSE = OrderedDict({})
REJECTION_CAUSE["MASKANGLE"]=1
REJECTION_CAUSE["DATA_GAP"]=2
REJECTION_CAUSE["MIN_SNR_F1"]=3
REJECTION_CAUSE["MIN_SNR_F2"]=4
REJECTION_CAUSE["MAX_PSR_OUTRNG_F1"]=5
REJECTION_CAUSE["MAX_PSR_OUTRNG_F2"]=6
REJECTION_CAUSE["MAX_PHASE_RATE_F1"]=7
REJECTION_CAUSE["MAX_PHASE_RATE_F2"]=8
REJECTION_CAUSE["MAX_PHASE_RATE_STEP_F1"]=9
REJECTION_CAUSE["MAX_PHASE_RATE_STEP_F2"]=10
REJECTION_CAUSE["MAX_CODE_RATE_F1"]=11
REJECTION_CAUSE["MAX_CODE_RATE_F2"]=12
REJECTION_CAUSE["MAX_CODE_RATE_STEP_F1"]=13
REJECTION_CAUSE["MAX_CODE_RATE_STEP_F2"]=14
REJECTION_CAUSE["CYCLE_SLIP"]=15

REJECTION_CAUSE_DESC = OrderedDict({})
REJECTION_CAUSE_DESC["1: Mask Angle"]=1
REJECTION_CAUSE_DESC["2: Data Gap"]=2
REJECTION_CAUSE_DESC["3: Minimum C/N0 in f1"]=3
REJECTION_CAUSE_DESC["4: Minimum C/N0 in f2"]=4
REJECTION_CAUSE_DESC["5: Maximum PR in f1"]=5
REJECTION_CAUSE_DESC["6: Maximum PR in f2"]=6
REJECTION_CAUSE_DESC["7: Maximum Phase Rate in f1"]=7
REJECTION_CAUSE_DESC["8: Maximum Phase Rate in f2"]=8
REJECTION_CAUSE_DESC["9: Maximum Phase Rate Step in f1"]=9
REJECTION_CAUSE_DESC["10: Maximum Phase Rate Step in f2"]=10
REJECTION_CAUSE_DESC["11: Maximum Code Rate in f1"]=11
REJECTION_CAUSE_DESC["12: Maximum Code Rate in f2"]=12
REJECTION_CAUSE_DESC["13: Maximum Code Rate Step in f1"]=13
REJECTION_CAUSE_DESC["14: Maximum Code Rate Step in f2"]=14
REJECTION_CAUSE_DESC["15: Cycle Slip"]=15

# CORR 
# Header
CorrHdr = "\
#SOD C PRN   ELEV    AZIM     FLAG  LEO-X         LEO-Y           LEO-Z  LEO-APO-X LEO-APO-Y LEO-APO-Z  SAT-X         SAT-Y           SAT-Z  SAT-APO-X SAT-APO-Y SAT-APO-Z         SAT-CLK CODE-BIA PHASE-BIA     TOF          DTR       CORR-CODE     CORR-PHASE       GEOM-RNGE   CODE-RES  PHASE-RES       RCVR-CLK       SUERE   \n"

# Line format
CorrFmt = "%05d %1s %2s %8.3f %8.3f %4d %14.3f %14.3f %14.3f %8.3f %8.3f %8.3f \
%14.3f %14.3f %14.3f %8.3f %8.3f %8.3f %14.3f %8.3f %8.3f %8.3f %14.3f \
%14.3f %14.3f %14.3f %10.4f %10.4f %14.3f \
%10.4f".split()

# File columns
CorrIdx = OrderedDict({})
CorrIdx["SOD"]=0
CorrIdx["CONST"]=1
CorrIdx["PRN"]=2
CorrIdx["ELEV"]=3
CorrIdx["AZIM"]=4
CorrIdx["FLAG"]=5
CorrIdx["LEO-X"]=6
CorrIdx["LEO-Y"]=7
CorrIdx["LEO-Z"]=8
CorrIdx["LEO-APO-X"]=9
CorrIdx["LEO-APO-Y"]=10
CorrIdx["LEO-APO-Z"]=11
CorrIdx["SAT-X"]=12
CorrIdx["SAT-Y"]=13
CorrIdx["SAT-Z"]=14
CorrIdx["SAT-APO-X"]=15
CorrIdx["SAT-APO-Y"]=16
CorrIdx["SAT-APO-Z"]=17
CorrIdx["SAT-CLK"]=18
CorrIdx["SAT-CODE-BIA"]=19
CorrIdx["SAT-PHASE-BIA"]=20
CorrIdx["FLIGHT-TIME"]=21
CorrIdx["DTR"]=22
CorrIdx["CORR-CODE"]=23
CorrIdx["CORR-PHASE"]=24
CorrIdx["GEOM-RNGE"]=25
CorrIdx["CODE-RES"]=26
CorrIdx["PHASE-RES"]=27
CorrIdx["RCVR-CLK"]=28
CorrIdx["SUERE"]=29

# POS-PVT
# Header
PosHdr = "\
# SOD     LONEST     LATEST        ALTEST       CLKEST      GGTO    SOL NSVIS NSV   HPE     VPE      EPE     NPE     UPE    HDOP      VDOP    PDOP   \n"

# Line format
PosFmt = "%05d %11.5f %10.5f %14.3f %10.3f %10.3f %4d %4d %4d "\
         "%7.3f %7.3f %7.3f %7.3f %7.3f %8.3f %8.3f %8.3f".split()

# File columns
PosIdx = OrderedDict({})
PosIdx["SOD"] = 0
PosIdx["LONEST"] = 1
PosIdx["LATEST"] = 2
PosIdx["ALTEST"] = 3
PosIdx["CLKEST"] = 4
PosIdx["GGTO"] = 5
PosIdx["SOL"] = 6
PosIdx["NSVVIS"] = 7
PosIdx["NSV"] = 8
PosIdx["HPE"] = 9
PosIdx["VPE"] = 10
PosIdx["EPE"] = 11
PosIdx["NPE"] = 12
PosIdx["UPE"] = 13
PosIdx["HDOP"] = 14
PosIdx["VDOP"] = 15
PosIdx["PDOP"] = 16

# PERF
# Header
PerfHdr = "\
#RCVR SAMPS NOSOL NSVMIN NSVMAX     HPERMS     VPERMS      HPE95      VPE95     HPEMAX     VPEMAX   PDOPMAX    HDOPMAX    VDOPMAX \n"

# Line format
PerfFmt = "%4s %6d %6d %6d %6d %10.3f %10.3f %10.3f %10.3f "\
         "%10.3f %10.3f %10.3f %10.3f %10.3f".split()

# File columns
PerfIdx = OrderedDict({})
PerfIdx["RCVR"] = 0
PerfIdx["SAMPS"] = 1
PerfIdx["NOSOL"] = 2
PerfIdx["NSVMIN"] = 3
PerfIdx["NSVMAX"] = 4
PerfIdx["HPERMS"] = 5
PerfIdx["VPERMS"] = 6
PerfIdx["HPE95"] = 7
PerfIdx["VPE95"] = 8
PerfIdx["HPEMAX"] = 9
PerfIdx["VPEMAX"] = 10
PerfIdx["PDOPMAX"] = 11
PerfIdx["HDOPMAX"] = 12
PerfIdx["VDOPMAX"] = 13

# Input functions
#----------------------------------------------------------------------
def checkConfParam(Key, Fields, MinFields, MaxFields, LowLim, UppLim):

    '''
    
    Purpose: check configuration parameter format, type and range

    Parameters
    ==========
    Key: str
            Configuration parameter key
    Fields: list
            Configuration parameter read from conf and split
    MinFields: int
            Minimum number of fields expected
    MaxFields: int
            Maximum number of fields expected
    LowLim: list
            List containing lower limit allowed for each of the fields
    UppLim: list
            List containing upper limit allowed for each of the fields

    Returns
    =======
    Values: str, int, float or list
            Configuration parameter value or list of values

    '''
    
    # Prepare output list
    Values = []

    # Get Fields length
    LenFields = len(Fields) - 1

    # Check that number of fields is not less than the expected minimum
    if(LenFields < MinFields):
        # Display an error
        sys.stderr.write("ERROR: Too few fields (%d) for configuration parameter %s. "\
        "Minimum = %d\n" % (LenFields, Key, MinFields))
        sys.exit(-1)
    # End if(LenFields < MinFields)

    # Check that number of fields is not greater than the expected minimum
    if(LenFields > MaxFields):
        # Display an error
        sys.stderr.write("ERROR: Too many fields (%d) for configuration parameter %s. "\
        "Maximum = %d\n" % (LenFields, Key, MaxFields))
        sys.exit(-1)
    # End if(LenFields > MaxFields)

    # Loop over fields
    for i, Field in enumerate(Fields[1:]):
        # If float
        try:
            # Convert to float and append to the outputs
            Values.append(float(Field))

        except:
            # isnumeric check
            try:
                Check = unicode(Field).isnumeric()

            except:
                Check = (Field).isnumeric()

            # If it is integer
            if(Check):
                # Convert to int and append to the outputs
                Values.append(int(Field))

            else:
                # Append to the outputs
                Values.append(Field)

    # End of for i, Field in enumerate(Fields[1:]):

    # Loop over values to check the range
    for i, Field in enumerate(Values):
        # If range shall be checked
        if(isinstance(LowLim[i], int) or \
            isinstance(LowLim[i], float)):
            # Try to check the range
            try:
                if(Field<LowLim[i] or Field>UppLim[i]):
                    # Out of range
                    sys.stderr.write("ERROR: Configuration parameter %s "\
                        "%f is out of range [%f, %f]\n" % 
                    (Key, Field, LowLim[i], UppLim[i]))

            except:
                # Wrong format
                sys.stderr.write("ERROR: Wrong type for configuration parameter %s\n" %
                Key)
                sys.exit(-1)

    # End of for i, Field in enumerate(Values):

    # If only one element, return the value directly
    if len(Values) == 1:
        return Values[0]

    # Otherwise, return the list
    else:
        return Values

# End of checkConfParam()

def readConf(CfgFile):

    '''
    Purpose: read the configuration file
    
    Parameters
    ==========
    CfgFile: str
            Path to conf file

    Returns
    =======
    Conf: Dict
            Conf loaded in a dictionary    
    '''

    # Function to check the format of the dates
    def checkConfDate(Key, Fields):
        # Split Fields
        FieldsSplit=Fields[1].split('/')

        # Set expected number of characters
        ExpectedNChar = [2,2,4]

        # Check the number of characters in each field
        for i, Field in enumerate(FieldsSplit):
            # if number of characters is incorrect
            if len(Field) != ExpectedNChar[i]:
                sys.stderr.write("ERROR: wrong format in configured %s\n" % Key)
                sys.exit(-1)

    # End of checkConfDate()

    # Initialize the variable to store the conf
    Conf = OrderedDict({})

    # Initialize the configuration parameters counter
    NReadParams = 0
    
    # Open the file
    with open(CfgFile, 'r') as f:
        # Read file
        Lines = f.readlines()

        # Parse each Line of configuration file
        for Line in Lines:
            
            # Check if Line is not a comment
            if Line[0]!='#':
                
                # Split Line in a list of Fields
                Fields=Line.rstrip('\n').split(' ')
                
                # Check if line is blank
                if '' in Fields:
                    Fields = list(filter(None, Fields))

                if Fields != None :
                    # if some parameter with its value missing, warn the user
                    if len(Fields) == 1:
                        sys.stderr.write("ERROR: Configuration file contains a parameter" \
                            "with no value: " + Line)
                        sys.exit(-1)

                    # if the line contains a conf parameter
                    elif len(Fields)!=0:
                        # Get conf parameter key
                        Key=Fields[0]

                        # Fill in Conf appropriately according to the configuration file Key
                        
                        # Scenario Start and End Dates [GPS time in Calendar format]
                        #--------------------------------------------------------------------
                        # Date format DD/MM/YYYY (e.g: 01/09/2019)
                        #--------------------------------------------------------------------
                        if Key=='INI_DATE':
                            # Check date format
                            checkConfDate(Key, Fields)

                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, [None], [None])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        elif Key=='END_DATE':
                            # Check date format
                            checkConfDate(Key, Fields)

                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, [None], [None])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Scenario Sampling Rate [SECONDS]
                        #-------------------------------------------
                        elif Key=='SAMPLING_RATE':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, [1], [Const.S_IN_D])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Navigation Solution Selection
                        #-----------------------------------------------
                        # Three Options:
                        #       GPS: SBAS GPS
                        #       GAL: SBAS Galileo
                        #       GPSGAL: SBAS GPS+Galileo
                        #-----------------------------------------------
                        elif Key=='NAV_SOLUTION':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, [None], [None])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Preprocessing outputs selection [0:OFF|1:ON]
                        #--------------------------------------------------------------------       
                        elif Key=='PREPRO_OUT':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, [0], [1])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Corrected outputs selection [0:OFF|1:ON]
                        #--------------------------------------------------------------------       
                        elif Key=='CORR_OUT':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, [0], [1])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # KML LEO Coordinates output selection [0:OFF|1:ON]
                        #--------------------------------------------------------------------       
                        elif Key=='KML_OUT':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, [0], [1])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Satellite ACRONYM
                        #-----------------------------------------------
                        elif Key=='SAT_ACRONYM':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, [None], [None])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Satellite Reference Positions
                        #-----------------------------------------------
                        elif Key=='SAT_POS':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, [None], [None])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # RIMS positions file Name  (if RCVR_INFO=STATIC)
                        #-----------------------------------------------
                        elif Key=='RCVR_FILE':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, [None], [None])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # RCVR mask Angle [DEG]
                        #-----------------------------------------------
                        elif Key== 'RCVR_MASK': 
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, 
                            [Const.MIN_MASK_ANGLE], [Const.MAX_MASK_ANGLE])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Minimum Carrier To Noise Ratio
                        #------------------------------
                        # p1: Check C/No [0:OFF|1:ON]
                        # p2: C/No Threshold [dB-Hz]
                        #------------------------------
                        elif Key== 'MIN_SNR':  
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 2, 2, 
                            [0, 0], [1, 80])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Check Cycle Slips 
                        #----------------------------------------
                        # p1: Check CS [0:OFF|1:ON]
                        # p2: CS threshold [cycles]
                        # p3: CS Nepoch
                        #----------------------------------------
                        elif Key== 'CYCLE_SLIPS':  
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 5, 5, 
                            [0, 0, 1, 0, 1], [1, 10, 10, 100, 10])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1
                        
                        # Check Pseudo-Range Measurement Out of Range
                        #-------------------------------------------
                        # p1: Check PSR Range [0:OFF|1:ON]
                        # p2: Max. Range [m]  (Default:330000000]
                        #-----------------------------------------------
                        elif Key== 'MAX_PSR_OUTRNG':  
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 2, 2, 
                            [0, 0], 
                            [1, 400000000])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1
                        
                        # Check Code Rate
                        #-----------------------------------------------
                        # p1: Check Code Rate [0:OFF|1:ON]
                        # p2: Max. Code Rate [m/s]  (Default: 952)
                        #-----------------------------------------------
                        elif Key== 'MAX_CODE_RATE':  
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 2, 2, 
                            [0, 0], 
                            [1, 9000])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Check Code Rate Step 
                        #-----------------------------------------------
                        # p1: Check Code Rate Step [0:OFF|1:ON]
                        # p2: Max. Code Rate Step [m/s**2]  (Default: 10)
                        #-----------------------------------------------
                        elif Key== 'MAX_CODE_RATE_STEP':  
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 2, 2, 
                            [0, 0], 
                            [1, 100])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Check Phase Measurement Rate 
                        #-----------------------------------------------
                        # p1: Check Phase Rate [0:OFF|1:ON]
                        # p2: Max. Phase Rate [m/s]  (Default: 952)
                        #-----------------------------------------------
                        elif Key== 'MAX_PHASE_RATE':  
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 2, 2, 
                            [0, 0], 
                            [1, 9000])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1
                            
                        # Check Phase Rate Step 
                        #-----------------------------------------------
                        # p1: Check Phase Rate Step [0:OFF|1:ON]
                        # p2: Max. Phase Rate Step [m/s**2]  (Default: 10 m/s**2)
                        #-----------------------------------------------
                        elif Key== 'MAX_PHASE_RATE_STEP':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 2, 2, 
                            [0, 0], 
                            [1, 100])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Max. DATA GAP for PSR Propagation reset [s]
                        #------------------------------------------
                        elif Key== 'MAX_DATA_GAP':  
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 2, 2, 
                            [0, 0], [1, 3600])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Hatch filter Smoothing time [s]
                        #----------------------------------
                        elif Key== 'HATCH_TIME':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, 
                            [0], 
                            [3600])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Hatch filter Steady State factor
                        #----------------------------------
                        elif Key== 'HATCH_STATE_F':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, 
                            [0], [10])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # LEO CoM position in SRF [m]
                        #----------------------------------
                        elif Key== 'LEO_COM_POS':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 3, 3, 
                            [-10]*3, [10]*3)

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # LEO ARP position in SRF [m]
                        #----------------------------------
                        elif Key== 'LEO_ARP_POS':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 3, 3, 
                            [-10]*3, [10]*3)

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # LEO GPS PCO in SRF [m]
                        #----------------------------------
                        elif Key== 'LEO_PCO_GPS':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 3, 3, 
                            [-10]*3, [10]*3)

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # LEO Galileo PCO in SRF [m]
                        #----------------------------------
                        elif Key== 'LEO_PCO_GAL':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 3, 3, 
                            [-10]*3, [10]*3)

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Satellite APO file
                        #-----------------------------------------------
                        elif Key== 'SAT_APO_FILE': 
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, [None], [None])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Satellite Signal Biases file
                        #-----------------------------------------------
                        elif Key== 'SAT_BIA_FILE': 
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, [None], [None])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # GPS satellites UERE [m]
                        #----------------------------------
                        elif Key== 'GPS_UERE':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, 
                            [0], 
                            [100])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Galileo satellites UERE [m]
                        #----------------------------------
                        elif Key== 'GAL_UERE':
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, 
                            [0], 
                            [100])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Max. Number of interations for Navigation Solution
                        #----------------------------------------------------
                        elif Key== 'MAX_LSQ_ITER': 
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, 
                            [0], [1e8])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

                        # Maximum PDOP Threshold for Solution [m]
                        # Default Value: 10000.0
                        #-----------------------------------------------
                        elif Key== 'PDOP_MAX': 
                            # Check parameter and load it in Conf
                            Conf[Key] = checkConfParam(Key, Fields, 1, 1, 
                            [0], [Const.MAX_PDOP_PVT])

                            # Increment number of read parameters
                            NReadParams = NReadParams + 1

    return Conf

# End of readConf()

def processConf(Conf):
    '''
    Purpose: process the configuration
    
    Parameters
    ==========
    Conf: dict
            Dictionary containing configuration

    Returns
    =======
    Conf: dict
            Dictionary containing configuration with
            Julian Days
    '''

    ConfCopy = Conf.copy()
    for Key in ConfCopy:
        Value = ConfCopy[Key]
        if Key == "INI_DATE" or Key == "END_DATE":
            ParamSplit = Value.split('/')

            # Compute Julian Day
            Conf[Key + "_JD"] = \
                int(round(
                    convertYearMonthDay2JulianDay(
                        int(ParamSplit[2]),
                        int(ParamSplit[1]),
                        int(ParamSplit[0]))
                    )
                )

    return Conf


def splitLine(Line):

    '''
    Purpose: split line
    
    Parameters
    ==========
    Line: str
            string containing line read from file

    Returns
    =======
    CfgFile: list
            line split using spaces as delimiter
    '''
    
    
    LineSplit = Line.split()

    return LineSplit

# End of splitLine()


def readObsEpoch(f):
    
    '''
    Purpose: read one epoch of OBS file (all the LoS)
    
    Parameters
    ==========
    f: file descriptor
            OBS file

    Returns
    =======
    EpochObsC: list
            list of the split lines for Codes
    EpochObsP: list
            list of the split lines for Phases
            EpochInfoC[1][1] is the second field of the 
            second line  
    '''

    EpochObsC = []
    EpochObsP = []
    
    # Read one line
    Line = f.readline()
    if(not Line):
        return []
    LineSplit = splitLine(Line)
    Sod = LineSplit[ObsIdxC["SOD"]]
    SodNext = Sod

    while SodNext == Sod:
        if (LineSplit[0] == 'C'):
            EpochObsC.append(LineSplit)
        else:
            EpochObsP.append(LineSplit)
        Pointer = f.tell()
        Line = f.readline()
        LineSplit = splitLine(Line)
        try: 
            SodNext = LineSplit[ObsIdxC["SOD"]]

        except:
            return EpochObsC, EpochObsP

    f.seek(Pointer)

    return EpochObsC, EpochObsP

# End of readObsEpoch()


def openInputFile(Path):

    '''
    Purpose: check existence and open input file
    
    Parameters
    ==========
    Path: str
            Path to file

    Returns
    =======
    f: File descriptor
            Descriptor of the input file
    
    '''

    # Display Message
    print("INFO: Reading file: %s..." %
    Path)

    # Try to open the file
    try:
        # Open PREPRO OBS file
        f = open(Path, 'r')

        # Read header line
        f.readline()

    # If file could not be opened
    except:
        # Display error
        sys.stderr.write("ERROR: In input file: %s...\n" %
        Path)

    return f

# End of openInputFile()

####### READ DATA FILES ###########

def readLeoDatFile(f):
    '''
    Purpose: read Leo data files
    
    Parameters
    ==========
    f: file descriptor

    Returns
    =======
    dict with info from file    
    '''

    data = {}

        # Read the data
    for line in f:
        # Split the line
        lineSplit = splitLine(line)

        # Convert the values to floats
        lineValues = [float(value) for value in lineSplit]

        # Save the values to the dictionary
        for i, value in enumerate(lineValues):
            if i not in data:
                data[i] = []
            data[i].append(value)

    # Close the file
    f.close()

    return data

# End of readLeoDatFile()

def readSatDatFile(f, prnindex):
    '''
    Purpose: read satellite data files

    Parameters
    ==========
    f: file descriptor
    prnindex: int
            index of the PRN column
    lastindex: int
            index of the last column to read

    Returns
    =======
    dict with info from file   
    '''

    data = {}

    # Read the data
    for line in f:
        # Split the line
        lineSplit = splitLine(line)

        # Convert the values to floats, except for the PRN column and the column after it
        lineValues = []
        for i, value in enumerate(lineSplit):
            if i == prnindex or i == prnindex + 1:
                lineValues.append(value)
            else:
                lineValues.append(float(value))

        # Add a SatLabel column that is the concatenation of the PRN column and the column after it
        satLabel = lineSplit[prnindex] + lineSplit[prnindex + 1]
        lineValues.append(satLabel)

        # Save the values to the dictionary
        for i, value in enumerate(lineValues):
            if i not in data:
                data[i] = []
            data[i].append(value)

    # Close the file
    f.close()

    return data
# End of readLeoDatFile()

####### GENERATE FILES ###########

def createOutputFile(Path, Hdr):
    '''
    Purpose: open output file and write its header
    
    Parameters
    ==========
    Path: str
            Path to file
    Hdr: str
            File header

    Returns
    =======
    f: File descriptor
            Descriptor of output file
    '''
    # Display Message
    print("INFO: Creating file: %s..." % Path)

    # Create output directory, if needed
    if not os.path.exists(os.path.dirname(Path)):
        os.makedirs(os.path.dirname(Path))

    # Open file
    f = open(Path, 'w')

    # Write header
    f.write(Hdr)

    return f

# End of createOutputFile()

def generatePreproFile(fpreprobs, PreproObsInfo):
    '''
    Purpose: generate output file with Preprocessing results

    Parameters
    ==========
    fpreprobs: file descriptor
            Descriptor for PREPRO OBS output file
    PreproObsInfo: dict
            Dictionary containing Preprocessing info for the 
            current epoch

    Returns
    =======
    Nothing    
    '''


    # Loop over satellites
    for SatLabel, SatPreproObs in PreproObsInfo.items():
        # Prepare outputs
        Outputs = OrderedDict({})
        Outputs["SOD"] = SatPreproObs["Sod"]
        Outputs["PRN"] = SatLabel
        Outputs["ELEV"] = SatPreproObs["Elevation"]
        Outputs["AZIM"] = SatPreproObs["Azimuth"]
        Outputs["VALID"] = SatPreproObs["Valid"]
        Outputs["REJECT"] = SatPreproObs["RejectionCause"]
        Outputs["STATUS"] = SatPreproObs["Status"]
        Outputs["C1"] = SatPreproObs["C1"]
        Outputs["C2"] = SatPreproObs["C2"]
        Outputs["L1"] = SatPreproObs["L1Meters"]
        Outputs["L2"] = SatPreproObs["L2Meters"]
        Outputs["S1"] = SatPreproObs["S1"]
        Outputs["S2"] = SatPreproObs["S2"]
        Outputs["CODE_RATE"] = SatPreproObs["RangeRateL1"]
        Outputs["CODE_RATE_STEP"] = SatPreproObs["RangeRateStepL1"]
        Outputs["PHASE_RATE"] = SatPreproObs["PhaseRateL1"]
        Outputs["PHASE_RATE_STEP"] = SatPreproObs["PhaseRateStepL1"]
        Outputs["CODE_IF"] = SatPreproObs["IF_C"]
        Outputs["PHASE_IF"] = SatPreproObs["IF_P"]
        Outputs["SMOOTH_IF"] = SatPreproObs["SmoothIF"]
        # Outputs["GEOMFREE"] = SatPreproObs["GeomFree_P"]

        # Write line
        for i, result in enumerate(Outputs):
            fpreprobs.write(((PreproFmt[i] + " ") % Outputs[result]))

        fpreprobs.write("\n")

# End of generatePreproFile()

def generateCorrFile(fcorr, CorrInfo): 
    '''
    # Purpose: generate output file with Corrected results

    # Parameters
    # ==========
    # fpreprobs: file descriptor
    #         Descriptor for CORR output file
    # PreproObsInfo: dict
    #         Dictionary containing Corrected info for the 
    #         current epoch

    # Returns
    # =======
    # Nothing    
    '''

    # Loop over satellites
    for SatLabel, SatCorr in CorrInfo.items():
        # Prepare outputs
        Outputs = OrderedDict({})
        Outputs["SOD"] = SatCorr["Sod"]
        Outputs["CONST"]= SatLabel[0]
        Outputs["PRN"]= SatLabel[1:3]
        Outputs["ELEV"]= SatCorr["Elevation"]
        Outputs["AZIM"]=SatCorr["Azimuth"]
        Outputs["FLAG"]=SatCorr["Flag"]
        Outputs["LEO-X"]=SatCorr["LeoX"]
        Outputs["LEO-Y"]=SatCorr["LeoY"]
        Outputs["LEO-Z"]=SatCorr["LeoZ"]
        Outputs["LEO-APO-X"]=SatCorr["LeoApoX"]
        Outputs["LEO-APO-Y"]=SatCorr["LeoApoY"]
        Outputs["LEO-APO-Z"]=SatCorr["LeoApoZ"]
        Outputs["SAT-X"]=SatCorr["SatX"]
        Outputs["SAT-Y"]=SatCorr["SatY"]
        Outputs["SAT-Z"]=SatCorr["SatZ"]
        Outputs["SAT-APO-X"]=SatCorr["SatApoX"]
        Outputs["SAT-APO-Y"]=SatCorr["SatApoY"]
        Outputs["SAT-APO-Z"]=SatCorr["SatApoZ"]
        Outputs["SAT-CLK"]=SatCorr["SatClk"]
        Outputs["SAT-CODE-BIA"]=SatCorr["SatCodeBia"]
        Outputs["SAT-PHASE-BIA"]=SatCorr["SatPhaseBia"]
        Outputs["FLIGHT-TIME"]=SatCorr["FlightTime"]
        Outputs["DTR"]=SatCorr["Dtr"]
        Outputs["CORR-CODE"]=SatCorr["CorrCode"]
        Outputs["CORR-PHASE"]=SatCorr["CorrPhase"]
        Outputs["GEOM-RNGE"]=SatCorr["GeomRange"]
        Outputs["CODE-RES"]=SatCorr["CodeResidual"]
        Outputs["PHASE-RES"]=SatCorr["PhaseResidual"]
        Outputs["RCVR-CLK"]=SatCorr["RcvrClk"]
        Outputs["SUERE"]=SatCorr["SigmaUere"]

        # Write line
        for i, result in enumerate(Outputs):
            fcorr.write(((CorrFmt[i] + " ") % Outputs[result]))

        fcorr.write("\n")

# End of generateCorrFile

def generatePosFile(fpos, PosInfo): 
    '''
    # Purpose: generate output file with PVT solution results

    # Parameters
    # ==========
    # fpos: file descriptor
    #         Descriptor for POS output file
    # PosInfo: dict
    #         Dictionary containing Position info for the 
    #         current epoch

    # Returns
    # =======
    # Nothing    
    '''

    # Prepare outputs
    Outputs = OrderedDict({})
    Outputs["SOD"] = PosInfo["SOD"]
    Outputs["LONEST"] = PosInfo["LONEST"]
    Outputs["LATEST"] = PosInfo["LATEST"]
    Outputs["ALTEST"] = PosInfo["ALTEST"]
    Outputs["CLKEST"] = PosInfo["CLKEST"]
    Outputs["GGTO"] = PosInfo["GGTO"]
    Outputs["SOL"] = PosInfo["SOL"]
    Outputs["NSVVIS"] = PosInfo["NSVVIS"]
    Outputs["NSV"] = PosInfo["NSV"]
    Outputs["HPE"] = PosInfo["HPE"]
    Outputs["VPE"] = PosInfo["VPE"]
    Outputs["EPE"] = PosInfo["EPE"]
    Outputs["NPE"] = PosInfo["NPE"]
    Outputs["UPE"] = PosInfo["UPE"]
    Outputs["HDOP"] = PosInfo["HDOP"]
    Outputs["VDOP"] = PosInfo["VDOP"]
    Outputs["PDOP"] = PosInfo["PDOP"]

    # Write line
    for i, result in enumerate(Outputs):
        fpos.write(((PosFmt[i] + " ") % Outputs[result]))

    fpos.write("\n")

# End of generatePosFile

def generatePerfFile(fperf, PerfInfo): 
    '''
    Purpose: generate output file with Performance results

    Parameters
    ==========
    fperf: file descriptor
            Descriptor for PERF output file
    PerfInfo: dict
            Dictionary containing Performance info 
    for the final computation
    
    Returns
    =======
    Nothing    
    '''

    # Prepare outputs
    Outputs = OrderedDict({})
    Outputs["RCVR"] = PerfInfo["RCVR"]
    Outputs["SAMPLES"] = PerfInfo["SAMPLES"]
    Outputs["SAMNOSOL"] = PerfInfo["SAMNOSOL"]
    Outputs["NSVMIN"] = PerfInfo["NSVMIN"]
    Outputs["NSVMAX"] = PerfInfo["NSVMAX"]
    Outputs["HPERMS"] = PerfInfo["HPERMS"]
    Outputs["VPERMS"] = PerfInfo["VPERMS"]
    Outputs["HPE95"] = PerfInfo["HPE95"]
    Outputs["VPE95"] = PerfInfo["VPE95"]
    Outputs["HPEMAX"] = PerfInfo["HPEMAX"]
    Outputs["VPEMAX"] = PerfInfo["VPEMAX"]
    Outputs["PDOPMAX"] = PerfInfo["PDOPMAX"]
    Outputs["HDOPMAX"] = PerfInfo["HDOPMAX"]
    Outputs["VDOPMAX"] = PerfInfo["VDOPMAX"]

    # Write line
    for i, result in enumerate(Outputs):
        fperf.write(((PerfFmt[i] + " ") % Outputs[result]))

    fperf.write("\n")

# End of generatePerfFile
