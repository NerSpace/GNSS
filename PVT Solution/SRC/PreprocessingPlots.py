#!/usr/bin/env python

########################################################################
# PreprocessingPlots.py:
# This is the PreprocessingPlots Module of SENTUS tool
#
#  Project:        SENTUS
#  File:           PreprocessingPlots.py
#
#   Author: GNSS Academy
#   Copyright 2024 GNSS Academy / Nerea Sánchez
#
########################################################################


import sys, os
from pandas import unique
from pandas import read_csv
from InputOutput import PreproIdx
from InputOutput import REJECTION_CAUSE_DESC
sys.path.append(os.getcwd() + '/' + \
    os.path.dirname(sys.argv[0]) + '/' + 'COMMON')
from COMMON import GnssConstants
import numpy as np
from collections import OrderedDict
from COMMON.Plots import generatePlot
import matplotlib.pyplot as plt


def initPlot(PreproObsFile, Title, Label):
    PreproObsFileName = os.path.basename(PreproObsFile)
    PreproObsFileNameSplit = PreproObsFileName.split('_')
    Rcvr = PreproObsFileNameSplit[2]
    DatepDat = PreproObsFileNameSplit[3]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]
    
    PlotConf = {}
    PlotConf["xLabel"] = "Hour of Day %s" % Doy

    PlotConf["Title"] = "%s from %s on Year %s"\
        " DoY %s" % (Title, Rcvr, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/PPVE/figures/' + \
        '%s_%s_Y%sD%s.png' % (Label, Rcvr, Year, Doy)
    
    return PlotConf


# Function to convert 'G01', 'G02', etc. to 1, 2, etc.
def convert_satlabel_to_prn(value):
    return int(value[1:])


# Function to convert 'G01', 'G02', etc. to 'G'
def convert_satlabel_to_const(value):
    return value[0]

# Function to convert 'G01', 'G02', etc. to  1, 2, ... integers
def convert_prn_to_int(prn):
            return int(prn[1:])


# Plot Satellite Visibility

# First status = 0 and status = 1 are separated
def separateDataByStatus(PreproObsData):
    status_0_data = PreproObsData[PreproObsData[PreproIdx["STATUS"]] == 0]
    status_1_data = PreproObsData[PreproObsData[PreproIdx["STATUS"]] == 1]
    return status_0_data, status_1_data


def plotSatVisibility(PreproObsFile, PreproObsData):
 
    Title = "Satellite Visibility"
    Label = "SAT_VIS"
    PlotConf = initPlot(PreproObsFile, Title, Label)

    PlotConf["Type"] = "LinesStatus"
    PlotConf["FigSize"] = (10,6)

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["MarkerSize"] = None
    PlotConf["LineWidth"] = 0.01

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    # List with every PRN
    all_prns = [f'E{str(i).zfill(2)}' for i in range(1, 37)] + [f'G{str(i).zfill(2)}' for i in range(1, 36)]
    PlotConf["yTicksLabels"] = all_prns
    PlotConf["FontSize"] = 6
    PlotConf["yLabel"] = "PRN"

    # Adding all the PRNs not included in the PreproObsData
    missing_prns = [prn for prn in all_prns if prn not in PreproObsData[PreproIdx["PRN"]].unique()]

    status_0_data, status_1_data = separateDataByStatus(PreproObsData)

    # Adding to status_0 and status_1 the missing PRNs
    for prn in missing_prns:
        status_0_data = status_0_data.append({PreproIdx["PRN"]: prn}, ignore_index=True)
        status_1_data = status_1_data.append({PreproIdx["PRN"]: prn}, ignore_index=True)

    status_0_data = status_0_data.sort_values(by=[PreproIdx["PRN"]])
    status_1_data = status_1_data.sort_values(by=[PreproIdx["PRN"]])

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    PlotConf["Zorder"] = {}

    # Prepare data for status = 0
    Label = "status_0"
    PlotConf["xData"][Label]  = status_0_data[PreproIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = status_0_data[PreproIdx["PRN"]]
    PlotConf["zData"][Label] = [0.5] * len(status_0_data)
    PlotConf["Zorder"][Label] = 2 # These measurements overlap the other ones

    # Prepare data for status = 1
    Label = "status_1"
    PlotConf["xData"][Label] = status_1_data[PreproIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = status_1_data[PreproIdx["PRN"]]
    PlotConf["zData"][Label] = status_1_data[PreproIdx["ELEV"]]
    PlotConf["Zorder"][Label] = 1

    generatePlot(PlotConf)
    plt.close()



# Plot Number of Satellites
def plotNumSats(PreproObsFile, PreproObsData):
    Constel = ["GPS+GAL", "GPS", "GAL"]
    for CONSTEL in Constel:

        if CONSTEL == "GPS":
            Data_Constel = PreproObsData[PreproObsData[PreproIdx["PRN"]].str.startswith("G")]
        elif CONSTEL == "GAL":
            Data_Constel = PreproObsData[PreproObsData[PreproIdx["PRN"]].str.startswith("E")]
        elif CONSTEL == "GPS+GAL":
            Data_Constel = PreproObsData

        PlotConf = {}
        Title = f"Number of {CONSTEL} satellites"
        Label = f"N_{CONSTEL}_SAT"
        PlotConf = initPlot(PreproObsFile, Title, Label)
    
        PlotConf["Type"] = "Lines"
        PlotConf["FigSize"] = (12,7)

        PlotConf["xTicks"] = range(0, 25)
        PlotConf["xLim"] = [0, 24]
        PlotConf["yTicks"] = range(0, 21)
        PlotConf["yLim"] = [0, 20]
        PlotConf["yLabel"] = Title

        PlotConf["Grid"] = 1

        PlotConf["Marker"] = '-'
        PlotConf["MarkerSize"] = None
        PlotConf["LineWidth"] = 1.0

        PlotConf["LineColor"] = {
            "raw": "orange",
            "smoothed": "green"
        }
    
        PlotConf["LegendLabels"] = {
            "raw": "Raw",
            "smoothed": "Smoothed"
        }

        PlotConf["Legend"] = True

        # Calculate number of Sats
        NSATS_raw = Data_Constel.groupby(PreproIdx["SOD"])[PreproIdx["PRN"]].nunique()
        status_0_data, status_1_data = separateDataByStatus(Data_Constel)
        NSATS_smoothed = status_1_data.groupby(PreproIdx["SOD"])[PreproIdx["PRN"]].nunique()

        PlotConf["xData"] = {}
        PlotConf["yData"] = {}
        PlotConf["zData"] = {}

        # Prepare data for raw
        Label = "raw"
        PlotConf["xData"][Label]  = Data_Constel[PreproIdx["SOD"]].unique()/ GnssConstants.S_IN_H
        PlotConf["yData"][Label] = NSATS_raw

        # Prepare data for smoothed
        Label = "smoothed"
        PlotConf["xData"][Label]  = status_1_data[PreproIdx["SOD"]].unique() / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = NSATS_smoothed
    
        generatePlot(PlotConf)
        plt.close()

# Plot Code IF - Code IF Smoothed
    
def plotIFIFSmoothed(PreproObsFile, PreproObsData):
    
    Constel = ["GPS", "GAL"]
    for CONSTEL in Constel:

        if CONSTEL == "GPS":
            Data_Constel = PreproObsData[PreproObsData[PreproIdx["PRN"]].str.startswith("G")]
        elif CONSTEL == "GAL":
            Data_Constel = PreproObsData[PreproObsData[PreproIdx["PRN"]].str.startswith("E")]

        Title = f"{CONSTEL} Code IF - IF Smoothed"
        Label = f"{CONSTEL}_IF-IFSmoothed"
        PlotConf = initPlot(PreproObsFile, Title, Label)

        PlotConf["Type"] = "LinesStatus"
        PlotConf["FigSize"] = (10,6)

        PlotConf["xTicks"] = range(0, 25)
        PlotConf["xLim"] = [0, 24]
        
        PlotConf["yLabel"] = f"{Title} [m]"

        if CONSTEL == "GPS":
            PlotConf["yLim"] = [-3, 3]
        elif CONSTEL == "GAL":
            PlotConf["yLim"] = [-1.5, 1.5]

        PlotConf["Grid"] = 1

        PlotConf["Marker"] = '.'
        PlotConf["MarkerSize"] = 10
        PlotConf["LineWidth"] = 1

        PlotConf["ColorBar"] = "gnuplot"
        PlotConf["ColorBarLabel"] = "Elevation [deg]"
        PlotConf["ColorBarMin"] = 0.
        PlotConf["ColorBarMax"] = 90.

        status_0_data, status_1_data = separateDataByStatus(Data_Constel)
        PlotConf["xData"] = {}
        PlotConf["yData"] = {}
        PlotConf["zData"] = {}
        PlotConf["Zorder"] = {}

        # Prepare data for status = 0
        Label = "status_0"
        y_data_0 = status_0_data[PreproIdx["CODE_IF"]] - status_0_data[PreproIdx["SMOOTH_IF"]]
        PlotConf["xData"][Label]  = status_0_data[PreproIdx["SOD"]] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = y_data_0
        PlotConf["zData"][Label] = [0.5] * len(status_0_data)
        PlotConf["Zorder"][Label] = 1

        # Prepare data for status = 1
        Label = "status_1"
        y_data_1 = status_1_data[PreproIdx["CODE_IF"]] - status_1_data[PreproIdx["SMOOTH_IF"]]
        PlotConf["xData"][Label] = status_1_data[PreproIdx["SOD"]] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = y_data_1
        PlotConf["zData"][Label] = status_1_data[PreproIdx["ELEV"]]
        PlotConf["Zorder"][Label] = 2

        generatePlot(PlotConf)
        plt.close()

# Plot C/N0
def plotCN0(PreproObsFile, PreproObsData, PlotTitle, PlotLabel):
    Constel = ["GPS", "GAL"]
    for CONSTEL in Constel:

        if CONSTEL == "GPS":
            Data_Constel = PreproObsData[PreproObsData[PreproIdx["PRN"]].str.startswith("G")]
        elif CONSTEL == "GAL":
            Data_Constel = PreproObsData[PreproObsData[PreproIdx["PRN"]].str.startswith("E")]

        Title = f"{CONSTEL} {PlotTitle}"
        Label = f"{CONSTEL}_{PlotTitle}"
        PlotConf = initPlot(PreproObsFile, Title, Label)

        PlotConf["Type"] = "Lines"
        PlotConf["FigSize"] = (10,6)

        PlotConf["xTicks"] = range(0, 25)
        PlotConf["xLim"] = [0, 24]
        
        PlotConf["yLabel"] = f"{Title} [Db-Hz]"

        PlotConf["Grid"] = 1

        PlotConf["Marker"] = '.'
        PlotConf["MarkerSize"] = 10
        PlotConf["LineWidth"] = 0
        
        PlotConf["ColorBar"] = "gnuplot"
        PlotConf["ColorBarLabel"] = "Elevation [deg]"
        PlotConf["ColorBarMin"] = 0.
        PlotConf["ColorBarMax"] = 90.

        PlotConf["xData"] = {}
        PlotConf["yData"] = {}
        PlotConf["zData"] = {}


        Label = 0
        PlotConf["xData"][Label]  = Data_Constel[PreproIdx["SOD"]] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = Data_Constel[PreproIdx[PlotLabel]]
        PlotConf["zData"][Label] = Data_Constel[PreproIdx["ELEV"]]

        generatePlot(PlotConf)
        plt.close()



# Plot Rejection Flags
def plotRejectionFlags(PreproObsFile, PreproObsData):
    Constel = ["GPS", "GAL"]
    for CONSTEL in Constel:

        if CONSTEL == "GPS":
            Data_Constel = PreproObsData[PreproObsData[PreproIdx["PRN"]].str.startswith("G")]
        elif CONSTEL == "GAL":
            Data_Constel = PreproObsData[PreproObsData[PreproIdx["PRN"]].str.startswith("E")]

        Title = f"{CONSTEL} Rejection Flags"
        Label = f"{CONSTEL}_REJFLAGS"
        PlotConf = initPlot(PreproObsFile, Title, Label)
        PlotConf["CONSTEL"] = CONSTEL

        PlotConf["Type"] = "Lines"
        PlotConf["FigSize"] = (10,6)

        PlotConf["xTicks"] = range(0, 25)
        PlotConf["xLim"] = [0, 24]
        
        PlotConf["yTicks"] = range(1, 16)
        PlotConf["yLim"] = [1, 16]
        PlotConf["yTicksLabels"] = REJECTION_CAUSE_DESC
        PlotConf["yLabel"] = f"{Title}"

        PlotConf["Grid"] = 1

        PlotConf["Marker"] = 'o'
        PlotConf["MarkerSize"] = 25
        PlotConf["LineWidth"] = 1

        PlotConf["ColorBar"] = "gist_ncar"
        if CONSTEL == "GPS":
            PlotConf["ColorBarLabel"] = "GPS-PRN"
            PlotConf["ColorBarMax"] = 32
            PlotConf["DiscreteColorBar"] = 32
            
        elif CONSTEL == "GAL":
            PlotConf["ColorBarLabel"] = "GAL-PRN"
            PlotConf["ColorBarMax"] = 36
            PlotConf["DiscreteColorBar"] = 36

        PlotConf["zTicks"] = 1
        PlotConf["ColorBarMin"] = 2

        PlotConf["xData"] = {}
        PlotConf["yData"] = {}
        PlotConf["zData"] = {}

        # Prepare data for reject != 0
        nonvalid_data = Data_Constel[Data_Constel[PreproIdx["REJECT"]] != 0].reset_index(drop=True)
        intPRN_list = [convert_prn_to_int(prn) for prn in nonvalid_data [PreproIdx["PRN"]]]
        Label = 0
        PlotConf["xData"][Label]  = nonvalid_data[PreproIdx["SOD"]] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = nonvalid_data[PreproIdx["REJECT"]]
        PlotConf["zData"][Label] = intPRN_list

        # For adding the ticks names
        PlotConf["TickName"] = nonvalid_data[PreproIdx["PRN"]]

        generatePlot(PlotConf)
        plt.close()


# Plot Rates
def plotRates(PreproObsFile, PreproObsData, PlotTitle, PlotLabel):
    Constel = ["GPS", "GAL"]
    for CONSTEL in Constel:

        if CONSTEL == "GPS":
            Data_Constel = PreproObsData[PreproObsData[PreproIdx["PRN"]].str.startswith("G")]
        elif CONSTEL == "GAL":
            Data_Constel = PreproObsData[PreproObsData[PreproIdx["PRN"]].str.startswith("E")]

        Title = f"{CONSTEL} {PlotTitle}"
        Label = f"{CONSTEL}_{PlotTitle}"
        PlotConf = initPlot(PreproObsFile, Title, Label)

        PlotConf["Type"] = "Lines"
        PlotConf["FigSize"] = (10,6)

        PlotConf["xTicks"] = range(0, 25)
        PlotConf["xLim"] = [0, 24]
        
        if PlotLabel.endswith("STEP"):
            PlotConf["yLabel"] = f"{Title} [m/s²]"
        else:
            PlotConf["yLabel"] = f"{Title} [m/s]"
        
        PlotConf["Grid"] = 1

        PlotConf["Marker"] = '.'
        PlotConf["MarkerSize"] = 10
        PlotConf["LineWidth"] = 0

        PlotConf["ColorBar"] = "gnuplot"
        PlotConf["ColorBarLabel"] = "Elevation [deg]"
        PlotConf["ColorBarMin"] = 0.
        PlotConf["ColorBarMax"] = 90.

        PlotConf["xData"] = {}
        PlotConf["yData"] = {}
        PlotConf["zData"] = {}

        # Prepare data for valid = 1
        valid_data = Data_Constel[Data_Constel[PreproIdx["VALID"]] == 1]

        Label = 0
        PlotConf["xData"][Label]  = valid_data[PreproIdx["SOD"]] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = valid_data[PreproIdx[PlotLabel]]
        PlotConf["zData"][Label] = valid_data[PreproIdx["ELEV"]]

        generatePlot(PlotConf)

        plt.close()


def generatePreproPlots(PreproObsFile, Conf):
    '''
    Purpose: generate output plots regarding Preprocessing results

    Parameters
    ==========
    PreproObsFile: str
            Path to PREPRO OBS output file

    Returns
    =======
    Nothing    
    '''

    # Satellite Visibility
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"],PreproIdx["PRN"],PreproIdx["STATUS"],PreproIdx["ELEV"]])
    
    print('INFO: Plot Satellite Visibility Periods ...')

    # Configure plot and call plot generation function
    plotSatVisibility(PreproObsFile, PreproObsData)


    # Number of satellites
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"],PreproIdx["PRN"],PreproIdx["STATUS"]])
    
    print('INFO: Plot Number of Satellites ...')

    # Configure plot and call plot generation function
    plotNumSats(PreproObsFile, PreproObsData)

    # Code IF - Code IF Smoothed
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"],PreproIdx["REJECT"],PreproIdx["STATUS"],PreproIdx["ELEV"],\
    PreproIdx["PRN"],PreproIdx["CODE_IF"],PreproIdx["SMOOTH_IF"]])
    
    print('INFO: Plot Code IF - Code IF Smoothed ...')

    # Configure plot and call plot generation function
    plotIFIFSmoothed(PreproObsFile, PreproObsData)


    # C/N0
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"],PreproIdx["S1"],PreproIdx["S2"],\
        PreproIdx["ELEV"], PreproIdx["PRN"]])
    
    print('INFO: Plot C/N0...')

    # Configure plot and call plot generation function
    plotCN0(PreproObsFile, PreproObsData, 'CN0_F1', 'S1')

    # Configure plot and call plot generation function
    plotCN0(PreproObsFile, PreproObsData, 'CN0_F2', 'S2')


    # Rejection Flags
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"],PreproIdx["PRN"],PreproIdx["REJECT"]])
    
    print('INFO: Plot Rejection Flags ...')

    # Configure plot and call plot generation function
    plotRejectionFlags(PreproObsFile, PreproObsData)


    # Code Rate
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"],PreproIdx["VALID"],PreproIdx["CODE_RATE"],\
        PreproIdx["ELEV"], PreproIdx["PRN"]])
    
    print('INFO: Plot Code Rate ...')

    # Configure plot and call plot generation function
    plotRates(PreproObsFile, PreproObsData, 'Code Rate', 'CODE_RATE')
    plotRates(PreproObsFile, PreproObsData, 'Code Rate', 'CODE_RATE')

    # Phase Rate
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"],PreproIdx["VALID"],PreproIdx["PHASE_RATE"],\
        PreproIdx["ELEV"], PreproIdx["PRN"]])

    print('INFO: Plot Phase Rate ...')

    # Configure plot and call plot generation function
    plotRates(PreproObsFile, PreproObsData, 'Phase Rate', 'PHASE_RATE')
    plotRates(PreproObsFile, PreproObsData, 'Phase Rate', 'PHASE_RATE')


    # Code Rate Step
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"],PreproIdx["VALID"],PreproIdx["CODE_RATE_STEP"],\
        PreproIdx["ELEV"], PreproIdx["PRN"]])
    
    print('INFO: Plot Code Rate Step...')

    # Configure plot and call plot generation function
    plotRates(PreproObsFile, PreproObsData, 'Code Rate Step', 'CODE_RATE_STEP')
    plotRates(PreproObsFile, PreproObsData, 'Code Rate Step', 'CODE_RATE_STEP')


    # Phase Rate Step
    # ----------------------------------------------------------
    # Read the cols we need from PREPRO OBS file
    PreproObsData = read_csv(PreproObsFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PreproIdx["SOD"],PreproIdx["VALID"],PreproIdx["PHASE_RATE_STEP"],\
        PreproIdx["ELEV"], PreproIdx["PRN"]])
    
    print('INFO: Plot Phase Rate Step...')

    # Configure plot and call plot generation function
    plotRates(PreproObsFile, PreproObsData, 'Phase Rate Step', 'PHASE_RATE_STEP')
    plotRates(PreproObsFile, PreproObsData, 'Phase Rate Step', 'PHASE_RATE_STEP')

