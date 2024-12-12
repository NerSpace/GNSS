#!/usr/bin/env python

########################################################################
# PreprocessingPlots.py:
# This is the PreprocessingPlots Module of SENTUS tool
#
#  Project:        SENTUS
#  File:           PreprocessingPlots.py
#
#   Author: Nerea Sánchez
#   Copyright 2024 GNSS Academy / Nerea Sánchez
#
########################################################################


import sys, os
from pandas import unique
from pandas import read_csv
from InputOutput import CorrIdx
sys.path.append(os.getcwd() + '/' + \
    os.path.dirname(sys.argv[0]) + '/' + 'COMMON')
from COMMON import GnssConstants
import numpy as np
from collections import OrderedDict
from COMMON.Plots import generatePlot
import matplotlib.pyplot as plt
from COMMON.Coordinates import xyz2llh


def initPlot(CorrFile, Title, Label):
    CorrFileName = os.path.basename(CorrFile)
    CorrFileNameSplit = CorrFileName.split('_')
    Rcvr = CorrFileNameSplit[1]
    DatepDat = CorrFileNameSplit[2]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]
    
    PlotConf = {}
    PlotConf["xLabel"] = "Hour of Day %s" % Doy

    PlotConf["Title"] = "%s from %s on Year %s"\
        " DoY %s" % (Title, Rcvr, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/CORR/figures/' + \
        '%s_%s_Y%sD%s.png' % (Label, Rcvr, Year, Doy)
    
    return PlotConf

def initPlotTracks(CorrFile, Title, Label):
    CorrFileName = os.path.basename(CorrFile)
    CorrFileNameSplit = CorrFileName.split('_')
    Rcvr = CorrFileNameSplit[1]
    DatepDat = CorrFileNameSplit[2]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]
    
    PlotConf = {}

    PlotConf["Title"] = "%s from %s on Year %s"\
        " DoY %s" % (Title, Rcvr, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/CORR/figures/' + \
        '%s_%s_Y%sD%s.png' % (Label, Rcvr, Year, Doy)
    
    return PlotConf

######### PLOT FUNCTIONS #########
##################################

# Satellite Tracks [Satellite Longitude and Latitude (Elev)]
def plotSatTracks(CorrFile, CorrData, Constel):

    # Get values with flag == 1
    CorrData = CorrData[CorrData[CorrIdx["FLAG"]]==1]
    CorrData = CorrData.reset_index(drop=True)

    for CONSTEL in Constel:

        if CONSTEL == "GPS":
            Data_Constel = CorrData[CorrData[CorrIdx["CONST"]]=="G"]
            Data_Constel = Data_Constel.reset_index(drop=True)
        elif CONSTEL == "GAL":
            Data_Constel = CorrData[CorrData[CorrIdx["CONST"]]=="E"]
            Data_Constel = Data_Constel.reset_index(drop=True)

        Title = f"{CONSTEL} Satellite Tracks"
        Label = f"{CONSTEL}_SATTRACKS"
        PlotConf = initPlotTracks(CorrFile, Title, Label)    

       
        PlotConf["Type"] = "Lines"
        PlotConf["FigSize"] = (16.8,15.2)

        PlotConf["LonMin"] = -180
        PlotConf["LonMax"] = 180
        PlotConf["LatMin"] = -60
        PlotConf["LatMax"] = 60
        PlotConf["LonStep"] = 15
        PlotConf["LatStep"] = 10

        PlotConf["MarkerSize"] = 10

        PlotConf["yTicks"] = range(PlotConf["LatMin"],PlotConf["LatMax"]+1,10)
        PlotConf["yLim"] = [PlotConf["LatMin"], PlotConf["LatMax"]]

        PlotConf["xTicks"] = range(PlotConf["LonMin"],PlotConf["LonMax"]+1,15)
        PlotConf["xLim"] = [PlotConf["LonMin"], PlotConf["LonMax"]]

        PlotConf["Grid"] = True

        PlotConf["Map"] = True

        PlotConf["Marker"] = '.'
        PlotConf["LineWidth"] = 0.01

        PlotConf["ColorBar"] = "gnuplot"
        PlotConf["ColorBarLabel"] = "Elevation [deg]"
        PlotConf["ColorBarMin"] = 0.
        PlotConf["ColorBarMax"] = 90.

        DataLen = len(Data_Constel[CorrIdx["SAT-X"]])
        Longitude = np.zeros(DataLen)
        Latitude = np.zeros(DataLen)
        
        for index in range(DataLen):
            x = Data_Constel[CorrIdx["SAT-X"]][index]
            y = Data_Constel[CorrIdx["SAT-Y"]][index]
            z = Data_Constel[CorrIdx["SAT-Z"]][index]
            Longitude[index], Latitude[index], h = xyz2llh(x, y, z)

        PlotConf["xData"] = {}
        PlotConf["yData"] = {}
        PlotConf["zData"] = {}
        Label = 0
        PlotConf["xData"][Label] = Longitude
        PlotConf["yData"][Label] = Latitude
        PlotConf["zData"][Label] = Data_Constel[CorrIdx["ELEV"]]


        generatePlot(PlotConf)
        plt.close()

# Plot Flight Time [ELEV]
    
def plotFlightTime(CorrFile, CorrData, Constel):

    # Get values with flag == 1
    CorrData = CorrData[CorrData[CorrIdx["FLAG"]]==1]
    CorrData = CorrData.reset_index(drop=True)
    
    for CONSTEL in Constel:

        if CONSTEL == "GPS":
            Data_Constel = CorrData[CorrData[CorrIdx["CONST"]]=="G"]
        elif CONSTEL == "GAL":
            Data_Constel = CorrData[CorrData[CorrIdx["CONST"]]=="E"]

        Title = f"{CONSTEL} Flight Time"
        Label = f"{CONSTEL}_TOF"
        PlotConf = initPlot(CorrFile, Title, Label)

        PlotConf["Type"] = "Lines"
        PlotConf["FigSize"] = (10,6)

        PlotConf["xTicks"] = range(0, 25)
        PlotConf["xLim"] = [0, 24]
        
        PlotConf["yLabel"] = "Flight Time [ms]"

        PlotConf["Grid"] = 1

        PlotConf["Marker"] = '.'
        PlotConf["MarkerSize"] = 6
        PlotConf["LineWidth"] = 1

        PlotConf["ColorBar"] = "gnuplot"
        PlotConf["ColorBarLabel"] = "Elevation [deg]"
        PlotConf["ColorBarMin"] = 0.
        PlotConf["ColorBarMax"] = 90.

        PlotConf["xData"] = {}
        PlotConf["yData"] = {}
        PlotConf["zData"] = {}

        # Prepare data for status = 0
        Label = "status_0"
        PlotConf["xData"][Label]  = Data_Constel[CorrIdx["SOD"]] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = Data_Constel[CorrIdx["FLIGHT-TIME"]]
        PlotConf["zData"][Label] = Data_Constel[CorrIdx["ELEV"]]

        generatePlot(PlotConf)
        plt.close()

# Plot DTR [ELEV]
    
def plotDtr(CorrFile, CorrData, Constel):
    
    # Get values with flag == 1
    CorrData = CorrData[CorrData[CorrIdx["FLAG"]]==1]
    CorrData = CorrData.reset_index(drop=True)

    for CONSTEL in Constel:

        if CONSTEL == "GPS":
            Data_Constel = CorrData[CorrData[CorrIdx["CONST"]]=="G"]
        elif CONSTEL == "GAL":
            Data_Constel = CorrData[CorrData[CorrIdx["CONST"]]=="E"]

        Title = f"{CONSTEL} Relativistic Correction (Dtr)"
        Label = f"{CONSTEL}_DTR"
        PlotConf = initPlot(CorrFile, Title, Label)

        PlotConf["Type"] = "Lines"
        PlotConf["FigSize"] = (10,6)

        PlotConf["xTicks"] = range(0, 25)
        PlotConf["xLim"] = [0, 24]
        
        PlotConf["yLabel"] = "Relativistic Correction (Dtr) [m]"

        # if CONSTEL == "GPS":
        #     PlotConf["yLim"] = [-3, 3]
        # elif CONSTEL == "GAL":
        #     PlotConf["yLim"] = [-1.5, 1.5]

        PlotConf["Grid"] = 1

        PlotConf["Marker"] = '.'
        PlotConf["MarkerSize"] = 10
        PlotConf["LineWidth"] = 1
        
        PlotConf["ColorBar"] = "gnuplot"
        PlotConf["ColorBarLabel"] = "Elevation [deg]"
        PlotConf["ColorBarMin"] = 0.
        PlotConf["ColorBarMax"] = 90.

        PlotConf["xData"] = {}
        PlotConf["yData"] = {}
        PlotConf["zData"] = {}

        Label = " "
        PlotConf["xData"][Label]  = Data_Constel[CorrIdx["SOD"]] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = Data_Constel[CorrIdx["DTR"]]
        PlotConf["zData"][Label] = Data_Constel[CorrIdx["ELEV"]]

        generatePlot(PlotConf)
        plt.close()

# Plot Code Residual [PRN]
    
def plotCodeRes(CorrFile, CorrData, Constel):
    
    # Get values with flag == 1
    CorrData = CorrData[CorrData[CorrIdx["FLAG"]]==1]
    CorrData = CorrData.reset_index(drop=True)
    
    for CONSTEL in Constel:

        if CONSTEL == "GPS":
            Data_Constel = CorrData[CorrData[CorrIdx["CONST"]]=="G"]
        elif CONSTEL == "GAL":
            Data_Constel = CorrData[CorrData[CorrIdx["CONST"]]=="E"]

        Title = f"{CONSTEL} Code Residual"
        Label = f"{CONSTEL}_CODRES"
        PlotConf = initPlot(CorrFile, Title, Label)

        PlotConf["Type"] = "Lines"
        PlotConf["FigSize"] = (10,6)

        PlotConf["xTicks"] = range(0, 25)
        PlotConf["xLim"] = [0, 24]
        
        PlotConf["yLabel"] = "Code Residuals [m]"

        PlotConf["Grid"] = 1

        PlotConf["Marker"] = '.'
        PlotConf["MarkerSize"] = 2
        PlotConf["LineWidth"] = 1

        PlotConf["ColorBar"] = "gnuplot"
        PlotConf["ColorBarLabel"] = "PRN"
        PlotConf["ColorBarMin"] = 1.
        PlotConf["ColorBarMax"] = GnssConstants.MAX_NUM_SATS_CONSTEL
        PlotConf["zTicks"] = [1., GnssConstants.MAX_NUM_SATS_CONSTEL]

        PlotConf["xData"] = {}
        PlotConf["yData"] = {}
        PlotConf["zData"] = {}

        Label = " "
        PlotConf["xData"][Label]  = Data_Constel[CorrIdx["SOD"]] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = Data_Constel[CorrIdx["CODE-RES"]]
        PlotConf["zData"][Label] = Data_Constel[CorrIdx["PRN"]]

        generatePlot(PlotConf)
        plt.close()


# Plot Phase Residual [PRN]
    
def plotPhaseRes(CorrFile, CorrData, Constel):
    
    # Get values with flag == 1
    CorrData = CorrData[CorrData[CorrIdx["FLAG"]]==1]
    CorrData = CorrData.reset_index(drop=True)

    for CONSTEL in Constel:

        if CONSTEL == "GPS":
            Data_Constel = CorrData[CorrData[CorrIdx["CONST"]]=="G"]
        elif CONSTEL == "GAL":
            Data_Constel = CorrData[CorrData[CorrIdx["CONST"]]=="E"]

        Title = f"{CONSTEL} Phase Residual"
        Label = f"{CONSTEL}_PHARES"
        PlotConf = initPlot(CorrFile, Title, Label)

        PlotConf["Type"] = "Lines"
        PlotConf["FigSize"] = (10,6)

        PlotConf["xTicks"] = range(0, 25)
        PlotConf["xLim"] = [0, 24]
        
        PlotConf["yLabel"] = "Phase Residuals [m]"

        PlotConf["Grid"] = 1

        PlotConf["Marker"] = '.'
        PlotConf["MarkerSize"] = 2
        PlotConf["LineWidth"] = 1

        PlotConf["ColorBar"] = "gnuplot"
        PlotConf["ColorBarLabel"] = "PRN"
        PlotConf["ColorBarMin"] = 1.
        PlotConf["ColorBarMax"] = GnssConstants.MAX_NUM_SATS_CONSTEL
        PlotConf["zTicks"] = [1., GnssConstants.MAX_NUM_SATS_CONSTEL]

        PlotConf["xData"] = {}
        PlotConf["yData"] = {}
        PlotConf["zData"] = {}

        Label = " "
        PlotConf["xData"][Label]  = Data_Constel[CorrIdx["SOD"]] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = Data_Constel[CorrIdx["PHASE-RES"]]
        PlotConf["zData"][Label] = Data_Constel[CorrIdx["PRN"]]

        generatePlot(PlotConf)
        plt.close()

# Receiver Clock
def plotRcvrClk(CorrFile, CorrData, Constel):

    # Get values with flag == 1
    CorrData = CorrData[CorrData[CorrIdx["FLAG"]]==1]
    CorrData = CorrData.reset_index(drop=True)

    for CONSTEL in Constel:

        if CONSTEL == "GPS":
            Data_Constel_GPS = CorrData[CorrData[CorrIdx["CONST"]]=="G"]
        elif CONSTEL == "GAL":
            Data_Constel_GAL = CorrData[CorrData[CorrIdx["CONST"]]=="E"]

    PlotConf = {}
    Title = "Receiver Clock Estimation [m]"
    Label = "RCVRCLK"
    PlotConf = initPlot(CorrFile, Title, Label)

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (12,7)

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["yTicks"] = range(0, 21)
    PlotConf["yLabel"] = Title

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '-'
    PlotConf["MarkerSize"] = 10
    PlotConf["LineWidth"] = 1.0

    PlotConf["LineColor"] = {
        "GPS": "red",
        "GAL": "blue"
    }

    PlotConf["LegendLabels"] = {
        "GPS": "GPS",
        "GAL": "GAL"
    }

    PlotConf["Legend"] = True

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    # Prepare data for raw
    Label = "GPS"
    PlotConf["xData"][Label]  = Data_Constel_GPS[CorrIdx["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = Data_Constel_GPS[CorrIdx["RCVR-CLK"]]

    # Prepare data for smoothed
    Label = "GAL"
    PlotConf["xData"][Label]  = Data_Constel_GAL[CorrIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = Data_Constel_GAL[CorrIdx["RCVR-CLK"]]

    generatePlot(PlotConf)
    plt.close()

def plotRcvrClkGPS(CorrFile, CorrData, Constel):

    # Get values with flag == 1
    CorrData = CorrData[CorrData[CorrIdx["FLAG"]]==1]
    CorrData = CorrData.reset_index(drop=True)

    Data_Constel_GPS = CorrData[CorrData[CorrIdx["CONST"]]=="G"]

    PlotConf = {}
    Title = "Receiver Clock Estimation [m]"
    Label = "RCVRCLK"
    PlotConf = initPlot(CorrFile, Title, Label)

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (12,7)

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["yTicks"] = range(0, 21)
    PlotConf["yLabel"] = Title

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '-'
    PlotConf["MarkerSize"] = 10
    PlotConf["LineWidth"] = 1.0

    PlotConf["LineColor"] = {
        "GPS": "red"
    }

    PlotConf["LegendLabels"] = {
        "GPS": "GPS"
    }

    PlotConf["Legend"] = True

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    # Prepare data for raw
    Label = "GPS"
    PlotConf["xData"][Label]  = Data_Constel_GPS[CorrIdx["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = Data_Constel_GPS[CorrIdx["RCVR-CLK"]]

    generatePlot(PlotConf)
    plt.close()

def plotRcvrClkGAL(CorrFile, CorrData, Constel):

    # Get values with flag == 1
    CorrData = CorrData[CorrData[CorrIdx["FLAG"]]==1]
    CorrData = CorrData.reset_index(drop=True)

    Data_Constel_GAL = CorrData[CorrData[CorrIdx["CONST"]]=="E"]

    PlotConf = {}
    Title = "Receiver Clock Estimation [m]"
    Label = "RCVRCLK"
    PlotConf = initPlot(CorrFile, Title, Label)

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (12,7)

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["yTicks"] = range(0, 21)
    PlotConf["yLabel"] = Title

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '-'
    PlotConf["MarkerSize"] = 10
    PlotConf["LineWidth"] = 1.0

    PlotConf["LineColor"] = {
        "GAL": "blue"
    }

    PlotConf["LegendLabels"] = {
        "GAL": "GAL"
    }

    PlotConf["Legend"] = True

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    # Prepare data for smoothed
    Label = "GAL"
    PlotConf["xData"][Label]  = Data_Constel_GAL[CorrIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = Data_Constel_GAL[CorrIdx["RCVR-CLK"]]

    generatePlot(PlotConf)
    plt.close()

######### GENERATE CORR FUNCT #########
#######################################

def generateCorrPlots(CorrFile, Conf):

    '''
    
    Purpose: generate output plots regarding Preprocessing results
    
    Parameters
    ==========
    CorrFile: str
            Path to PREPRO OBS output file

    Returns
    =======
    Nothing

    '''

    # Get which constellation is being corrected
    if Conf["NAV_SOLUTION"] == "GPS":
        Constel = ["GPS"]
    elif Conf["NAV_SOLUTION"] == "GAL":
        Constel = ["GAL"]
    elif Conf["NAV_SOLUTION"] == "GPSGAL":
        Constel = ["GPS", "GAL"]


    # Satellite Tracks
    # ----------------------------------------------------------
    # Read the cols we need from CORR file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["ELEV"], CorrIdx["SAT-X"], CorrIdx["SAT-Y"], CorrIdx["SAT-Z"], CorrIdx["CONST"],  CorrIdx["FLAG"]])
    
    print('INFO: Plot Satellite Tracks ...')

    # Configure plot and call plot generation function
    plotSatTracks(CorrFile, CorrData, Constel)

    # Flight Time
    # ----------------------------------------------------------
    # Read the cols we need from CORR file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["ELEV"], CorrIdx["SOD"], CorrIdx["FLIGHT-TIME"] , CorrIdx["CONST"],  CorrIdx["FLAG"]])
    
    print('INFO: Plot Flight Time ...')

    # Configure plot and call plot generation function
    plotFlightTime(CorrFile, CorrData, Constel)

    # DTR
    # ----------------------------------------------------------
    # Read the cols we need from CORR file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["ELEV"], CorrIdx["SOD"], CorrIdx["DTR"] , CorrIdx["CONST"],  CorrIdx["FLAG"]])
    
    print('INFO: Plot DTR ...')

    # Configure plot and call plot generation function
    plotDtr(CorrFile, CorrData, Constel)

    # Code Residual
    # ----------------------------------------------------------
    # Read the cols we need from CORR file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["SOD"], CorrIdx["CODE-RES"] , CorrIdx["CONST"], CorrIdx["PRN"],  CorrIdx["FLAG"]])
    
    print('INFO: Plot Code Residual ...')

    # Configure plot and call plot generation function
    plotCodeRes(CorrFile, CorrData, Constel)

    # Phase Residual
    # ----------------------------------------------------------
    # Read the cols we need from CORR file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["SOD"], CorrIdx["PHASE-RES"] , CorrIdx["CONST"], CorrIdx["PRN"],  CorrIdx["FLAG"]])
    
    print('INFO: Plot Phase Residual ...')

    # Configure plot and call plot generation function
    plotPhaseRes(CorrFile, CorrData, Constel)
    
    # Receiver Clock Estimation
    # ----------------------------------------------------------
    # Read the cols we need from CORR file
    CorrData = read_csv(CorrFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[CorrIdx["SOD"], CorrIdx["RCVR-CLK"] , CorrIdx["CONST"],  CorrIdx["FLAG"]])
    
    print('INFO: Plot Receiver Clock Estimation ...')

    # Configure plot and call plot generation function
    if Conf["NAV_SOLUTION"] == "GPS":
        plotRcvrClkGPS(CorrFile, CorrData, Constel)
    elif Conf["NAV_SOLUTION"] == "GAL":
        plotRcvrClkGAL(CorrFile, CorrData, Constel)
    elif Conf["NAV_SOLUTION"] == "GPSGAL":
        plotRcvrClk(CorrFile, CorrData, Constel)
    


