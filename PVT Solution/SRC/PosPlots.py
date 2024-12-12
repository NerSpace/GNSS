#!/usr/bin/env python

########################################################################
# PosPlots.py:
# This is the Positioning Plots Module of SENTUS tool
#
#  Project:        SENTUS
#  File:           PosPlots.py
#
#   Author: Nerea Sánchez
#   Copyright 2024 GNSS Academy / Nerea Sánchez
#
########################################################################


import sys, os
from pandas import unique
from pandas import read_csv
from InputOutput import PosIdx
sys.path.append(os.getcwd() + '/' + \
    os.path.dirname(sys.argv[0]) + '/' + 'COMMON')
from COMMON import GnssConstants
import numpy as np
from collections import OrderedDict
from COMMON.Plots import generatePlot
import matplotlib.pyplot as plt
from COMMON.Coordinates import xyz2llh


def initPlot(PosFile, Title, Label):
    PosFileName = os.path.basename(PosFile)
    PosFileNameSplit = PosFileName.split('_')
    Rcvr = PosFileNameSplit[1]
    DatepDat = PosFileNameSplit[2]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]
    
    PlotConf = {}
    PlotConf["xLabel"] = "Hour of Day %s" % Doy

    PlotConf["Title"] = "%s from %s on Year %s"\
        " DoY %s" % (Title, Rcvr, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/PVTS/figures/' + \
        '%s_%s_Y%sD%s.png' % (Label, Rcvr, Year, Doy)
    
    return PlotConf

def initPlotTracks(PosFile, Title, Label):
    PosFileName = os.path.basename(PosFile)
    PosFileNameSplit = PosFileName.split('_')
    Rcvr = PosFileNameSplit[1]
    DatepDat = PosFileNameSplit[2]
    Date = DatepDat.split('.')[0]
    Year = Date[1:3]
    Doy = Date[4:]
    
    PlotConf = {}

    PlotConf["Title"] = "%s from %s on Year %s"\
        " DoY %s" % (Title, Rcvr, Year, Doy)

    PlotConf["Path"] = sys.argv[1] + '/OUT/PVTS/figures/' + \
        '%s_%s_Y%sD%s.png' % (Label, Rcvr, Year, Doy)
    
    return PlotConf

######### PLOT FUNCTIONS #########
##################################

# Number of Satellites
def plotNumSats(PosFile, PosData):

    PlotConf = {}
    Title = "Number of satellites"
    Label = "NUM_SATS"
    PlotConf = initPlot(PosFile, Title, Label)

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (12,7)

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["yTicks"] = range(0, 21)
    PlotConf["yLim"] = [0, 20]
    PlotConf["yLabel"] = Title

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '-'
    PlotConf["MarkerSize"] = 10
    PlotConf["LineWidth"] = 1.0

    PlotConf["LineColor"] = {
        "Raw": "orange",
        "Used": "green"
    }

    PlotConf["LegendLabels"] = {
        "Raw": "Raw",
        "Used": "Used"
    }

    PlotConf["Legend"] = True

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    # Prepare data for Raw
    Label = "Raw"
    PlotConf["xData"][Label]  = PosData[PosIdx["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PosData[PosIdx["NSVVIS"]]

    # Prepare data for Used
    Label = "Used"
    PlotConf["xData"][Label]  = PosData[PosIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PosData[PosIdx["NSV"]]

    generatePlot(PlotConf)
    plt.close()

# DOP Plot
def plotDOPs(PosFile, PosData):

    PlotConf = {}
    Title = "DOP"
    Label = "DOP"
    PlotConf = initPlot(PosFile, Title, Label)

    PlotConf["Type"] = "DualYAxis"
    PlotConf["FigSize"] = (12,7)

    # Left y-axis    
    PlotConf["yLabelLeft"] = Title
    # PlotConf["yLim"] = [0.5, 3]

    # Right y-axis
    PlotConf["yLabelRight"] = "Number of Satellites"
    PlotConf["yLimR"] = [0, 19]
    PlotConf["yTicksR"] = range(0, 20)    

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]


    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '-'
    PlotConf["MarkerSize"] = 10
    PlotConf["LineWidth"] = 1.0

    PlotConf["LineColor"] = {
        "PDOP": "blue",
        "VDOP": "green",
        "HDOP": "cyan",
        "NSATS": "orange"
    }

    PlotConf["LegendLabels"] = {
        "PDOP": "PDOP",
        "VDOP": "VDOP",
        "HDOP": "HDOP",
        "NSATS": "Num SV"
    }

    PlotConf["Legend"] = True

    PlotConf["xData"] = {}
    PlotConf["yDataLeft"] = {}
    PlotConf["yDataRight"] = {}

    # Prepare data for PDOP
    Label = "PDOP"
    PlotConf["xData"][Label]  = PosData[PosIdx["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yDataLeft"][Label] = PosData[PosIdx["PDOP"]]

    # Prepare data for VDOP
    Label = "VDOP"
    PlotConf["xData"][Label]  = PosData[PosIdx["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yDataLeft"][Label] = PosData[PosIdx["VDOP"]]

    # Prepare data for HDOP
    Label = "HDOP"
    PlotConf["xData"][Label]  = PosData[PosIdx["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yDataLeft"][Label] = PosData[PosIdx["HDOP"]]

    # Prepare data for Num Sats
    Label = "NSATS"
    PlotConf["xData"][Label]  = PosData[PosIdx["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yDataRight"][Label] = PosData[PosIdx["NSV"]]   


    generatePlot(PlotConf)
    plt.close()

# LEO Tracks [Satellite Longitude and Latitude (Elev)]
def plotSatTracks(PosFile, PosData):

        Title = "LEO Satellite Tracks"
        Label = "LEO_TRACKS"
        PlotConf = initPlotTracks(PosFile, Title, Label)    

       
        PlotConf["Type"] = "Lines"
        PlotConf["FigSize"] = (16.8,15.2)

        PlotConf["LonMin"] = -180
        PlotConf["LonMax"] = 180
        PlotConf["LatMin"] = -90
        PlotConf["LatMax"] = 90
        PlotConf["LonStep"] = 15
        PlotConf["LatStep"] = 10

        PlotConf["MarkerSize"] = 6

        PlotConf["yTicks"] = range(PlotConf["LatMin"],PlotConf["LatMax"]+1,10)
        PlotConf["yLim"] = [PlotConf["LatMin"], PlotConf["LatMax"]]

        PlotConf["xTicks"] = range(PlotConf["LonMin"],PlotConf["LonMax"]+1,15)
        PlotConf["xLim"] = [PlotConf["LonMin"], PlotConf["LonMax"]]

        PlotConf["Grid"] = True

        PlotConf["Map"] = True

        PlotConf["Marker"] = '.'
        PlotConf["LineWidth"] = 0.01

        PlotConf["ColorBar"] = "gnuplot"
        PlotConf["ColorBarLabel"] = "Second of Day"
        PlotConf["ColorBarMin"] = 0.
        PlotConf["ColorBarMax"] = 86400.

        PlotConf["xData"] = {}
        PlotConf["yData"] = {}
        PlotConf["zData"] = {}
        Label = 0
        PlotConf["xData"][Label] = PosData[PosIdx["LONEST"]]
        PlotConf["yData"][Label] = PosData[PosIdx["LATEST"]]
        PlotConf["zData"][Label] = PosData[PosIdx["SOD"]]


        generatePlot(PlotConf)
        plt.close()

# Receiver Clock
def plotRcvrClk(PosFile, PosData):

    PlotConf = {}
    Title = "Estimated Receiver Clock wrt GPST"
    Label = "RCVR_CLK"
    PlotConf = initPlot(PosFile, Title, Label)

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (12,7)

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["yLabel"] = f"{Title} [m]"

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["MarkerSize"] = 6
    PlotConf["LineWidth"] = 1.0

    PlotConf["LineColor"] = {
        "CLKEST": "green"
    }

    PlotConf["LegendLabels"] = {
        "CLKEST": "CLKEST"
    }

    PlotConf["Legend"] = True

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}

    # Prepare data for CLKEST
    Label = "CLKEST"
    PlotConf["xData"][Label]  = PosData[PosIdx["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PosData[PosIdx["CLKEST"]]

    generatePlot(PlotConf)
    plt.close()

# GGTO

def plotGGTO(PosFile, PosData):

    PlotConf = {}
    Title = "Estimated GGTO"
    Label = "GGTO"
    PlotConf = initPlot(PosFile, Title, Label)

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (12,7)

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["yLabel"] = f"{Title} [ns]"

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["MarkerSize"] = 6
    PlotConf["LineWidth"] = 1.0

    PlotConf["LineColor"] = {"GGTO": "red"}

    PlotConf["Legend"] = False

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}

    Label = "GGTO"
    PlotConf["xData"][Label]  = PosData[PosIdx["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PosData[PosIdx["GGTO"]]

    generatePlot(PlotConf)
    plt.close()

# ENU PE

def plotENUPE(PosFile, PosData):

    PlotConf = {}
    Title = "East North Up Position Errors"
    Label = "ENU_PE"
    PlotConf = initPlot(PosFile, Title, Label)

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (12,7)

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["yLabel"] = f"{Title} [m]"

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["MarkerSize"] = 6
    PlotConf["LineWidth"] = 1.0

    PlotConf["LineColor"] = {
        "EPE": "orange",
        "NPE": "red",
        "UPE": "green"
    }

    PlotConf["LegendLabels"] = {
        "EPE": "EPE",
        "NPE": "NPE",
        "UPE": "UPE"
    }

    PlotConf["Legend"] = True

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}

    # Prepare data for EPE
    Label = "EPE"
    PlotConf["xData"][Label]  = PosData[PosIdx["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PosData[PosIdx["EPE"]]

    # Prepare data for NPE
    Label = "NPE"
    PlotConf["xData"][Label]  = PosData[PosIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PosData[PosIdx["NPE"]]

    # Prepare data for UPE
    Label = "UPE"
    PlotConf["xData"][Label]  = PosData[PosIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PosData[PosIdx["UPE"]]

    generatePlot(PlotConf)
    plt.close()

# EPE vs NPE [HDOP]
    
def plotHPEScatter(PosFile, PosData):

    Title = f"Horizontal Position Error vs DOP"
    Label = f"EPE_VS_NPE"
    PlotConf = initPlot(PosFile, Title, Label)

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    
    PlotConf["yLabel"] = "NPE [m]"
    PlotConf["xLabel"] = "EPE [m]"

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["MarkerSize"] = 6
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "HDOP"
    PlotConf["ColorBarMin"] = 0.6
    PlotConf["ColorBarMax"] = 1.6

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    # Prepare data for status = 0
    Label = "status_0"
    PlotConf["xData"][Label]  = PosData[PosIdx["EPE"]]
    PlotConf["yData"][Label] = PosData[PosIdx["NPE"]]
    PlotConf["zData"][Label] = PosData[PosIdx["HDOP"]]

    generatePlot(PlotConf)
    plt.close()

# H/VPE

def plotHVPE(PosFile, PosData):

    PlotConf = {}
    Title = "Horizontal and Vertical Position Errors"
    Label = "HV_PE"
    PlotConf = initPlot(PosFile, Title, Label)

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (12,7)

    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["yLabel"] = Title

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["MarkerSize"] = 6
    PlotConf["LineWidth"] = 1.0

    PlotConf["LineColor"] = {
        "VPE": "green",
        "HPE": "red"
    }

    PlotConf["LegendLabels"] = {
        "VPE": "VPE",
        "HPE": "HPE"
    }

    PlotConf["Legend"] = True

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}

    # Prepare data for VPE
    Label = "VPE"
    PlotConf["xData"][Label]  = PosData[PosIdx["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PosData[PosIdx["VPE"]]

    # Prepare data for HPE
    Label = "HPE"
    PlotConf["xData"][Label]  = PosData[PosIdx["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PosData[PosIdx["HPE"]]

    generatePlot(PlotConf)
    plt.close()


######### GENERATE Pos FUNCT #########
#######################################

def generatePosPlots(PosFile, Conf):

    '''
    Purpose: generate output plots regarding Preprocessing results
    
    Parameters
    ==========
    PosFile: str
            Path to PREPRO OBS output file

    Returns
    =======
    Nothing    
    '''

    # Get which constellation is being Posected
    if Conf["NAV_SOLUTION"] == "GPS":
        Constel = ["GPS"]
    elif Conf["NAV_SOLUTION"] == "GAL":
        Constel = ["GAL"]
    elif Conf["NAV_SOLUTION"] == "GPSGAL":
        Constel = ["GPS", "GAL"]

    # Number of Satellites
    # ----------------------------------------------------------
    # Read the cols we need from Pos file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PosIdx["NSV"], PosIdx["NSVVIS"], PosIdx["SOD"]])

    print('INFO: Plot Number of Satellites ...')

    # Configure plot and call plot generation function
    plotNumSats(PosFile, PosData)

    # DOP Plot
    # ----------------------------------------------------------
    # Read the cols we need from Pos file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PosIdx["NSV"], PosIdx["PDOP"], PosIdx["VDOP"], PosIdx["HDOP"], PosIdx["SOD"]])

    print('INFO: Plot DOP ...')

    # Configure plot and call plot generation function
    plotDOPs(PosFile, PosData)
    
    # LEO Tracks
    # ----------------------------------------------------------
    # Read the cols we need from Pos file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PosIdx["SOD"], PosIdx["LONEST"], PosIdx["LATEST"], PosIdx["ALTEST"]])
    
    print('INFO: Plot LEO Tracks ...')

    # Configure plot and call plot generation function
    plotSatTracks(PosFile, PosData)

    # Receiver Clock
    # ----------------------------------------------------------
    # Read the cols we need from Pos file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PosIdx["SOD"], PosIdx["CLKEST"]])
    
    print('INFO: Plot Receiver Clock ...')

    # Configure plot and call plot generation function
    plotRcvrClk(PosFile, PosData)

    # GGTO
    # ----------------------------------------------------------
    # Read the cols we need from Pos file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PosIdx["SOD"], PosIdx["GGTO"]])
    
    print('INFO: Plot GGTO ...')

    # Configure plot and call plot generation function
    plotGGTO(PosFile, PosData)

    # ENU PE
    # ----------------------------------------------------------
    # Read the cols we need from Pos file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PosIdx["SOD"], PosIdx["EPE"], PosIdx["NPE"], PosIdx["UPE"]])
    
    print('INFO: Plot ENU PE ...')

    # Configure plot and call plot generation function
    plotENUPE(PosFile, PosData)

    # EPE vs NPE [HDOP]
    # ----------------------------------------------------------
    # Read the cols we need from Pos file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PosIdx["EPE"], PosIdx["NPE"], PosIdx["HDOP"]])
    
    print('INFO: Plot Horizontal PE vs DOP ...')

    # Configure plot and call plot generation function
    plotHPEScatter(PosFile, PosData)

    # HVPE
    # ----------------------------------------------------------
    # Read the cols we need from Pos file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[PosIdx["SOD"], PosIdx["VPE"], PosIdx["HPE"]])
    
    print('INFO: Plot H/VPE ...')

    # Configure plot and call plot generation function
    plotHVPE(PosFile, PosData)






    


