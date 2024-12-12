#!/usr/bin/env python


########################################################################
# SwRcvr.py:
# Main Script of GNSS Academy Software Defined Receiver
#
#  Project:        sw-rcvr-python
#  File:           SwRcvr.py
#
#   Author: Nerea Sánchez / GNSS Academy
#   Copyright 2024 Nerea Sánchez / GNSS Academy
#
########################################################################

# Dependencies:
# * python 3
# * numpy
# * scipy
# * matplotlib

''' Uncomment depending on which constellation you want to acquire/track'''
# from SettingsGps import Settings as Settings
from SettingsGal import Settings as Settings

import os
import numpy as np
from Acquisition import acquisitionGpsL1C, acquisitionGalE1B
from Tracking import tracking, loadChannels
from Measurements import buildMeas

# Init settings
settings = Settings()

# Read input file
npyFile = settings.inputf + '.' + str(settings.msToProcess)
npyFilenpy = settings.inputf + '.' + str(settings.msToProcess) + '.npy'
if(os.path.isfile(npyFilenpy)):
    inputSignal = np.load(npyFilenpy)
else:
    inputSignal = np.loadtxt(settings.inputf, 
    dtype=int, 
    max_rows=round(settings.msToProcess * 1e-3 * settings.samplingFreq))

    np.save(npyFile, inputSignal)

# Run Acquisition
if "gpsl1c" in settings.signal:
    acqResultsGpsL1C = acquisitionGpsL1C(settings, inputSignal)
if "gale1b" in settings.signal:
    acqResultsGalE1B = acquisitionGalE1B(settings, inputSignal)

# Create Channels class
class Channels:
    PRN                = np.zeros(settings.numberOfChannels, dtype=int)
    acquiredFreq       = np.zeros(settings.numberOfChannels)
    acquiredDoppler    = np.zeros(settings.numberOfChannels)
    acquiredDelay      = np.zeros(settings.numberOfChannels)
    acquiredDelayChips = np.zeros(settings.numberOfChannels)
    SNR                = np.zeros(settings.numberOfChannels)
    status             = np.zeros(settings.numberOfChannels, dtype=int)

channelsGps = Channels()
channelsGal = Channels()

print("Done!")