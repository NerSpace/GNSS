#!/usr/bin/env python

########################################################################
# SENTUS-V3/SRC/Perf.py:
# This is the Performance Module of SENTUS tool
#
#  Project:        SENTUS
#  File:           Perf.py
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
import numpy as np

def initializePerfInfo(Satname):

    """Initialize performance information for the satellite 
    receiver."""

        # Initialize output
    PerfInfo = {
        "RCVR": Satname,      # Rcvr Acronym with 4 characters
        "SAMPLES": 0.0,       # Number of Total Samples processed
        "SAMNOSOL": 0.0,      # Number of samples with no PVT Solution
        "NSVMIN": 40,         # Min nsats used in PVT Sol
        "NSVMAX": 0.0,        # Max nsats used in PVT Sol
        "HPERMS": 0.0,        # RMS of HPE
        "VPERMS": 0.0,        # RMS of VPE
        "HPE95": 0.0,         # 95% Percentile of HPE
        "VPE95": 0.0,         # 95% Percentile of VPE
        "HPEMAX": 0.0,        # Max reached HPE
        "VPEMAX": 0.0,        # Max reached VPE
        "PDOPMAX": 0.0,       # Max PDOP
        "HDOPMAX": 0.0,       # Max HDOP
        "VDOPMAX": 0.0,       # Max VDOP
        "HPE": [],
        "VPE": []
    } # End of PerfInfo

    return PerfInfo

def updatePerfEpoch(PosInfo, PerfInfo):

    """Update performance metrics based on position 
    information for the current epoch."""

    # Min and max values
    
    if PerfInfo["NSVMIN"] > PosInfo["NSV"] and PosInfo["NSV"]!=0:
        PerfInfo["NSVMIN"] = PosInfo["NSV"]

    if PerfInfo["NSVMAX"] < PosInfo["NSV"]:
        PerfInfo["NSVMAX"] = PosInfo["NSV"]

    if PerfInfo["PDOPMAX"] < PosInfo["PDOP"]:
        PerfInfo["PDOPMAX"] = PosInfo["PDOP"]

    if PerfInfo["HDOPMAX"] < PosInfo["HDOP"]:
        PerfInfo["HDOPMAX"] = PosInfo["HDOP"]

    if PerfInfo["VDOPMAX"] < PosInfo["VDOP"]:
        PerfInfo["VDOPMAX"] = PosInfo["VDOP"]

    PerfInfo["HPE"].append(PosInfo["HPE"])
    PerfInfo["VPE"].append(PosInfo["VPE"])

    return PerfInfo

def computeFinalPerf(PerfInfo):

    """Compute final performance metrics from accumulated data."""

    HPE = PerfInfo["HPE"] 
    VPE = PerfInfo["VPE"] 

    PerfInfo["HPERMS"] = 0
    PerfInfo["VPERMS"] = 0

    # HPE MAX
    PerfInfo["HPEMAX"] = max(HPE)
    # VPE MAX
    PerfInfo["VPEMAX"] = max(VPE)

    # HPE RMS
    squared_sum_HPE = sum([x**2 for x in HPE])
    mean_squared_HPE = squared_sum_HPE / len(HPE)
    PerfInfo["HPERMS"] = np.sqrt(mean_squared_HPE)
    # VPE RMS
    squared_sum_VPE = sum([x**2 for x in VPE])
    mean_squared_VPE = squared_sum_VPE / len(VPE)
    PerfInfo["VPERMS"] = np.sqrt(mean_squared_VPE)

    # HPE Percentile 95%
    # Step 1: Define bins (e.g., with a 0.001 meter resolution)
    bin_resolution_HPE = 0.001
    bins_HPE = np.arange(np.min(HPE), np.max(HPE), bin_resolution_HPE)

    # Step 2: Count the number of samples in each bin
    histogram_HPE, bin_edges_HPE = np.histogram(HPE, bins=bins_HPE)

    # Step 3: Compute the ratios (R_i = N_i / N_total)
    total_samples_HPE = len(HPE)
    ratios_HPE = histogram_HPE / total_samples_HPE

    # Step 4: Compute the cumulative sum of the ratios
    cumulative_sum_HPE = np.cumsum(ratios_HPE)

    # Step 5: Find the bin where the cumulative sum exceeds 0.95 (95th percentile)
    percentile_index_HPE = np.where(cumulative_sum_HPE > 0.95)[0][0]
    PerfInfo["HPE95"] = bin_edges_HPE[percentile_index_HPE + 1]  # Use the upper edge of the bin

    # VPE Percentile 95%
    # Step 1: Define bins (e.g., with a 0.001 meter resolution)
    bin_resolution_VPE = 0.001
    bins_VPE = np.arange(np.min(VPE), np.max(VPE), bin_resolution_VPE)

    # Step 2: Count the number of samples in each bin
    histogram_VPE, bin_edges_VPE = np.histogram(VPE, bins=bins_VPE)

    # Step 3: Compute the ratios (R_i = N_i / N_total)
    total_samples_VPE = len(VPE)
    ratios_VPE = histogram_VPE / total_samples_VPE

    # Step 4: Compute the cumulative sum of the ratios
    cumulative_sum_VPE = np.cumsum(ratios_VPE)

    # Step 5: Find the bin where the cumulative sum exceeds 0.95 (95th percentile)
    percentile_index_VPE = np.where(cumulative_sum_VPE > 0.95)[0][0]
    PerfInfo["VPE95"] = bin_edges_VPE[percentile_index_VPE + 1]  # Use the upper edge of the bin

    return PerfInfo




