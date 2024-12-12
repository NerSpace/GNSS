

########################################################################
# Measurements.py:
# Measurements module
#
#  Project:        SW-RCVR-PYTHON
#  File:           Measurements.py
#
#   Author: Nerea Sánchez / GNSS Academy
#   Copyright 2024 Nerea Sánchez / GNSS Academy
#
########################################################################
    

import os
import numpy as np
from Ephemeris import decodeEphemeris

SPEED_OF_LIGHT = 299792458


def navParityCheck(bits):
    # This function is called to compute and status the parity bits on GPS word
    # Based on the flowchart in Figure 2-10 in the 2nd Edition of the GPS-SPS
    # Signal Spec.

    # status = navParityCheck(bits)

    #   Inputs:
    #       bits        - an array (1x32) of 32 bits represent a GPS navigation
    #                   word which is 30 bits plus two previous bits used in
    #                   the parity calculation (-2 -1 0 1 2 ... 28 29)

    #   Outputs:
    #       status      - the test value which equals EITHER +1 or -1 if parity
    #                   PASSED or 0 if parity fails. The +1 means bits #1-24
    #                   of the current word have the correct polarity, while -1
    #                   means the bits #1-24 of the current word must be
    #                   inverted.

    # In order to accomplish the exclusive or operation using multiplication
    # this program represents a '0' with a '-1' and a '1' with a '1' so that
    # the exclusive or table holds true for common data operations

    #	a	b	xor 			a	b	product
    #  --------------          -----------------
    #	0	0	 1			   -1  -1	   1
    #	0	1	 0			   -1   1	  -1
    #	1	0	 0			    1  -1	  -1
    #	1	1	 1			    1   1	   1

    # --- Check if the data bits must be inverted ------------------------------
    if bits[1] != 1:
        bits[2:26] *= (-1)

    # --- Calculate 6 parity bits ----------------------------------------------
    # The elements of the bits array correspond to the bits showed in the table
    # 20-XIV (ICD-200C document) in the following way:
    # The first element in the bits is the D29* bit and the second - D30*.
    # The elements 3 - 26 are bits d1-d24 in the table.
    # The elements 27 - 32 in the bits array are the received bits D25-D30.
    # The array "parity" contains the computed D25-D30 (parity) bits.
    parity = np.zeros(6)
    parity[0] = bits[0] * bits[2] * bits[3] * bits[4] * bits[6] * \
                bits[7] * bits[11] * bits[12] * bits[13] * bits[14] * \
                bits[15] * bits[18] * bits[19] * bits[21] * bits[24]

    parity[1] = bits[1] * bits[3] * bits[4] * bits[5] * bits[7] * \
                bits[8] * bits[12] * bits[13] * bits[14] * bits[15] * \
                bits[16] * bits[19] * bits[20] * bits[22] * bits[25]

    parity[2] = bits[0] * bits[2] * bits[4] * bits[5] * bits[6] * \
                bits[8] * bits[9] * bits[13] * bits[14] * bits[15] * \
                bits[16] * bits[17] * bits[20] * bits[21] * bits[23]

    parity[3] = bits[1] * bits[3] * bits[5] * bits[6] * bits[7] * \
                bits[9] * bits[10] * bits[14] * bits[15] * bits[16] * \
                bits[17] * bits[18] * bits[21] * bits[22] * bits[24]

    parity[4] = bits[1] * bits[2] * bits[4] * bits[6] * bits[7] * \
                bits[8] * bits[10] * bits[11] * bits[15] * bits[16] * \
                bits[17] * bits[18] * bits[19] * bits[22] * bits[23] * \
                bits[25]

    parity[5] = bits[0] * bits[4] * bits[6] * bits[7] * bits[9] * \
                bits[10] * bits[11] * bits[12] * bits[14] * bits[16] * \
                bits[20] * bits[23] * bits[24] * bits[25]

    # --- Compare if the received parity is equal to the calculated parity --------
    if (parity == bits[26:]).sum() == 6:
        # Parity is OK. Function output is -1 or 1 depending if the data bits
        # must be inverted or not. The "bits[2]" is D30* bit - the last bit of
        # previous subframe.
        status = -1 * bits[1]

    else:
        # Parity failure
        status = 0

    return status



def findPreambles(settings, trackResults):
    # findPreambles finds the first preamble occurrence in the bit stream of
    # each channel. The preamble is verified by check of the spacing between
    # preambles (6sec) and parity checking of the first two words in a
    # subframe. At the same time function returns list of channels, that are in
    # tracking state and with valid preambles in the nav data stream.

    # firstSubFrameIdx = findPreambles(trackResults, settings)

    #   Inputs:
    #       trackResults    - output from the tracking function
    #       settings        - Receiver settings.

    #   Outputs:
    #       firstSubframe   - the array contains positions of the first
    #                       preamble in each channel. The position is ms count
    #                       since start of tracking. Corresponding value will
    #                       be set to 0 if no valid preambles were detected in
    #                       the channel.

    # Preamble search can be delayed to a later point in the tracking results
    # to avoid noise due to tracking loop transients
    searchStartOffset = 0

    # --- Initialize the firstSubFrameIdx array -----------------------------------
    firstSubFrameIdx = np.zeros(settings.numberOfChannels, dtype=int)

    # --- Generate the preamble pattern ----------------------------------------
    preamble_bits = np.r_[1, - 1, - 1, - 1, 1, - 1, 1, 1]

    # "Upsample" the preamble - make 20 values per bit. The preamble must be
    # found with precision of a sample.
    preamble_ms = np.kron(preamble_bits, np.ones(20))

    # === For all tracking channels ...
    for channelNr in range(len(trackResults)):
        if trackResults[channelNr].status != 0:
            # Correlate tracking output with preamble ================================
            # Read output from tracking. It contains the navigation bits. The start
            # of record is skipped here to avoid tracking loop transients.
            bits = trackResults[channelNr].I_P[searchStartOffset:].copy()

            bits[bits > 0] = 1

            bits[bits <= 0] = - 1

            # zero-pad the preamble so that they are the same length
            tlmXcorrResult = np.correlate(bits,
                                          np.pad(preamble_ms, (0, bits.size - preamble_ms.size), 'constant'),
                                          mode='full')

            # Find all starting points of all preamble like patterns ================
            xcorrLength = int((len(tlmXcorrResult) + 1) / 2)
            IdxPreamble = (np.abs(tlmXcorrResult[xcorrLength - 1:xcorrLength * 2]) > 153).nonzero()[0] + searchStartOffset

            # Analyze detected preamble like patterns ================================
            for i in range(len(IdxPreamble)):
                # --- Find distances in time between this occurrence and the rest of
                # preambles like patterns. If the distance is 6000 milliseconds (one
                # subframe), then do further verifications by validating the parities
                # of two GPS words
                IdxPreambleDelta = IdxPreamble - IdxPreamble[i]

                if (IdxPreambleDelta == 6000).any():
                    # === Re-read bit values for preamble verification ==============
                    # Preamble occurrence is verified by checking the parity of
                    # the first two words in the subframe. Now it is assumed that
                    # bit boundaries are known. Therefore the bit values over 20ms are
                    # combined to increase receiver performance for noisy signals.
                    # in total, 62 bits must be read:
                    # 2 bits from previous subframe are needed for parity checking
                    # 60 bits for the first two 30 bit words (TLM and HOW words)
                    # The IdxPreamble is pointing at the start of TLM word.
                    bits = trackResults[channelNr].I_P[IdxPreamble[i] - 40 : IdxPreamble[i] + 20 * 60].copy()

                    bits = bits.reshape(20, -1, order='F')

                    bits = bits.sum(0)

                    bits[bits > 0] = 1

                    bits[bits <= 0] = - 1

                    if navParityCheck(bits[:32]) != 0 and navParityCheck(bits[30:62]) != 0:
                        # Parity was OK. Record the preamble start position. Skip
                        # the rest of preamble pattern checking for this channel
                        # and process next channel.
                        firstSubFrameIdx[channelNr] = IdxPreamble[i]

                        break
            
            # Exclude channel from the active channel list if no valid preamble was
            # detected
            if firstSubFrameIdx[channelNr] == 0:
                # Exclude channel from further processing. It does not contain any
                # valid preamble and therefore nothing more can be done for it.
                trackResults[channelNr].status = 0
                print('Could not find valid preambles in channel %2d!' % channelNr)

    return firstSubFrameIdx


def buildPseudoranges(settings, msOfTheSignal, trackResults):
    # buildPseudoranges finds relative pseudoranges for all satellites

    # [pseudoranges] = buildPseudoranges(settings, msOfTheSignal, trackResults)

    #   Inputs:
    #       trackResults    - output from the tracking function
    #       msOfTheSignal   - pseudorange measurement point (millisecond) in
    #                       the trackResults structure
    #       settings        - receiver settings

    #   Outputs:
    #       pseudoranges    - relative pseudoranges to the satellites.

    # --- Set initial travel time to infinity ----------------------------------
    # Later in the code a shortest pseudorange will be selected
    # Pseudoranges from non-tracking channels must be the longest => infinite.
    travelTime = np.Inf * np.ones(settings.numberOfChannels)

    # Find number of samples per spreading code
    samplesPerCode = round((settings.samplingFreq * settings.codeLength) / \
                        settings.codeFreqBasis)

    # === For all tracking channels ...
    for channelNr in range(len(trackResults)):
        if trackResults[channelNr].status != 0:        
            # --- Compute the travel times -----------------------------------------
            travelTime[channelNr] = trackResults[channelNr].absoluteSample[
                                        np.int(msOfTheSignal[channelNr])] / samplesPerCode

    # --- Truncate the travelTime and compute pseudoranges ---------------------
    minTravelTime = np.floor(travelTime.min())

    travelTime = settings.zenithTravelTime + (travelTime - minTravelTime)

    # --- Convert travel time to a distance ------------------------------------
    # The speed of light must be converted from meters per second to meters
    # per millisecond.
    pseudoranges = travelTime * SPEED_OF_LIGHT / 1000
    return pseudoranges


def buildMeas(settings, trackResults):
    # Function calculates the pseudorange measurements
    # [eph] = postNavigation(settings, trackResults)

    #   Inputs:
    #       trackResults    - results from the tracking function (structure
    #                       array).
    #       settings        - receiver settings.
    #   Outputs:
    #       eph             - received ephemerides of all SV (structure array).

    # Check is there are enough data to obtain any navigation solution ===========
    # It is necessary to have at least three subframes (number 1, 2 and 3) to
    # find satellite coordinates. Then receiver position can be found too.
    # The function requires all 5 subframes, because the tracking starts at an
    # arbitrary point. Therefore the first received subframes can be any three
    # from the 5.
    # One subframe length is 6 seconds, therefore we need at least 30 sec long
    # record (5 * 6 = 30 sec = 30000ms). We add extra seconds for the cases,
    # when tracking has started in a middle of a subframe.

    print('Building measurements...')

    if (settings.msToProcess < 36000 or \
        sum([trackResults[i].status != 0 for i in range(settings.numberOfChannels)]) < 4):
        # Show the error message and exit
        print('Record is too short (<36 seconds) or too few satellites tracked. Exiting!')
        return

    # Find preamble start positions ==========================================
    subFrameStartIdx = findPreambles(settings, trackResults)

    # Decode ephemerides =====================================================
    nActiveChn = 0
    field_str = 'weekNumber,accuracy,health,T_GD,IODC,t_oc,a_f2,a_f1,a_f0,'
    field_str += 'IODE_sf2,C_rs,deltan,M_0,C_uc,e,C_us,sqrtA,t_oe,'
    field_str += 'C_ic,omega_0,C_is,i_0,C_rc,omega,omegaDot,IODE_sf3,iDot'
    eph = np.recarray((32,), formats=['O'] * 27, names=field_str)
    for channelNr in range(len(trackResults)):
        if trackResults[channelNr].status != 0:
            # === Convert tracking output to navigation bits =======================
            # --- Copy 5 sub-frames from tracking output ---------------------------
            navBitsSamples = trackResults[channelNr].I_P[subFrameStartIdx[channelNr] - 20:
                                                            subFrameStartIdx[channelNr] + 1500 * 20].copy()

            navBitsSamples = navBitsSamples.reshape(20, -1, order='F')

            navBits = navBitsSamples.sum(0)

            # The expression (navBits > 0) returns an array with elements set to 1
            # if the condition is met and set to 0 if it is not met.
            navBits = (navBits > 0) * 1
            navBitsBin = navBits.astype(str)

            # Decode ephemeris
            eph[trackResults[channelNr].PRN - 1], TOW = decodeEphemeris(navBitsBin[1:], navBitsBin[0])

            if eph[trackResults[channelNr].PRN - 1].IODC is None or \
                    eph[trackResults[channelNr].PRN - 1].IODE_sf2 is None or \
                    eph[trackResults[channelNr].PRN - 1].IODE_sf3 is None:
                # --- Exclude channel from the list (from further processing) ------
                trackResults[channelNr].status = 0

            else:
                nActiveChn = nActiveChn + 1

    # Check if the number of satellites is still above 3 =====================
    if nActiveChn < 4:
        # Show error message and exit
        print('Too few satellites with ephemeris data for position calculations. Exiting!')
        return

    # Initialization =========================================================
    TransmissionTime = TOW

    ###########################################################################
    #   Build measurements                                                    #
    ###########################################################################
    filepath = os.path.dirname(settings.inputf) + \
    '/SW-RCVR-PYTHON/MEAS_%s' % os.path.basename(settings.inputf)
    with open(filepath, 'w') as f:
        for currMeasNr in range(np.int(np.fix(settings.msToProcess - subFrameStartIdx.max()) / settings.navSolPeriod)):
            # Find pseudoranges ======================================================
            rawP = buildPseudoranges(settings,
                                     subFrameStartIdx + settings.navSolPeriod * currMeasNr,
                                     trackResults)

            # dump measurements to file
            for channelNr in range(len(trackResults)):
                if trackResults[channelNr].status != 0:
                    f.write("%5d %4s %02d %16.4f\n" % \
                        (currMeasNr, "G", trackResults[channelNr].PRN, rawP[channelNr])) 
        
    return eph
