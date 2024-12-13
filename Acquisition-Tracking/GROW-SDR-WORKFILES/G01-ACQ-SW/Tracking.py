

########################################################################
# Tracking.py:
# Tracking module
#
#  Project:        SW-RCVR-PYTHON
#  File:           Tracking.py
#  Date(YY/MM/DD): 30/03/22
#
#   Author: Nerea Sánchez
#   Copyright 2024 Nerea Sánchez / GNSS Academy
#
########################################################################


import numpy as np

from GoldCodes import generateGoldCode
from GalileoCodes import gale1bGeneratePrnCode, gale1bModulatePrnCode


def loadChannels(settings, acqResults, channels):
    #--- Sort peaks to find strongest signals, keep the peak index information
    indices = np.argsort(acqResults.peakMetric)[::-1]

    #--- Load information about each satellite --------------------------------
    # Maximum number of initialized channels depends on the number of 
    # detected signals, it's not greater than the number of channels 
    # specified in the settings.
    for ii in range(min([settings.numberOfChannels, sum(acqResults.carrFreq > 0)])):
        channels.PRN[ii]           = indices[ii]
        channels.acquiredFreq[ii]  = acqResults.carrFreq[indices[ii]]
        channels.acquiredDoppler[ii]  = (acqResults.carrFreq[indices[ii]] \
                                                    - settings.IF)/1000
        channels.acquiredDelay[ii] = acqResults.codeDelay[indices[ii]]
        channels.acquiredDelayChips[ii] = acqResults.codeDelay[indices[ii]] \
                                            * settings.codeFreqBasis / settings.samplingFreq
        channels.SNR[ii]  = acqResults.SNR[indices[ii]]
        # Set tracking mode
        channels.status[ii]       = 1

    # end of for ii in range(min([settings.numberOfChannels, sum(acqResults.carrFreq > 0)]))

    # Display Channels Status
    print("%4s %02s %19s %10s %7s %9s" % ("ChId", "PRN", "Doppler [kHz] ", "Delay [chip]", "SNR", "Status"))
    print("-"*60)
    for ichannel, line in enumerate(
        np.array([channels.PRN, channels.acquiredDoppler, channels.acquiredDelayChips, channels.SNR, channels.status]).T):
        print("%4d " % (ichannel) + " %02d %19.3f %10.3f %10.3f %9d" % tuple(line))
    print()


def calcLoopCoef(LBW, zeta, k):
    # Solve natural frequency
    Wn = LBW*8*zeta / (4*zeta**2 + 1)

    # solve for t1 & t2
    tau1 = k / (Wn * Wn)
    tau2 = 2.0 * zeta / Wn

    return [tau1, tau2]

# end of calcLoopCoef

def tracking(settings, rawSignal, channels):

    class TrackingResults():
        def __init__(self, i):
            self.ChId = i

            # PRN
            self.PRN            =  0

            # Channel status
            self.status         =  0       # No tracked signal, or lost lock

            # The absolute sample in the record of the C/A code start:
            self.absoluteSample = np.zeros(settings.msToProcess)

            # Freq of the C/A code:
            self.codeFreq       = np.inf * np.ones(settings.msToProcess)

            # Frequency of the tracked carrier wave:
            self.carrFreq       = np.inf * np.ones(settings.msToProcess)

            # Outputs from the correlators (In-phase):
            self.I_P            = np.inf * np.ones(settings.msToProcess)
            self.I_E            = np.inf * np.ones(settings.msToProcess)
            self.I_L            = np.inf * np.ones(settings.msToProcess)

            # Outputs from the correlators (Quadrature-phase):
            self.Q_E            = np.inf * np.ones(settings.msToProcess)
            self.Q_P            = np.inf * np.ones(settings.msToProcess)
            self.Q_L            = np.inf * np.ones(settings.msToProcess)

            # Loop discriminators
            self.dllDiscr       = np.inf * np.ones(settings.msToProcess)
            self.dllDiscrFilt   = np.inf * np.ones(settings.msToProcess)
            self.pllDiscr       = np.inf * np.ones(settings.msToProcess)
            self.pllDiscrFilt   = np.inf * np.ones(settings.msToProcess)

    # end of TrackingResults

    def plot(settings, trackResults):
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gs
        import os

        # This function plots the tracking results for the given channel list.

        # plotTracking(settings, trackResults)

        #   Inputs:
        #       trackResults    - tracking results from the tracking function.
        #       settings        - receiver settings.

        mpl.rcdefaults()
        # mpl.rcParams['font.sans-serif']
        # mpl.rcParams['font.family'] = 'serif'
        mpl.rc('savefig', bbox='tight', transparent=False, format='png')
        mpl.rc('axes', grid=True, linewidth=1.5, axisbelow=True)
        mpl.rc('lines', linewidth=1.5, solid_joinstyle='bevel')
        mpl.rc('figure', figsize=[25, 16], autolayout=False, dpi=240)
        # mpl.rc('text', usetex=True)
        mpl.rc('font', size=10)
        mpl.rc('mathtext', fontset='cm')

        # mpl.rc('font', size=16)
        # mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')

        channelList = range(settings.numberOfChannels)

        # Protection - if the list contains incorrect channel numbers
        channelList = np.intersect1d(channelList, range(settings.numberOfChannels))

        # === For all listed channels ==============================================
        for channelNr in channelList:
            # Select (or create) and clear the figure ================================
            # The number 200 is added just for more convenient handling of the open
            # figure windows, when many figures are closed and reopened.
            # Figures drawn or opened by the user, will not be "overwritten" by
            # this function.
            f = plt.figure(channelNr + 200)
            f.suptitle('Channel ' + str(channelNr) +
                        ' (PRN%02d' % trackResults[channelNr].PRN + ') results',
                        size=20)
            # Draw axes ==============================================================
            # Row 1
            spec = gs.GridSpec(3, 3)
            h11 = plt.subplot(spec[0, 0])

            h12 = plt.subplot(spec[0, 1:])

            h21 = plt.subplot(spec[1, 0])

            h22 = plt.subplot(spec[1, 1:])

            h31 = plt.subplot(spec[2, 0])

            h32 = plt.subplot(spec[2, 1])

            h33 = plt.subplot(spec[2, 2])

            # Plot all figures =======================================================
            timeAxisInSeconds = np.arange(settings.msToProcess) / 1000.0

            h11.plot(trackResults[channelNr].I_P, trackResults[channelNr].Q_P, '.')
            h11.grid()
            h11.axis('equal')
            h11.set(title='Discrete-Time Scatter Plot', xlabel='I prompt', ylabel='Q prompt')
            
            h12.plot(timeAxisInSeconds, trackResults[channelNr].I_P)
            h12.grid()
            h12.set(title='Bits of the navigation message', xlabel='Time (s)')
            h12.axis('tight')
            
            h21.plot(timeAxisInSeconds, trackResults[channelNr].pllDiscr, 'r')
            h21.grid()
            h21.axis('tight')
            h21.set(xlabel='Time (s)', ylabel='Amplitude', title='Raw PLL discriminator')
            
            h22.plot(timeAxisInSeconds,
                     np.sqrt(trackResults[channelNr].I_E ** 2 + trackResults[channelNr].Q_E ** 2),
                     timeAxisInSeconds,
                     np.sqrt(trackResults[channelNr].I_P ** 2 + trackResults[channelNr].Q_P ** 2),
                     timeAxisInSeconds,
                     np.sqrt(trackResults[channelNr].I_L ** 2 + trackResults[channelNr].Q_L ** 2), '-*')
            h22.grid()
            h22.set(title='Correlation results', xlabel='Time (s)')
            h22.axis('tight')
            h22.legend(['$\sqrt{I_{E}^2 + Q_{E}^2}$', '$\sqrt{I_{P}^2 + Q_{P}^2}$',
                        '$\sqrt{I_{L}^2 + Q_{L}^2}$'])

            h31.plot(timeAxisInSeconds, trackResults[channelNr].pllDiscrFilt, 'b')
            h31.grid()
            h31.axis('tight')
            h31.set(xlabel='Time (s)',
                    ylabel='Amplitude',
                    title='Filtered PLL discriminator')
            
            h32.plot(timeAxisInSeconds, trackResults[channelNr].dllDiscr, 'r')
            h32.grid()
            h32.axis('tight')
            h32.set(xlabel='Time (s)',
                    ylabel='Amplitude',
                    title='Raw DLL discriminator')
            
            h33.plot(timeAxisInSeconds, trackResults[channelNr].dllDiscrFilt, 'b')
            h33.grid()
            h33.axis('tight')
            h33.set(xlabel='Time (s)',
                    ylabel='Amplitude',
                    title='Filtered DLL discriminator')
            
            figspath = os.path.dirname(settings.inputf) + '/SW-RCVR-PYTHON/'
            plt.savefig(figspath + 'TRACKING_CH%02d.png' % channelNr)

        # end of for channelNr in channelList

        plt.close('all')

    # end of plot(settings, trackResults)


    ## Initialize result structure ============================================
    trackResults = []

    #--- Copy initial settings for all channels -------------------------------
    for i in range(settings.numberOfChannels):
        trackResults.append(TrackingResults(i))

    ## Initialize tracking variables ==========================================
    nCodesToProcess = settings.msToProcess     # For GPS L1 C/A: one C/A code is one ms

    #--- DLL variables --------------------------------------------------------
    # Define early-late semi-distance (in chips)
    earlyLateSemi = settings.dllCorrelatorSpacing / 2

    # Summation interval
    PDIcode = settings.Nc

    # Calculate filter coefficient values
    [tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth,
                                        settings.dllDampingRatio,
                                        1.0)

    #--- PLL variables --------------------------------------------------------
    # Summation interval
    PDIcarr = settings.Nc

    # Calculate filter coefficient values
    [tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth,
                                        settings.pllDampingRatio,
                                        0.25)
    
    print('Tracking...')

    ## Start processing channels ==============================================
    for channelNr in range(settings.numberOfChannels):
        
        # Only process PRN if acquisition was successful
        if (channels.status[channelNr] != 0):
            # Save information on each channels's tracked PRN
            trackResults[channelNr].PRN = channels.PRN[channelNr]

            # Skip through the data file to start at the appropriate sample 
            # (corresponding to the beginning of the code in the incoming signal)
            # => the Code Phase
            dataPtr = int(channels.acquiredDelay[channelNr])

            # Generate code replica
            codeReplica = generateGoldCode(channels.PRN[channelNr])
            
            # Then make it possible to do early and late versions
            codeReplica = np.concatenate(([codeReplica[-1]], codeReplica, [codeReplica[1]]))

            #--- Perform various initializations ------------------------------

            # define initial code frequency basis of NCO
            codeFreq      = settings.codeFreqBasis
            
            # define residual code phase (in chips)
            remCodeDelay  = 0.0
            
            # define carrier frequency which is used over whole tracking period
            carrFreq      = channels.acquiredFreq[channelNr]
            carrFreqBasis = channels.acquiredFreq[channelNr]
            
            # define residual carrier phase
            remCarrPhase  = 0.0

            # code tracking loop (DLL) parameters
            oldCodeNco   = 0.0
            oldCodeError = 0.0

            # carrier/Costas loop (PLL) parameters
            oldCarrNco   = 0.0
            oldCarrError = 0.0

            #=== Process the number of specified code periods =================
            for loopCnt in range(nCodesToProcess):

                # Update the Code Phase step based on code freq (variable) 
                # and sampling frequency (constant)
                codePhaseStep = codeFreq / settings.samplingFreq
                
                blksize = np.ceil(
                    (settings.codeLength-remCodeDelay) / codePhaseStep
                    ).astype(int)

                # If we ran out of data, continue to next channel
                endDataBlock = dataPtr + blksize
                if (endDataBlock > len(rawSignal)):
                    continue
                
                # Read in the appropriate number of samples to process this
                # iteration
                samplesRead = rawSignal[dataPtr:endDataBlock]
                # DBG
                # if channels.PRN[channelNr]==21: 
                #     print("samplesRead PRN", channels.PRN[channelNr], dataPtr, blksize, 
                #     samplesRead[:10],samplesRead[-10:],)
                # DBG

                dataPtr = endDataBlock

                ## Set up all the code delay tracking information -------------------------
                # Define index into *** EARLY *** code vector
                tcode = np.linspace(remCodeDelay - earlyLateSemi,
                                    blksize * codePhaseStep + remCodeDelay - earlyLateSemi,
                                    blksize, 
                                    endpoint=False)

                tcode2      = np.ceil(tcode).astype(int)
                earlyCode   = codeReplica[tcode2]
                
                # Define index into *** LATE *** code vector
                tcode = np.linspace(remCodeDelay + earlyLateSemi,
                                    blksize * codePhaseStep + remCodeDelay + earlyLateSemi,
                                    blksize, 
                                    endpoint=False)

                tcode2      = np.ceil(tcode).astype(int)
                lateCode    = codeReplica[tcode2]
                
                # Define index into *** PROMPT *** code vector
                tcode = np.linspace(remCodeDelay,
                                    blksize * codePhaseStep + remCodeDelay,
                                    blksize, endpoint=False)

                tcode2      = np.ceil(tcode).astype(int)
                promptCode  = codeReplica[tcode2]
                
                remCodeDelay = (tcode[blksize-1] + codePhaseStep) - settings.codeLength
                # DBG
                # if channels.PRN[channelNr]==21: print("remCodeDelay", earlyCode[:3],lateCode[:3],promptCode[:3],remCodeDelay)
                # DBG

                ## Generate the carrier frequency -----------------------------------------
                time    = np.arange(0, blksize+1) / settings.samplingFreq

                # Get the argument to sin/cos functions
                trigarg = ((carrFreq * 2.0 * np.pi) * time) + remCarrPhase
                remCarrPhase = trigarg[blksize] % (2 * np.pi)

                # Finally, compute the replica of the Carrier
                carrCosReplica = np.cos(trigarg[:blksize])
                carrSinReplica = np.sin(trigarg[:blksize])
                
                ## Remove the Carrier ---------------------------
                qCodeData = carrCosReplica * samplesRead
                iCodeData = carrSinReplica * samplesRead

                # Now get early, late, and prompt values for each
                I_E = np.sum(earlyCode  * iCodeData)
                Q_E = np.sum(earlyCode  * qCodeData)
                I_P = np.sum(promptCode * iCodeData)
                Q_P = np.sum(promptCode * qCodeData)
                I_L = np.sum(lateCode   * iCodeData)
                Q_L = np.sum(lateCode   * qCodeData)

                ## Find PLL error and update carrier NCO ----------------------------------

                # Implement carrier loop discriminator
                carrError = np.arctan(Q_P / I_P) / (2.0 * np.pi)
                
                # Implement carrier loop filter and generate NCO command
                carrNco = oldCarrNco + \
                    (PDIcarr/tau1carr) * carrError + \
                    (tau2carr/tau1carr) * (carrError - oldCarrError)

                oldCarrNco   = carrNco
                oldCarrError = carrError

                # Modify carrier freq based on NCO command
                carrFreq = carrFreqBasis + carrNco

                trackResults[channelNr].carrFreq[loopCnt] = carrFreq

                ## Find DLL error and update code NCO -------------------------------------
                # Implement code loop discriminator
                codeError = (np.sqrt(I_E * I_E + Q_E * Q_E) - np.sqrt(I_L * I_L + Q_L * Q_L)) /\
                    (np.sqrt(I_E * I_E + Q_E * Q_E) + np.sqrt(I_L * I_L + Q_L * Q_L))

                # Implement code loop filter and generate NCO command
                codeNco = oldCodeNco + \
                    (PDIcode/tau1code) * codeError + \
                    (tau2code/tau1code) * (codeError - oldCodeError) 

                oldCodeNco   = codeNco
                oldCodeError = codeError
                
                # Modify code freq based on NCO command
                codeFreq = settings.codeFreqBasis - codeNco
                
                trackResults[channelNr].codeFreq[loopCnt] = codeFreq

                ## Record various measurements to show in postprocessing ----------------------
                trackResults[channelNr].absoluteSample[loopCnt] = dataPtr

                trackResults[channelNr].dllDiscr[loopCnt]       = codeError
                trackResults[channelNr].dllDiscrFilt[loopCnt]   = codeNco
                trackResults[channelNr].pllDiscr[loopCnt]       = carrError
                trackResults[channelNr].pllDiscrFilt[loopCnt]   = carrNco

                trackResults[channelNr].I_E[loopCnt] = I_E
                trackResults[channelNr].I_P[loopCnt] = I_P
                trackResults[channelNr].I_L[loopCnt] = I_L
                trackResults[channelNr].Q_E[loopCnt] = Q_E
                trackResults[channelNr].Q_P[loopCnt] = Q_P
                trackResults[channelNr].Q_L[loopCnt] = Q_L

                # DBG
                # print("trackResults",
                # "PRN", trackResults[channelNr].PRN,
                # "loopCnt", loopCnt,
                # "dllDiscr", trackResults[channelNr].dllDiscr[loopCnt],
                # "pllDiscr", trackResults[channelNr].pllDiscr[loopCnt],
                # "blksize", blksize,
                # )
                # DBG

            # end # for loopCnt

            # If we got so far, this means that the tracking was successful
            # Now we only copy status, but it can be update by a lock detector
            # if implemented
            trackResults[channelNr].status  = channels.status[channelNr]

        # end # if a PRN is assigned
    # end # for channelNr 

    plot(settings, trackResults)

    return trackResults
