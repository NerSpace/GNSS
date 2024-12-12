

########################################################################
# Acquisition.py:
# Acquisition module
#
#  Project:        sw-rcvr-python
#  File:           Acquisition.py
#
#   Author: Nerea Sánchez / GNSS Academy
#   Copyright 2024 Nerea Sánchez / GNSS Academy
#
########################################################################


import numpy as np
from scipy.fftpack import fft, ifft

np.set_printoptions(threshold=10000000)

from GoldCodes import generateGoldCode
from GalileoCodes import gale1bGeneratePrnCode, gale1bModulatePrnCode,\
    gale1cGeneratePrnCode, gale1cModulatePrnCode
from Misc import nextPowerOf2


class AcqResults():
    def initialize(self, settings):
        # Find number of samples per Ranging Code
        samplesPerCode = round((settings.samplingFreq * settings.codeLength) / \
                            settings.codeFreqBasis)

        # Get Number of frequency bins for the given acquisition band (125Hz steps)
        nFrqBins = round(settings.acqFreqRangekHz * 8) + 1

        # Search space
        self.searchSpace   = np.zeros((len(settings.satMask)+1, nFrqBins, samplesPerCode))

        # Carrier frequencies of detected signals
        self.carrFreq      = np.zeros(len(settings.satMask)+1)

        # Code delays of detected signals
        self.codeDelay     = np.zeros(len(settings.satMask)+1)

        # Correlation peak ratios of the detected signals
        self.peakMetric    = np.zeros(len(settings.satMask)+1)

        # Signal to Noise ratio
        self.SNR    = np.zeros(len(settings.satMask)+1)

    def plot(self, settings, samplesPerCode, nFrqBins):
        # assert isinstance(self._results, np.recarray)
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        # from scipy.io.matlab import loadmat
        import os
        from mpl_toolkits.mplot3d import Axes3D

        figspath = os.path.dirname(settings.inputf) + '/SW-RCVR-PYTHON/'
        if not os.path.isdir(figspath):
            os.makedirs(figspath)

        mpl.rcdefaults()
        mpl.rc('savefig', bbox='tight', transparent=False, format='png')
        mpl.rc('axes', linewidth=1.5, axisbelow=True)
        mpl.rc('lines', linewidth=1.5, solid_joinstyle='bevel')
        mpl.rc('figure', figsize=[16, 9], autolayout=False, dpi=240)
        mpl.rc('font', size=10)
        mpl.rc('mathtext', fontset='cm')

        #   Inputs:
        #       settings      - Tool settings.
        #       acqResults    - Acquisition results from function acquisition.

        # Plot all results =======================================================
        # Plot Acquisition search-space
        if(True):
            for PRN in range(1, len(settings.satMask)+1):
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')

                frequencies = (-settings.acqFreqRangekHz/2) + \
                                    0.125 * np.arange(nFrqBins)
                delays = np.arange(samplesPerCode) * settings.codeFreqBasis / settings.samplingFreq
                frequencies, delays = np.meshgrid(frequencies, delays)

                surf = ax.plot_surface(frequencies, delays, self.searchSpace[PRN].T,
                cmap=mpl.cm.coolwarm,
                antialiased=False)
                plt.title('PRN %02d Search Space' % PRN)
                plt.xlabel('Doppler Frequency [kHz]')
                plt.ylabel('Code Delay')
                oldAxis = plt.axis()

                plt.grid(linestyle='--', linewidth=0.5)

                plt.savefig(figspath + 'SEARCH_SPACE_PRN%02d.png' % PRN)
                plt.close('all')

            # end for prn in range(1, len(settings.satMask)+1):

        # Plot Acquisition metric
        f, hAxes = plt.subplots()

        plt.bar(range(1, len(settings.satMask)+1), self.peakMetric[1:])
        plt.title('Acquisition results')
        plt.xlabel('PRN number (no bar - SV is not in the acquisition list)')
        plt.ylabel('Acquisition Metric')
        oldAxis = plt.axis()

        plt.axis([0, len(settings.satMask)+1, 0, oldAxis[-1]])
        plt.xticks(range(1, len(settings.satMask)+1), size=12)
        # plt.minorticks_on()
        plt.grid(linestyle='--', linewidth=0.5)
        # Mark acquired signals ==================================================

        acquiredSignals = self.peakMetric[1:] * (self.carrFreq[1:] > 0)

        plt.bar(range(1, len(settings.satMask)+1), acquiredSignals)
        plt.legend(['Not acquired signals', 'Acquired signals'])
        plt.savefig(figspath + 'ACQUISITION_METRIC.png')

        # Plot C/N0
        f, hAxes = plt.subplots()

        plt.bar(range(1, len(settings.satMask)+1), self.SNR[1:])
        plt.title('C/N0')
        plt.xlabel('PRN number (no bar - SV is not in the acquisition list)')
        plt.ylabel('Signal to Noise Ratio [dB-Hz]')
        oldAxis = plt.axis()

        plt.axis([0, len(settings.satMask)+1, 0, oldAxis[-1]])
        plt.xticks(range(1, len(settings.satMask)+1), size=12)
        # plt.minorticks_on()
        plt.grid(linestyle='--', linewidth=0.5)
        # Mark acquired signals ==================================================

        acquiredSignals = self.SNR[1:] * (self.carrFreq[1:] > 0)

        plt.bar(range(1, len(settings.satMask)+1), acquiredSignals)

        ylim = plt.gca().get_ylim()
        plt.ylim([20, ylim[1]])

        plt.legend(['Not acquired signals', 'Acquired signals'])
        plt.savefig(figspath + 'SNR.png')

        plt.close('all')



def acquisitionGpsL1C(settings, inputSignal):
    # Find number of samples per Ranging Code
    samplesPerCode = round((settings.samplingFreq * settings.codeLength) / \
                        settings.codeFreqBasis)

    # Find number of samples per Ranging Code Chip
    samplesPerCodeChip = round(settings.samplingFreq / settings.codeFreqBasis)

    # Create two 1ms-long vectors of data to correlate with and one with zero DC
    signal1 = inputSignal[:samplesPerCode]
    signal2 = inputSignal[samplesPerCode : 2*samplesPerCode]
    signal0DC = inputSignal - np.mean(inputSignal) 

    # Find sampling periods
    ts = 1 / settings.samplingFreq
    tc = 1 / settings.codeFreqBasis

    # Find phase points of the local carrier wave at 2*k*pi*ts
    phasePoints = 2 * np.arange(samplesPerCode, dtype=int) * np.pi * ts

    # Get Number of frequency bins for the given acquisition band 
    # (125Hz steps => 8 steps per kHz)
    # Recall: to avoid correlation colapse => step < 1/(2T) = 1/(2*1ms) = 500Hz
    nFrqBins = round(settings.acqFreqRangekHz * 8) + 1

    #--- Initialize arrays to speed up the code -------------------------------
    # Search results of all frequency bins and code shifts (for one satellite)
    results = np.zeros((nFrqBins, samplesPerCode))

    # Carrier frequencies of the frequency bins
    frqBins     = np.zeros(nFrqBins)

    #--- Make index array to read C/A code values -------------------------
    # The length of the index array depends on the sampling frequency -
    # number of samples per millisecond (because one C/A code period is one
    # millisecond).
    codeOversampIdx = np.ceil((ts * np.arange(1, samplesPerCode+1)) / tc).astype(int) - 1

    #--- Correct the last index (due to number rounding issues) -----------
    codeOversampIdx[-1] = settings.codeLength-1

    #--- Get indices to generate 10msec long C/A codes sequence for a given PRN --------
    codeOversampIdx10 = np.floor(
        (ts * settings.codeFreqBasis) * np.arange(1,10*samplesPerCode+1)
        ).astype(int)

    #--- Initialize acqResults ------------------------------------------------
    acqResults = AcqResults()
    acqResults.initialize(settings)

    delayBinIndex = 0
    frequencyBinIndex = 0
    peakSize = np.nan

    print('Acquiring GPS L1C ...\n(', end='')

    #=== For all satellite PRN-s \
    for PRN in settings.satMask:
        # Generate CA code for given PRN
        caCodeReplica = generateGoldCode(PRN)

        # ====================== Digitizing =====================

        # Make the digitized version of the C/A code (oversampling)
        # The unsampled code is made by selecting values from the CA code 
        # chip array (caCodeReplica) for the time instances of each sample.
        caCodeReplicaOversamp = caCodeReplica[codeOversampIdx]

        # ====================== Correlate signals =====================
        # Perform DFT of C/A code
        caCodeReplicaFreqDom = np.conj(fft(caCodeReplicaOversamp, ))

        #--- Make the correlation for whole frequency band (for all freq. bins)
        for frqBinIndex in range(nFrqBins):
            #--- Generate carrier wave frequency grid (125 Hz step)
            frqBins[frqBinIndex] = settings.IF - (settings.acqFreqRangekHz/2)*1000 + 125*frqBinIndex

            #--- Generate local sine and cosine
            sinCarrReplica = np.sin(frqBins[frqBinIndex] * phasePoints)
            cosCarrReplica = np.cos(frqBins[frqBinIndex] * phasePoints)

            #--- Remove carrier from signal (demodulation)
            I1 = sinCarrReplica * signal1
            Q1 = cosCarrReplica * signal1
            I2 = sinCarrReplica * signal2
            Q2 = cosCarrReplica * signal2

            #--- Convert the signal to frequency domain 
            IQfreqDom1 = fft(I1+Q1*1j)
            IQfreqDom2 = fft(I2+Q2*1j)

            #--- Multiplication in the frequency domain (correlation in time domain) --> "Time loop"
            convCodeIQ1 = IQfreqDom1 * caCodeReplicaFreqDom
            convCodeIQ2 = IQfreqDom2 * caCodeReplicaFreqDom

            #--- Perform inverse DFT and store correlation results
            acqRes1 = abs(ifft(convCodeIQ1))**2
            acqRes2 = abs(ifft(convCodeIQ2))**2

            #--- Check which msec had the greater power and save that one 
            # it will "blend" 1st and 2nd msec but will correct data bit issues
            if (max(acqRes1) > max(acqRes2)):
                results[frqBinIndex] = acqRes1
            else:
                results[frqBinIndex] = acqRes2
        # End of for frqBinIndex in range(nFrqBins)

        # ===================== Look for correlation peaks in the results =====================
        # Find the highest peak and compare it to the second highest peak 
        # The second highest peak is chosen not closer than 1 chip to the highest peak

        #--- Store result ----------------------------------------------------------
        acqResults.searchSpace[PRN] = results

        #--- Find the correlation peak ---------------------------------------------
        peakSize = results.max()
        indices = np.where(results == peakSize)
        frequencyBinIndex = indices[0][0]

        delayBinIndex = indices[1][0]

        #--- Find 1 chip wide C/A delay to exclude 2 chips around the peak ----
        excludeRangeIndex1 = delayBinIndex - samplesPerCodeChip
        excludeRangeIndex2 = delayBinIndex + samplesPerCodeChip

        #--- Correct C/A code phase excluding range if the range includes array
        # boundaries
        if excludeRangeIndex1 < 1:
            NoiseIdx = range(excludeRangeIndex2,
                            (samplesPerCode + excludeRangeIndex1 + 1))
                            
        elif excludeRangeIndex2 > samplesPerCode:
            NoiseIdx = range((excludeRangeIndex2 - samplesPerCode),
                            excludeRangeIndex1 + 1)
        else:
            NoiseIdx = list(range(excludeRangeIndex1)) + \
                            list(range(excludeRangeIndex2,samplesPerCode))
        
        # end of if excludeRangeIndex1 < 1

        #--- Find the second highest correlation peak in the same freq. bin ---
        secondPeakSize = max(results[frequencyBinIndex, NoiseIdx])

        #--- Store result -----------------------------------------------------
        acqResults.peakMetric[PRN] = peakSize/secondPeakSize
        
        #--- C/N0 -------------------------------------------------------------
        # Estimate Noise level
        noiseLevel = np.mean(results[frequencyBinIndex, NoiseIdx])

        # Compute C/N0
        CNo = (peakSize-noiseLevel)/noiseLevel/1e-3
        acqResults.SNR[PRN] = 10 * np.log10(CNo)

        # If the result is above threshold, then the satellite is acquired
        if (peakSize/secondPeakSize) > settings.acqTh:
            ## Fine resolution frequency search =======================================
            
            #--- Indicate PRN number of the detected signal -------------------
            print(' %02d ' % PRN, end='')
            
            #--- Generate 10msec long C/A codes sequence for given PRN --------
            longCaCodeReplica = caCodeReplica[codeOversampIdx10 % settings.codeLength]
        
            #--- Remove C/A code modulation from the original signal ----------
            # (Using detected C/A code phase)
            IncomingCarrier = signal0DC[range(delayBinIndex,(delayBinIndex + 10*samplesPerCode))] \
                * longCaCodeReplica
            
            #--- Find the next highest power of two and increase by 8x --------
            fftNumPts = 8*(nextPowerOf2(len(IncomingCarrier)))
            
            #--- Compute the magnitude of the FFT, find maximum and the
            # associated carrier frequency 
            fftxc = abs(fft(IncomingCarrier, fftNumPts)) 
            
            uniqFftPts = np.ceil((fftNumPts + 1) / 2).astype(int)
            fftMaxIndex = np.argmax(fftxc[4 : uniqFftPts-6])
            # DBG
            # print("fftxc PRN", PRN, fftxc.shape, fftxc[:5])
            # print("fftMaxIndex PRN", PRN, fftMaxIndex)
            # DBG
            
            fftFreqBins = (settings.samplingFreq/fftNumPts) * np.arange(uniqFftPts)
            
            #--- Save properties of the detected satellite signal -------------
            acqResults.carrFreq[PRN]  = fftFreqBins[fftMaxIndex]
            acqResults.codeDelay[PRN] = delayBinIndex

        else:
            #--- No signal with this PRN --------------------------------------
            print(' _ ', end='')
        # end   # if (peakSize/secondPeakSize) > settings.acqThreshold

    # end    # for PRN = satelliteList

    #=== Acquisition is over ==================================================
    print(')\n')

    #=== Plot results =========================================================
    acqResults.plot(settings, samplesPerCode, nFrqBins)

    return acqResults


def acquisitionGalE1B(settings, inputSignal):
    # Find number of samples per Ranging Code
    samplesPerCode = round((settings.samplingFreq * settings.codeLength) / \
                        settings.codeFreqBasis)

    # Find number of samples per Ranging Code Chip
    samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis)

    # Create two 4ms-long vectors of data to correlate with and one with zero DC
    signal1 = inputSignal[:samplesPerCode]
    signal2 = inputSignal[samplesPerCode : 2*samplesPerCode]

    # Find sampling periods
    ts = 1 / settings.samplingFreq
    # correct to modulated code period
    tc = 1 / (12 * settings.codeFreqBasis)

    # Find phase points of the local carrier wave at 2*k*pi*ts
    phasePoints = 2 * np.arange(samplesPerCode, dtype=int) * np.pi * ts

    # Get Number of frequency bins for the given acquisition band (125Hz steps)
    # Recall: to avoid correlation colapse => step < 1/(2T) = 1/(2*4ms) = 125Hz
    nFrqBins = round(settings.acqFreqRangekHz * 8) + 1


    #--- Initialize arrays to speed up the code -------------------------------
    # Search results of all frequency bins and code shifts (for one satellite)
    results = np.zeros((nFrqBins, samplesPerCode))

    # Carrier frequencies of the frequency bins
    frqBins     = np.zeros(nFrqBins)

    #--- Make index array to read E1B code values -------------------------
    # The length of the index array depends on the sampling frequency -
    # number of samples per millisecond (because one E1B code period is one
    # millisecond).
    codeOversampIdx = np.ceil((ts * np.arange(1, samplesPerCode+1)) / tc).astype(int) - 1

    #--- Correct the last index (due to number rounding issues) -----------
    codeOversampIdx[-1] = settings.codeLength-1

    #--- Initialize acqResults ------------------------------------------------
    acqResults = AcqResults()
    acqResults.initialize(settings)

    print('Acquiring Galileo E1B ...\n(', end='')

    #=== For all satellite PRN-s \
    for PRN in settings.satMask:

        # -- Generate Ranging codes for a given PRN
        galCodeE1Brep = gale1bGeneratePrnCode(PRN)
        galCodeE1Brep = gale1bModulatePrnCode(settings, galCodeE1Brep)
        galCodeE1Crep = gale1cGeneratePrnCode(PRN)
        galCodeE1Crep = gale1cModulatePrnCode(settings, galCodeE1Crep)

        # ====================== Digitizing =====================

        # Make the digitized version of the E1B code (oversampling)
        # The unsampled code is made by selecting values from the Ranging Code
        # chip array (galCodeE1Brep) for the time instances of each sample.
        galCodeE1BrepOversamp = galCodeE1Brep[codeOversampIdx]

        # Make the digitized version of the E1C code (oversampling)
        galCodeE1CrepOversamp = galCodeE1Crep[codeOversampIdx]

        # ====================== Correlate signals =====================
        # Perform DFT of E1B code 
        galCodeE1BrepFreqDom = np.conj(fft(galCodeE1BrepOversamp))

        # Perform DFT of E1C code 
        galCodeE1CrepFreqDom = np.conj(fft(galCodeE1CrepOversamp))

        #--- Make the correlation for whole frequency band (for all freq. bins)
        for frqBinIndex in range(nFrqBins):
            #--- Generate carrier wave frequency grid (125 Hz step)
            frqBins[frqBinIndex] = settings.IF - (settings.acqFreqRangekHz/2)*1000 + 125*frqBinIndex

            #--- Generate local sine and cosine
            sinCarrReplica = np.sin(frqBins[frqBinIndex] * phasePoints)
            cosCarrReplica = np.cos(frqBins[frqBinIndex] * phasePoints)

            #--- Remove carrier from signal (demodulation)
            I1 = sinCarrReplica * signal1
            Q1 = cosCarrReplica * signal1
            I2 = sinCarrReplica * signal2
            Q2 = cosCarrReplica * signal2

            #--- Convert the signal to frequency domain 
            IQfreqDom1 = fft(I1+Q1*1j)
            IQfreqDom2 = fft(I2+Q2*1j)

            #--- Multiplication in the frequency domain (correlation in time domain) --> "Time loop"
            # E1B
            convCodeE1BIQ1 = IQfreqDom1 * galCodeE1BrepFreqDom
            convCodeE1BIQ2 = IQfreqDom2 * galCodeE1BrepFreqDom

            # E1C
            convCodeE1CIQ1 = IQfreqDom1 * galCodeE1CrepFreqDom
            convCodeE1CIQ2 = IQfreqDom2 * galCodeE1CrepFreqDom

            #--- Perform inverse DFT and store correlation results
            # E1B
            E1BacqRes1 = abs(ifft(convCodeE1BIQ1))**2
            E1BacqRes2 = abs(ifft(convCodeE1BIQ2))**2

            # E1C
            E1CacqRes1 = abs(ifft(convCodeE1CIQ1))**2
            E1CacqRes2 = abs(ifft(convCodeE1CIQ2))**2

            #--- Store results of non-coherent integration (integrate for two periods separately and sum)
            results[frqBinIndex] = E1BacqRes1 + E1BacqRes2 + E1CacqRes1 + E1CacqRes2

        # End of for frqBinIndex in range(nFrqBins)
        
        ## Look for correlation peaks in the results ==============================

        #--- Store result -----------------------------------------------------
        acqResults.searchSpace[PRN] = results
        
        #--- Find the correlation peak and the carrier frequency --------------
        peakSize = results.max()
        indices = np.where(results == peakSize)
        frequencyBinIndex = indices[0][0]
        # DBG
        # print("frequencyBinIndex", frequencyBinIndex, peakSize)
        # DBG

        delayBinIndex = indices[1][0]
        # DBG
        # print("delayBinIndex", delayBinIndex, peakSize)
        # DBG

        #--- Find 1 chip wide E1B delay to exclude chip around the peak ----
        excludeRangeIndex1 = delayBinIndex - samplesPerCodeChip
        excludeRangeIndex2 = delayBinIndex + samplesPerCodeChip

        #--- Correct E1B code phase exclude range if the range includes array
        #boundaries
        if excludeRangeIndex1 < 1:
            NoiseIdx = range(excludeRangeIndex2,
                            (samplesPerCode + excludeRangeIndex1 + 1))
                            
        elif excludeRangeIndex2 > samplesPerCode:
            NoiseIdx = range((excludeRangeIndex2 - samplesPerCode),
                            excludeRangeIndex1 + 1)
        else:
            NoiseIdx = list(range(excludeRangeIndex1)) + \
                            list(range(excludeRangeIndex2,samplesPerCode))
        
        # end of if excludeRangeIndex1 < 1

        # Estimate Noise level
        noiseLevel = np.mean(results[frequencyBinIndex, NoiseIdx])
        
        # Estimate Variance
        variance = np.std(results[frequencyBinIndex, NoiseIdx])

        #--- Store result -----------------------------------------------------
        acqResults.peakMetric[PRN] = (peakSize - noiseLevel)/variance if peakSize else 0
        
        #--- C/N0 -------------------------------------------------------------
        # Compute C/N0
        CNo = ((peakSize-noiseLevel)/noiseLevel/1e-3) if peakSize else 1
        acqResults.SNR[PRN] = (10 * np.log10(CNo)) if peakSize else 0

        # If the result is above threshold, then the satellite is acquired
        if (acqResults.peakMetric[PRN]) > settings.acqTh:
            ## Fine resolution frequency search =======================================
            
            #--- Indicate PRN number of the detected signal -------------------
            print(' %02d ' % PRN, end='')
            
            #--- Save properties of the detected satellite signal -------------
            acqResults.carrFreq[PRN]  = frqBins[frequencyBinIndex]
            acqResults.codeDelay[PRN] = delayBinIndex

        else:
            #--- No signal with this PRN --------------------------------------
            print(' _ ', end='')
        # end of if (acqResults.peakMetric[PRN]) > settings.acqTh

    # end of for PRN = satelliteList

    #=== Acquisition is over ==================================================
    print(')\n')

    #=== Plot results =========================================================
    acqResults.plot(settings, samplesPerCode, nFrqBins)

    return acqResults
