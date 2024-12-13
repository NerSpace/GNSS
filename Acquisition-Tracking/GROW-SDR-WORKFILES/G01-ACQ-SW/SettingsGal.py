
class Settings:
    # Input file containing Raw Signal
    inputf = "../DATA/GAL_RECORDED_RAW_SIGNAL_37000ms.dat"

    # signal
    signal = "gale1b"

    # Intermediate Frequency
    IF            = 6390000   # [Hz]

    # Sampling Frequency
    samplingFreq  = 26e6  # [Hz]

    # E1B code frequency
    codeFreqBasis = 1.023e6   # [Hz]

    # E1B code length [chips]
    codeLength    = 4092
    
    # GAL Mask
    satMask = range(1,36+1)

    # Number of Frequency bins in Acquisition
    acqFreqRangekHz = 12      # [kHz]
    acqTh           = 15

    # Number of channels of the Receiver
    numberOfChannels = 10

    # Number of milliseconds to be processed used 36000 + any transients (see
    # below - in Nav parameters) to ensure nav subframes are provided
    msToProcess        = 20        #[ms]


# end of class Settings
