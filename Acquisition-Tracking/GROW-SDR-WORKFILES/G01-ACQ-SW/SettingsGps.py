
class Settings:
    # Input file containing Raw Signal
    inputf = "../DATA/GPS_RECORDED_RAW_SIGNAL_37000ms.dat"

    # signal
    signal = "gpsl1c"

    # Intermediate Frequency
    IF            = 9.548e6   # [Hz]

    # Sampling Frequency
    samplingFreq  = 38.192e6  # [Hz]

    # C/A code frequency
    codeFreqBasis = 1.023e6   # [Hz]

    # C/A code length [chips]
    codeLength    = 1023
    
    # GPS Mask
    satMask = range(1,32+1)

    # Number of Frequency bins in Acquisition
    acqFreqRangekHz = 14      # [kHz]
    acqTh           = 2.5

    # Number of channels of the Receiver
    numberOfChannels = 10

    # Number of milliseconds to be processed used 36000 + any transients (see
    # below - in Nav parameters) to ensure nav subframes are provided
    msToProcess        = 3000       #[ms]
    # msToProcess        = 36000       #[ms]

    ## Tracking loops settings ================================================
    # Code tracking loop parameters
    # dllDampingRatio         = 0.3
    # dllDampingRatio         = 2
    dllDampingRatio         = 0.7
    dllNoiseBandwidth       = 2       #[Hz]
    dllCorrelatorSpacing    = 1       #[chips]
    Nc                      = 0.001   #integration time in seconds 

    # Carrier tracking loop parameters
    pllDampingRatio         = 0.7
    # pllNoiseBandwidth       = 10      #[Hz]
    # pllNoiseBandwidth       = 60      #[Hz]
    pllNoiseBandwidth       = 25      #[Hz]

    ## Navigation filter settings =============================================
    zenithTravelTime  = 68.802

    navSolPeriod = 2

# end of class Settings
