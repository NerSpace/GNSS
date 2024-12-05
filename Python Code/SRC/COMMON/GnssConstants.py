#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PHYSICAL CONSTANTS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Speed of light (meters/sec)
SPEED_OF_LIGHT=299792458.0 

# Earth's radius (meters) 
EARTH_RADIUS=6378136.3

# Semi mayor axis of the earth (meters) 
EARTH_SEMIAXIS=6378137.0

# Earth Flattening 
FLATTENING=1.0/298.257223563

# Constant E2 (FLATTENING*(2.0 - FLATTENING))
E2=0.0066943799901

# Constant E12=(E2/(1.0 - E2))
E12=0.0067394967423

# Earth's rotation rate (rad/sec)
OMEGA_EARTH=7.2921151467e-5

# Earth's gravitational 
MU_EARTH=3.986004415e+14

# Earth's eccentricity
ECCENTRICITY_EARTH=8.2e-2

# Estimated latitude  of the Magnetic Dipole [deg] 
MAGNETIC_DIPOLE_LATITUDE=80.37

# Estimated longitude  of the Magnetic Dipole [deg] 
MAGNETIC_DIPOLE_LONGITUDE=287.37

# GPS L1 frequency (MHz)
GPS_L1_MHZ=1575.42

# GPS L2 frequency (MHz) 
GPS_L2_MHZ=1227.60

# GPS L5 frequency (MHz) 
GPS_L5_MHZ=1176.45

# GPS L1 wave length (meters) = SPEED_OF_LIGHT/(GPS_L1_MHZ*1e6) 
GPS_L1_WAVE=0.19029367279836487

# GPS L2 wave length (meters) = SPEED_OF_LIGHT/(GPS_L2_MHZ*1e6) 
GPS_L2_WAVE=0.24421021342456825

# GPS L5 wave length (meters) = SPEED_OF_LIGHT/(GPS_L5_MHZ*1e6) 
GPS_L5_WAVE=0.254828049

# Gamma value as GPS_L1_MHZ^2/GPS_L2_MHZ^2 
GPS_GAMMA_L1L2=1.646944444444444

# Gamma value as GPS_L1_MHZ^2/GPS_L5_MHZ^2 
GPS_GAMMA_L1L5=1.793270321

# Gamma value as GPS_L2_MHZ^2/GPS_L5_MHZ^2 
GPS_GAMMA_L2L5=1.088846881

# GPS K1 value = GAMMA/(1-GAMMA) 
GPS_K1_L1L2=-2.545727780163161

# GPS K2 value = 1/(1-GAMMA) 
GPS_K2_L1L2=-1.545727780163161

# GPS K1 value = GAMMA/(1-GAMMA) 
GPS_K1_L1L5=-2.260604328

# GPS K2 value = 1/(1-GAMMA) 
GPS_K2_L1L5=-1.260604328

# Galileo E1 frequency (MHz)
GAL_E1_MHZ     = 1575.42

# Galileo E1 wave length (meters)   SPEED_OF_LIGHT/(GAL_E1_MHZ*1e6)
GAL_E1_WAVE    = 0.19029367279836487

# Galileo E5A frequency (MHz)
GAL_E5A_MHZ    = 1176.45

# Galileo E5B frequency (MHz)
GAL_E5B_MHZ    = 1206.14

# Galileo E5A wave length (meters)   SPEED_OF_LIGHT/(GAL_E5A_MHZ*1e6)
GAL_E5A_WAVE   = 0.25482804879085386

# Galileo E5B wave length (meters)   SPEED_OF_LIGHT/(GAL_E5B_MHZ*1e6)
GAL_E5B_WAVE   = 0.24855527384880693

# Gamma value as GAL_E1_MHZ**2/GAL_E5A_MHZ**2
GAL_GAMMA_E1E5A = 1.7932703213610588

# Gamma value as GAL_E1_MHZ**2/GAL_E5B_MHZ**2
GAL_GAMMA_E1E5B = 1.7060716571123251

# Galileo K1 value   GAMMA/(1-GAMMA)
GAL_K1_E1E5A    = -2.262277204112014

# Galileo K1 value   GAMMA/(1-GAMMA)
GAL_K1_E1E5B    = -2.420010419144861

# Galileo K2 value   1/(1-GAMMA)
GAL_K2_E1E5A    = -1.262277204112014

# Galileo K2 value   1/(1-GAMMA)
GAL_K2_E1E5B    = -1.4200104191448613

# Convert TECU to Meters in L1
TEC_TO_METERS_L1=0.1624

# Convert TECU to Meters in L2
TEC_TO_METERS_L2=0.2674

# Convert TECU to Meters in L5
TEC_TO_METERS_L5=0.2911

# Height of the Iono layer (meters)
IONO_HEIGHT=350000.0

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TIME CONSTANTS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Seconds in one hour
S_IN_H = 3600.0

# Seconds in one hour
S_IN_D = 86400

# Days in one week
D_IN_W = 7

# Julian date for GPS start epoch (1980 January 6)
JD_0 = 2444244.5

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# IONO CONSTANTS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Iono regions latitude upper limit [DEG]
MIDHI_LAT = 61.0
MID_LAT = 50.0
MIDLO_LAT = 32.0
EQUAT_LAT = 20.0

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# OTHER CONSTANTS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Not a Number
NAN = -9999.9

# Minimum number of satellites to compute PVT solution
MIN_NUM_SATS_PVT = 4

# Maximum PDOP to compute PVT solution
MAX_PDOP_PVT = 10000

# Minimum norm of the Receiver Position Delta
# (cut condition of LSQ)
LSQ_DELTA_EPS = 1e-4

# MOPS Kh factor in PA
MOPS_KH_PA = 6.0

# MOPS Kv factor in PA
MOPS_KV_PA = 5.33

# APV-I Horizontal Alarm Limit (HAL) req [m]
APVI_HAL = 40

# APV-I Vertical Alarm Limit (VAL) req [m]
APVI_VAL = 50

# Minimum of longitudes [deg]
MIN_LON = -180.0

# Maximum of longitudes [deg]
MAX_LON = 180.0

# Minimum of latitudes [deg]
MIN_LAT = -90.0

# Maximum of latitudes [deg]
MAX_LAT = 90.0

# Minimum Masking angle [deg]
MIN_MASK_ANGLE = 0

# Maximum Masking angle [deg]
MAX_MASK_ANGLE = 90.0

# Maximum number of satellites per constellation
MAX_NUM_SATS_CONSTEL = 36

# Minimum PRN of a GEO
MIN_GEO_PRN = 120

# Maximum PRN of a GEO
MAX_GEO_PRN = 158

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DIMENSIONING CONSTANTS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Maximum number of IGPs in IGP Mask (dimensioning constant)
MAX_NUM_IGPS = 500

# Maximum number of Users in User Grid (dimensioning constant)
MAX_NUM_USRS = 1000

# Maximum number of Receivers in RCVR (dimensioning constant)
MAX_NUM_RCVR = 1000
