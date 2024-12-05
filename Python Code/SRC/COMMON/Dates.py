
import sys, os
from math import fmod

# Ref.: ESA GNSS Book TM-23 Vol I Section A.1.4 in Appendix A
def convertYearMonthDay2JulianDay(Year, Month, Day):
    # Case where month number is greater than 2
    if Month > 2:

        # Set year for algorithm
        NewYear = Year

        # Set month for algorithm
        NewMonth = Month

    # Case where month number is less than 2
    else:

        # Fix year for algorithm
        NewYear = Year - 1

        # Fix month for algorithm
        NewMonth = Month + 12

    # Compute A variable
    A = int(NewYear / 100)

    # Compute B variable
    B = 2 - A + int(A / 4)

    # Compute Julian date
    JulianDay = int(365.25 * NewYear) + \
                int(30.6001 * (NewMonth + 1)) + \
                Day + 1720981.5 + .5

    # Return Julian date
    return JulianDay

# Ref.: ESA GNSS Book TM-23 Vol I Section A.1.4 in Appendix A
def convertJulianDay2YearMonthDay(JulianDay):
    Jd2 = (JulianDay + 0.5)
    Z = int(Jd2)
    F = int(Jd2 - Z)
    Alpha = int((Z - 1867216.25) / 36524.25)
    A = ((Z + 1 + Alpha) - int((Alpha) / 4.0))
    B = (A + 1524)
    C = int(((B) - 122.1) / 365.25)
    D = int(365.25 * (C))
    E = int((B - D) / 30.6001)

    Day = ((B - D) - int(30.6001 * (E))) + F

    if ((E) < 13.5):
        Month = (E - 1)

    else:
        Month = (E - 13)

    if ((Month) > 2.5):
        Year = (C - 4716)

    else:
        Year = (C - 4715)

    return Year, Month, Day


def convertYearMonthDay2Doy(Year, Month, Day):
    # Do modulo 4 leap year check
    Modulo4Check = int(fmod(fmod((Year),4.)+4.,4.))

    # Do modulo 100 leap year check
    Modulo100Check = int(fmod(fmod((Year),100.)+100.,100.))

    # Do modulo 400 leap year check
    Modulo400Check = int(fmod(fmod((Year),400.)+400.,400.))

    # If year con't be divided by 4
    if (Modulo4Check != 0):
        # Current year isn't a leap year
        LeapYear = 0

    # Else: if year can't be divided by 100
    elif (Modulo100Check != 0):
        # Current year is a leap year
        LeapYear = 1

    # Else: if year can't be divided by 400
    elif (Modulo400Check != 0):
        # Current year isn't a leap year
        LeapYear = 0

    # Else:
    else:
        # Current year is a leap year
        LeapYear = 1

    # If current year is a leap year:
    if (LeapYear != 0):
        # Compute day of year using leap year formula
        DayOfYear = ((int((275 * Month) / 9.0) -\
                    int((Month + 9) / 12.0)) + Day) - 30
    # If current day isn't a leap year:
    else:
        # Compute day of year using non-leap year formula
        DayOfYear = ((int((275 * Month) / 9.0) - (2 *\
                    int((Month + 9) / 12.0))) + Day) - 30

    # Return output
    return DayOfYear


def convertJulianDay2EgnosEpoch(Jd):
    # Check that JD is int
    if not isinstance(Jd, int):
        sys.stderr.write("In convertJulianDay2EgnosEpoch: Jd not integer\n")
        sys.exit()

    InputYear, Month, Day = convertJulianDay2YearMonthDay(Jd)

    # If input year has two digits and is greater than 80:
    # assume that year is 19XX
    if ( (InputYear < 100) and (InputYear >= 80) ):
        # Correct year from YY to 19YY
        CorrectedYear = (InputYear + 1900)

    # If input year has two digits and is less than 80:
    # assume that year is 20XX
    elif ( InputYear < 80 ):
        # Correct year from YY to 20YY
        CorrectedYear = (InputYear + 2000)

    # If input year has more than two digits and is greater than 80:
    else:
        # Do not modify input year
        CorrectedYear = InputYear

    # Compute Julian Day
    CorrectedJd = convertYearMonthDay2JulianDay(CorrectedYear, Month, Day)

    # Compute EGNOS epoch
    EgnosEpoch = (CorrectedJd - 2444244.5 - (1024.0 * 7.0)) * 86400.0

    return EgnosEpoch

def convertYearDoy2JulianDay(Year, Doy, Sod):
    Jd = int( 365.25*(Year-1) ) + 428 + Doy + 1720981.5 + (Sod/86400)
    return Jd
