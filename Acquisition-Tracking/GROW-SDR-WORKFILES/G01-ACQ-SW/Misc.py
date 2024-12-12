

########################################################################
# Misc.py:
# Miscelaneous functions
#
#  Project:        sw-rcvr-python
#  File:           Misc.py
#
#   Author: Nerea Sánchez / GNSS Academy
#   Copyright 2024 Nerea Sánchez / GNSS Academy
#
########################################################################


import math

def nextPowerOf2(x):
    return 1 if x == 0 else 2**math.ceil(math.log2(x))