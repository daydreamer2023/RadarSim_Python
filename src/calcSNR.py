"""
calcSNR.py

Description: Define function to calculate SNR
"""

import math
import simInit as si

def calcSNR(rng, theta):

    # assume transmitGain = receiveGain
    gain = 316.2 * math.cos(math.radians(20 - theta))

    # calculate SNR (Pt/Pr)
    snr = si.pt * gain**2 * si.sigma * si.lambda0**2
    snr = snr / ((4*math.pi)**3 * rng**4 * si.k * si.tn * si.bw)

    return snr
