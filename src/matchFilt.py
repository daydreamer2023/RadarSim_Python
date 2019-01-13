"""
matchFilt.py

Description: Define function Matched Filter function
"""

# NEEDS CORRECTION

import numpy as np
from scipy import signal
import cmath
import math
import simInit as si
from calcSNR import calcSNR

def matchFilt(rng, val):

    A = math.sqrt( 2 * calcSNR(rng,math.atan2(500,rng)) )
    td = 2*rng/si.c0

    xmtPulse = np.zeros(si.bc.size-1)
    recPulse = np.zeros(si.bc.size-1)

    for ii in range(0, si.bc.size):
        t2 = 1
        xmtChip = si.bc[ii]
        recChip = xmtChip * cmath.exp(-1j*4*math.pi*si.pri*val*si.v*t2/si.lambda0)

        if ii == 0:
            xmtPulse[ii] = xmtChip
            recPulse[ii] = recChip
        else:
            xmtPulse = np.concatenate((xmtPulse, xmtChip), axis=1)
            recPulse = np.concatenate((recPulse, recChip), axis=1)

    recPulse = np.array(A*recPulse, np.zeros(1,numSamps-len(recPulse)))
    numShifts = math.ceil(td/si.timeS)
    qNoise = np.random.randn(1, si.numSamps)
    iNoise = np.random.randn(1, si.numSamps)
    noiseSamp = math.sqrt(qnoise**2+inoise**2) * cmath.exp(-1j*math.atan2(qnoise,inoise))

    recNoisySamp = recPulse + noiseSamp

    tempShiftSampNoisy = np.transpose(recNoisySamp)
    tempShiftSampNoisy = np.conjugate(tempShiftSampNoisy)
    tempRngSampNoisy = np.roll(tempShiftSampNoisy, numShifts)
    rngSampNoisy = np.transpose(tempRngSampNoisy)

    priCut = np.convolve(xmtPulse, np.fliplr(rngSampNoisy))

    return priCut
