"""
matchFilt.py

Description: Define function Matched Filter function
"""

import numpy as np
from scipy import signal
import math
import simInit as si
from calcSNR import calcSNR
import sys

def matchFilt(rng, val):

    A = np.sqrt( 2 * calcSNR(rng,np.arctan2(500,rng)) )
    td = 2*rng/si.c0

    xmtChip = np.zeros(3)
    recChip= np.zeros(3)
    xmtPulse = np.zeros(si.bc.size-1)
    recPulse = np.zeros(si.bc.size-1)
    print('xmtPulse = ', xmtPulse)
    print('recPulse = ', recPulse)
    print('xmtChip = ', xmtChip)
    

    for ii in range(0, si.bc.size):
        t2 = 1
        xmtChip[ii] = si.bc[ii]
        recChip = xmtChip * np.exp(-1j*4*np.pi*si.pri*val*si.v*t2/si.lambda0)
        
        if ii == 0:
            xmtPulse[ii] = xmtChip
            recPulse = recChip
        else:
        
            print('ii = ', ii)
            print(si.bc.size)
            print('xmtChip = ', xmtChip)
            print('recChip = ', recChip)
            print('recChip.size = ', recChip.size)
            print('xmtPulse = ', xmtPulse)
            print('recPulse = ', recPulse)
            xmtPulse = np.concatenate((xmtPulse, xmtChip[:]), axis=0)
            recPulse = np.concatenate((recPulse, recChip[:,None]), axis=0)
            sys.exit()
            
    

    recPulse = np.array(A*recPulse, np.zeros(1,numSamps-len(recPulse)))
    numShifts = np.ceil(td/si.timeS)
    qNoise = np.random.randn(1, si.numSamps)
    iNoise = np.random.randn(1, si.numSamps)
    noiseSamp = np.sqrt(qnoise**2+inoise**2) * np.exp(-1j*np.arctan2(qnoise,inoise))

    recNoisySamp = recPulse + noiseSamp

    tempShiftSampNoisy = np.transpose(recNoisySamp)
    tempShiftSampNoisy = np.conjugate(tempShiftSampNoisy)
    tempRngSampNoisy = np.roll(tempShiftSampNoisy, numShifts)
    rngSampNoisy = np.transpose(tempRngSampNoisy)

    priCut = np.convolve(xmtPulse, np.fliplr(rngSampNoisy))

    return priCut
