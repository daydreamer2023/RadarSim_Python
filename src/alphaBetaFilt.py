"""
alphaBetaFilt.py

Description: Define Alpha Beta Filter function
"""

import numpy as np
import simInit as si

def alphaBetaFilt(alpha, beta, inArr):

    smoothOut = np.zeros(inArr.size-2)
    smoothOutRate = np.zeros(inArr.size-2)

    for ii in range(1, inArr.size):

        if ii == 1:
            xHatNow = inArr[ii-1]
            xDotBarNowMinus = ( inArr[ii-1] - inArr[ii] ) / si.ts
        else:
            measXNow = inArr[ii]
            xBarNow = xHatNow + alpha * (measXNow - xHatNow)
            xDotBarNow = xDotBarNowMinus + beta * ( (measXNow - xHatNow)/si.ts )
            xHatNowPlus = xHatNow + xDotBarNow*si.ts
            xHatNow = xHatNowPlus

            smoothOut[ii-2] = xBarNow
            smoothOutRate[ii-2] = xDotBarNow

    return smoothOut, smoothOutRate
