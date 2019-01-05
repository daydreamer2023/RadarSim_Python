# 725_main.py
#
# Project Path: C:\Users\Joseph\PycharmProjects\RadarSim
# Description: EE 725 Radar Project
#              (originally implemented using MATLAB, converted to Python)
#
#   Program Outline:
#       1) Initialize program parameters
#           a) Set up constants
#           b) Set up radar/waveform parameters
#       2) Generate target truth data for both radar channels
#       3) Generate range-pulse/range-Doppler matrices for both radar channels
#       4) Generate smoothed position with AB Filter
#       5) Generate smoothed elevation angle data with AB Filter
#       6) Plot results

# import libraries
import numpy as np
import math as m
import sys
from scipy import signal
import simInit as si
from calcSNR import calcSNR
from matchFilt import matchFilt
from alphaBetaFilt import alphaBetaFilt

# main program
def main():

    end = 40
    t = np.arange(0,end,si.ts)
    numel = t.size
    
    # initialize target truth data arrays for both radar channels
    trueRtCh1 = np.zeros((numel, 2))
    trueRtCh2 = np.zeros((numel, 2))
    trueLosRCh1 = np.zeros((numel, 2))
    trueLosRCh2 = np.zeros((numel, 2))
    trueElAngCh1 = np.zeros((numel, 2))
    trueElAngCh2 = np.zeros((numel, 2))
    trueLosRdotCh1 = np.zeros((numel, 2))
    trueLosRdotCh2 = np.zeros((numel, 2))
    fdCh1 = np.zeros((numel, 2))
    fdCh2 = np.zeros((numel, 2))
    PrCh1 = np.zeros((numel, 2))
    PrCh2 = np.zeros((numel, 2))

    # generate target truth data for both radar channels
    for ii in range(0,numel):
        trueRtCh1[ii,:] = np.subtract(si.RinitCh1, si.v*t[ii])
        trueRtCh2[ii,:] = np.subtract(si.RinitCh2, si.v*t[ii])

        trueLosRCh1[ii,:] = m.sqrt(trueRtCh1[ii,1]**2 + trueRtCh1[ii,2]**2)
        trueLosRCh2[ii,:] = m.sqrt(trueRtCh2[ii,1]**2 + trueRtCh2[ii,2]**2)

        trueElAngCh1[ii,:] = m.atan2(trueRtCh1[ii,2]/trueRtCh1[ii,1])
        trueElAngCh2[ii,:] = m.atan2(trueRtCh2[ii,2]/trueRtCh2[ii,1])

        trueLosRdotCh1[ii,:] = m.cos(trueElAngCh1[ii,1]) * si.v
        trueLosRdotCh2[ii,:] = m.cos(trueElAngCh2[ii,1]) * si.v

        fdCh1[ii,:] = 2 * trueLosRdotCh1[ii,:]/si.lambda0
        fdCh2[ii,:] = 2 * trueLosRdotCh2[ii,:]/si.lambda0

        PrCh1[ii,:] = calcSNR(trueLosRCh1[ii,1],trueElAngCh1[ii,1])
        PrCh2[ii,:] = calcSNR(trueLosRCh2[ii,1],trueElAngCh2[ii,1])

    # generate range-pulse/range-Doppler matrices for both radar channels
    for jj in range(0,numel):

        for kk in range(0,si.numPulses):
            priCutCh1[kk,:] = matchFilt(trueLosRCh1[jj,1],kk)
            priCutCh2[kk,:] = matchFilt(trueLosRCh2[jj,1],kk)

        priSize = max(priCutCh1.shape)

        for kk in range(0,priSize):
            priCh1Cheb = signal.chebwin(np.size(priCutCh1,kk),si.atten)
            priCh2Cheb = signal.chebwin(np.size(priCutCh2,kk),si.atten)

            priCutWinCh1[:,kk] = priCutCh1[:,kk]*priCh1Cheb
            priCutWinCh2[:,kk] = priCutCh2[:,kk]*priCh2Cheb

        rngPulseMatCh1 = priCutWinCh1
        rngPulseMatCh2 = priCutWinCh2

        priShape = priCutCh1.shape
        priIndex = priShape[1]

        for kk in range(0,priIndex):
            rngDopMatTransposeCh1[kk,:] = np.fft.fft(rngPulseMatCh1[:,kk])
            rngDopMatTransposeCh2[kk,:] = np.fft.fft(rngPulseMatCh2[:,kk])

        tempRngDopMatCh1 = np.transpose(rngDopMatTransposeCh1)
        tempRngDopMatCh2 = np.transpose(rngDopMatTransposeCh2)
        rngDopMatCh1 = np.conjugate(tempRngDopMatCh1)
        rngDopMatCh2 = np.conjugate(tempRngDopMatCh2)

        rngDopMatCh1dB = 10*np.log10(abs(rngDopMatCh1))
        rngDopMatCh2dB = 10*np.log10(abs(rngDopMatCh2))
        maxResponseCh1 = max(rngDopMatCh1dB)
        [indexY,indexX] = np.where(rngDopMatCh1dB == maxResponseCh1)

        measRng[jj] = (536 - indexX) * si.timeS * si.c0 / 2
        measDopFreq[jj] = indexY * (1/(si.pri*numPulses))
        measPhaseDiff[jj] = np.angle(rngDopMatCh1(indexY,indexX)-rngDopMatCh2(indexY,indexX))
        measElAng[jj] = m.asin(measPhaseDiff(index4)/(2*m.pi*0.78))

    # generate smoothed position data with alpha-beta filter
    smoothPos, smoothVel = alphaBetaFilt(0.5, 0.3, measRng)

    # generate smoothed elevation angle data with alpha-beta filter
    smoothElAng, smoothElAngRate = alphaBetaFilt(0.55, 0.3, measElAng)

    return None


if __name__ == "__main__":
    main()
