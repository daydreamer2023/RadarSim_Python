# simInit.py
# Description: Initialize simulation parameters

import numpy as np

# global constants
c0 = 3e8                                                # Speed of light (m/s)
k = 1.38e-23                                            # Boltzmann's Constant (J/K)

# waveform parameters
tn = 273                                                # Temperature (Kelvin)
freq = 20e9                                             # Center Frequency (Hz)
lambda0 = c0/freq                                       # Wavelength (m)
freqS = 26e6                                            # Sampling Frequency (Hz)
timeS = 1/freqS                                         # Sampling Time (s)
ts = 0.5                                                # Time Step (s)
prf = 50e3                                              # Pulse Repetition Frequency (Hz)
pri = 1/prf                                             # Pulse Repetition Interval (s)
numSamps = pri/timeS                                    # Number of samples (#)
numPulses = 128                                         # Number of pulses (#)
atten = 50                                              # Chebyshev window attenuation (dB)
bc = np.array([1, 1, 1, 1, 1, -1, -1,                   # Barker Code (Length = 13)
               1, 1, -1, 1, -1, 1])

# radar parameters
pt = 40e3                                               # Transmitted power (W)
tau = 1e-6                                              # Pulse width (s)
bw = 1/(tau/bc.size)                                    # Bandwidth (Hz)
RinitCh1 = np.array([10e3, 500])                        # Initialize radar channel 1
RinitCh2 = np.array([10e3, 500-0.0117])                 # Initialize radar channel 2

# target parameters
sigma = 0.1                                             # Target radar cross section (m**2)
v = np.array([250, 0])                                  # Target velocity vector (m/s)
