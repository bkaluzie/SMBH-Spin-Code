# -*- coding: utf-8 -*-
"""
Author: Brooke Kaluzienski
Created: Fri Feb 16 13:33:26 2024

"""
import numpy as np
import matplotlib.pyplot as plt
#%% import tables
pgspins = np.loadtxt('pgspins.csv', delimiter = ',')
rgspins = np.loadtxt('rgspins.csv', delimiter = ',')

dj = 0.01
dm = 0.001
j0 = np.arange(0, 1 + dj, dj)
m0 = np.arange(0, 1 + dm, dm)
#%% bilinear interpolation 
def bilinear(j, M, m, spins):
    mfrac = m/M
    if int(j/dj) + 1 == len(j0):
        jind1 = int(j/dj) - 1
        jind2 = int(j/dj)
    else: 
        jind1 = int(j/dj)
        jind2 = int(j/dj) + 1
    if int(mfrac/dm) + 1 == len(m0):
        mind1 = int(mfrac/dm) - 1
        mind2 = int(mfrac/dm)
    else: 
        mind1 = int(mfrac/dm)
        mind2 = int(mfrac/dm) + 1
    
    j1 = j0[jind1]
    j2 = j0[jind2]
    m1 = m0[mind1]
    m2 = m0[mind2]
    
    N = 1/((j2 - j1) * (m2 - m1))
    F11 = spins[jind1][mind1]
    F12 = spins[jind1][mind2]
    F21 = spins[jind2][mind1]
    F22 = spins[jind2][mind2]
    A = F11 * (j2 - j) * (m2 - mfrac)
    B = F21 * (j - j1) * (m2 - mfrac)
    C = F12 * (j2 - j) * (mfrac - m1)
    D = F22 * (j - j1) * (mfrac - m1)
    F = N * (A + B + C + D)
    return F
#%%
spins = []
for n in range(1000):
    j = 0
    M0 = 10**3
    m = 10**3
    while M0 < 10**6:
        mfrac = m/M0
        p = np.random.randint(0, 2)
        if p == 0: # prograde accretion
            if float(int(mfrac/dm)) == mfrac/dm:
                j = pgspins[int(j/dj)][int(mfrac/dm)]
            elif float(int(mfrac/dm)) != mfrac/dm:
                j = bilinear(j, M0, m, pgspins)
        elif p != 0: # retrograde accretion
            if float(int(mfrac/dm)) == mfrac/dm:
                j = rgspins[int(j/dj)][int(mfrac/dm)]
            elif float(int(mfrac/dm)) != mfrac/dm:
                j = bilinear(j, M0, m, rgspins)
        M0 = M0 + m
    spins.append(np.round(j, 2))
#np.savetxt(f'data/spindata3in{int(np.log10(m))}-2.csv', spins)
#%%
plt.figure()
plt.hist(spins)