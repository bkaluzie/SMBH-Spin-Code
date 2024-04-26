# -*- coding: utf-8 -*-
"""
Author: Brooke Kaluzienski
Created: Fri Dec  1 10:02:12 2023

"""
import numpy as np
#%% equations 
def radius_p(j): # note: missing mass factor M
    Z1 = 1 + ((1-j**2)**(1/3) * ((1+j)**(1/3) + (1-j)**(1/3)))
    Z2 = np.sqrt(3*j**2 + Z1**2)
    r = 3 + Z2 - np.sqrt((3 - Z1)*(3 + Z1 + 2*Z2))
    return r

def momentum_p(j): # note: missing mass factor M (carries over from radius)
    r = radius_p(j)
    if j == 1: 
        return 1.16 # 2E
    else:   
        u = (1/np.sqrt(r)) * (r**2 - 2*j*np.sqrt(r) + j**2) / np.sqrt(r**2 - 3*r + 2*j*np.sqrt(r))
        return u

def energy_p(j): # note: INDEPENDENT of M
    r = radius_p(j)
    if j == 1:
        return 0.58 # energy not radiated away
    else:
        u = (1/r) * (r**2 - 2*r + j*np.sqrt(r)) / np.sqrt(r**2 - 3*r + 2*j*np.sqrt(r))
        return u
    
def pg(Mfrac, j0, M0, J0, d):
    J = J0 + ((Mfrac/d) * M0**2 * momentum_p(j0))
    M = M0 + ((Mfrac/d) * M0 * energy_p(j0))
    return M, J

def pg_final(m, j0, M0, J0):
    J = J0 + (m * momentum_p(j0))
    M = M0 + (m/M0 * energy_p(j0))
    return M, J

def radius_r(j):
    Z1 = 1 + (1-j**2)**(1/3)*((1+j)**(1/3) + (1-j)**(1/3))
    Z2 = (3*j**2 + Z1**2)**(1/2)
    r = (3 + Z2 + ((3 - Z1)*(3 + Z1 + 2*Z2))**(1/2))
    return r

def momentum_r(j):
    r = radius_r(j)
    u = -(np.sqrt(r) * (r**2 + 2*j*np.sqrt(r) + j**2)) / (r * np.sqrt(r**2 - 3*r - 2*j*np.sqrt(r)))
    return u

def energy_r(j):
    r = radius_r(j)
    u = (r**2 - 2*r - j*np.sqrt(r)) / (r * np.sqrt(r**2 - 3*r - 2*j*np.sqrt(r)))
    return u

def rg(Mfrac, j0, M0, J0, d):
    J = J0 + ((Mfrac/d) * M0**2 * momentum_r(j0))
    M = M0 + ((Mfrac/d) * M0 * energy_r(j0))
    return M, J

def rg_final(m, j0, M0, J0):
    J = J0 + (m * momentum_r(j0))
    M = M0 + (m/M0 * energy_r(j0))
    return M, J
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
#%% generating the tables
dj = 0.01
dm = 0.001
j0 = np.arange(0, 1 + dj, dj)
m0 = np.arange(0, 1 + dm, dm)
mfrac = np.arange(0, 1 + dm, dm)

pgspins = np.zeros(len(mfrac)*len(j0)).reshape(len(j0), len(mfrac))
for i in range(len(j0)):
    pgspins[i][0] = j0[i]
    for k in range(1, len(mfrac)):
        j = j0[i]
        M = 1
        J = j * M**2
        while M < (1 + mfrac[k]):
            M_prev = M # need for last little bit
            J_prev = J
            j = J / M**2
            M, J = pg(mfrac[k], j, M, J, 500) # mass after NEXT accretion
        # add last bit of accretion here
        m = (1 + mfrac[k] - M_prev) / energy_p(j)
        M, J = pg_final(m, j, M_prev, J_prev)
        j = J / M**2
        pgspins[i][k] = np.round(j, 2)
#np.savetxt('pgspins.csv', pgspins, delimiter = ',')

rgspins = np.zeros(len(mfrac)*len(j0)).reshape(len(j0), len(mfrac))
for i in range(len(j0)):
    rgspins[i][0] = j0[i]
    for k in range(1, len(mfrac)):
        j = j0[i]
        M = 1
        J = j * M**2
        ret = True
        while M < (1 + mfrac[k]):
            M_prev = M # need for last little bit
            J_prev = J
            j = J / M**2
            if j != abs(j): ret = False
            if ret == True:
                M, J = rg(mfrac[k], j, M, J, 500)
            elif ret == False:
                M, J = pg(mfrac[k], j, M, J, 500)
        if ret == True:
            m = (1 + mfrac[k] - M_prev) / energy_r(j)
            M, J = rg_final(m, j, M_prev, J_prev)
            j = J / M**2
        elif ret == False:
            m = (1 + mfrac[k] - M_prev) / energy_p(j)
            M, J = pg_final(m, j, M_prev, J_prev)
            j = J / M**2
        rgspins[i][k] = np.round(j, 2)
#np.savetxt('rgspins.csv', rgspins, delimiter = ',')
#%% testing accuracy
m = (10**5)

j_tabs = []
j_calcs = []
for n in range(1000):
    j_tab = 0
    j_calc = 0
    M0 = 10**5
    while M0 < 10**6:  
        mfrac = m/M0
        M = M0
        j = j_calc
        J = j * M0**2
        
        p = np.random.randint(0, 2)
        if p == 0: # prograde accretion
            while M < M0 * (1 + mfrac):
                M_prev = M # need for last little bit
                J_prev = J
                j = J / M**2
                M, J = pg(mfrac, j, M, J, 500)
            mf = (M0 + m - M_prev) / energy_p(j)
            M, J = pg_final(mf, j, M_prev, J_prev)
            j = J / M**2
            j_calc = np.round(j, 2)
        
            if float(int(mfrac/dm)) == mfrac/dm:
                j_tab = pgspins[int(j_tab/dj)][int(mfrac/dm)]
            elif float(int(mfrac/dm)) != mfrac/dm:
                j_tab = bilinear(j_tab, M0, m, pgspins)
                
        elif p != 0: # retrograde accretion
            ret = True
            while M < M0 * (1 + mfrac):
                M_prev = M # need for last little bit
                J_prev = J
                j = J / M**2
                if j != abs(j): ret = False
                if ret == True:
                    M, J = rg(mfrac, j, M, J, 500)
                elif ret == False:
                    M, J = pg(mfrac, j, M, J, 500)
            if j != abs(j): ret = False
            if ret == True:
                mf = (M0 + m - M_prev) / energy_r(j)
                M, J = rg_final(mf, j, M_prev, J_prev)
                j = J / M**2
            elif ret == False:
                mf = (M0 + m - M_prev) / energy_p(j)
                M, J = pg_final(mf, j, M_prev, J_prev)
                j = J / M**2 
            j_calc = np.round(j, 2)
            
            if float(int(mfrac/dm)) == mfrac/dm:
                j_tab = rgspins[int(j_tab/dj)][int(mfrac/dm)]
            elif float(int(mfrac/dm)) != mfrac/dm:
                j_tab = bilinear(j_tab, M0, m, rgspins)
        M0 = M0 + m
        J = j_calc * M0**2
    j_tabs.append(np.round(j_tab, 2))
    j_calcs.append(j_calc)

diff = np.array(j_tabs) - np.array(j_calcs)