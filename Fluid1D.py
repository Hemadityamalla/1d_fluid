#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 17:03:48 2019

@author: hemadity
"""


import numpy as np
import matplotlib.pyplot as plt

#Domain setup
x0 = 0.0
xL = 1024.0
dx = 0.25
nCells = int((xL-x0)/dx)
xCF = np.linspace(x0,xL,nCells+1)
xCC = xCF[:-1]+dx*0.5
#parameters
D = 0.1
dt = 0.1
Eb = 1
xb = 512
Tfinal = 200

#Initialize fields
neCC = np.zeros(len(xCC))
ECC = np.zeros(len(xCC))
ECF = np.zeros(len(xCF))
#Initial conditions
neCC = 0.01*np.exp(-0.01*(xCC - xb)**2)
npCC = 0.01*np.exp(-0.01*(xCC - xb)**2)
ECF = -Eb*np.ones(len(xCF))
ECC = 0.5*(ECF[:-1]+ECF[1:])

t = 0
while t < Tfinal:
    

    eastFlux = np.multiply(np.append(neCC[1:],neCC[-1]),np.append(ECC[1:],-Eb))
    westFlux = np.multiply(neCC,ECC)
    diffusionMatrix = -2.0*np.diag(np.ones(len(neCC)))
    diffusionMatrix[0,0] = -3.0
    diffusionMatrix[-1,-1] = 0
    diffusionMatrix[-1,-2] = -1.0
    Ghalf = (eastFlux - westFlux)/dx + (D/dx**2)*(np.append(neCC[1:],neCC[-1]) + np.diag(diffusionMatrix*neCC) + np.append(np.append(0.0,neCC[1:-1]),0.0))
    Shalf = np.multiply(np.multiply(neCC, np.abs(ECC)), np.exp(-np.divide(np.ones(len(ECC)),np.abs(ECC))))
    neCChalf = neCC + dt*(Ghalf + Shalf)
    npCChalf = npCC + dt*(Shalf)
    ECFhalf = np.append(ECF[1:] + dx*(neCChalf - npCChalf),-Eb)
    ECChalf = 0.5*(ECF[:-1]+ECF[1:])
    
    
    eastFlux = np.multiply(np.append(neCChalf[1:],neCChalf[-1]),np.append(ECChalf[1:],-Eb))
    westFlux = np.multiply(neCChalf,ECChalf)
    #G = (eastFlux - westFlux)/dx + (D/dx**2)*(np.append(neCC[1:],neCC[-1]) - 2*neCC + np.append(0.0,neCC[:-1]))
    G = (eastFlux - westFlux)/dx + (D/dx**2)*(np.append(neCChalf[1:],neCChalf[-1]) + np.diag(diffusionMatrix*neCChalf) + np.append(np.append(0.0,neCChalf[1:-1]),0.0))
    S = np.multiply(np.multiply(neCChalf, np.abs(ECChalf)), np.exp(-np.divide(np.ones(len(ECChalf)),np.abs(ECChalf))))

    neCC = neCC + 0.5*dt*(G + S + Ghalf + Shalf)
    npCC = npCC + 0.5*dt*(S+ Shalf)
    ECF = np.append(ECF[1:] + dx*(neCC - npCC),-Eb)
    ECC = 0.5*(ECF[:-1]+ECF[1:])
    plt.plot(xCC,ECC)
    plt.show()
    
    if np.isnan(neCC).any():
        print('Solution diverges\n')
        break
    dt = 0.1*dx/np.max(np.abs(ECC))
    t = t + dt
    print('Time= ',t,'\n')
    