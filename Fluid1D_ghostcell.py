#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 16:45:56 2019

@author: hemadity
"""


#Code using ghost cells
import numpy as np
import matplotlib.pyplot as plt

#Domain setup
x0 = 0.0
xL = 512.0
dx = 0.25
nCells = int((xL-x0)/dx)
xCF = np.linspace(x0,xL,nCells+1)
xCC = xCF[:-1]+dx*0.5
#parameters
D = 0.1
dt = 0.1
Eb = 1
xb = 31
Tfinal = 200

#Initialize fields
neCC = np.zeros(len(xCC))
npCC = np.zeros(len(xCC))
ECC = np.zeros(len(xCC))
ECF = np.zeros(len(xCF))
#Initial conditions
neCC = 0.01*np.exp(-(xCC - xb)**2)
npCC = 0.01*np.exp(-(xCC - xb)**2)
ECF = -Eb*np.ones(len(xCF))
ECC = 0.5*(ECF[:-1]+ECF[1:])

t = 0
iiter =1
while t < Tfinal:
    

    eastFlux = np.multiply(neCC,ECC) 
    westFlux = np.multiply(np.append(0.0, neCC[:-1]), np.append(0.0, ECC[:-1]))
    diffusionMatrix = np.diag(np.ones(len(xCC)-1),1) - 2.0*np.diag(np.ones(len(xCC))) + np.diag(np.ones(len(xCC)-1),-1)
    diffusionMatrix[-1,-1] = -1.0
    Ghalf = (eastFlux - westFlux)/dx + (D/dx**2)*(diffusionMatrix.dot(neCC))
    Shalf = np.multiply(np.multiply(neCC, np.abs(ECC)), np.exp(-np.divide(np.ones(len(ECC)),np.abs(ECC))))
    neCChalf = neCC + dt*(Ghalf + Shalf)
    npCChalf = npCC + dt*(Shalf)
    ECFhalf = np.append(ECF[1:] + dx*(neCChalf - npCChalf),-Eb)
    ECChalf = 0.5*(ECFhalf[:-1]+ECFhalf[1:])
    
    
    eastFlux = np.multiply(neCChalf, ECChalf)
    westFlux = np.multiply(np.append(0.0, neCChalf[:-1]), np.append(0.0, ECChalf[:-1]))
    #G = (eastFlux - westFlux)/dx + (D/dx**2)*(np.append(neCC[1:],neCC[-1]) - 2*neCC + np.append(0.0,neCC[:-1]))
    G = (eastFlux - westFlux)/dx + (D/dx**2)*(diffusionMatrix.dot(neCChalf))
    S = np.multiply(np.multiply(neCChalf, np.abs(ECChalf)), np.exp(-np.divide(np.ones(len(ECChalf)),np.abs(ECChalf))))

    neCC = neCC + 0.5*dt*(G + S + Ghalf + Shalf)
    npCC = npCC + 0.5*dt*(S+ Shalf)
    ECF = np.append(ECFhalf[1:] + dx*(neCC - npCC),-Eb)
    ECC = 0.5*(ECF[:-1]+ECF[1:])
    
    if (np.mod(iiter,10)) == 0:
        print('Saving')
        np.savetxt(str(iiter)+('.txt'),np.asarray([xCC, neCC,npCC, neCC-npCC, ECC]).T,delimiter=',')
    
    #plt.plot(xCC,ECC)
    #plt.show()
    
    if np.isnan(neCC).any():
        print('Solution diverges\n')
        break
    #dt = 0.3*dx/np.max(np.abs(ECC))
    t = t + dt
    iiter = iiter + 1
    #print('Time= ',t,'\n')
    