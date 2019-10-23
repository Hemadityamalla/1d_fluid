#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 16:26:14 2019

Description: Code to propagate a negative streamer with a constant background electric field.
             The equations and parameters used are all in their non-dimensional form
	     We use a finite volume discretization method as follows:
	     Advection term: Upwind discretization
	     Diffusion term: 2nd order central dscretization
             The time integration is performed using the trapezoidal rule (two stage)

@author: hemaditya
"""



import numpy as np
import matplotlib.pyplot as plt

#Domain setup
x0 = 0.0
xL = 512.0
dx = 0.0625
nCells = int((xL-x0)/dx) #Number of finite volumes (cells)
xCF = np.linspace(x0,xL,nCells+1) #Coordinates of cell faces
xCC = xCF[:-1]+dx*0.5 #Coordinates of cell centers
#parameters
D = 0.1
dt = 0.2*dx
Eb = 1
xb = 31 #Position of the initial seed
Tfinal = 200

#Initialize fields
neCC = np.zeros(len(xCC)+2) #Electron density	
npCC = np.zeros(len(xCC)+2) #Positive ion density
ECC = np.zeros(len(xCC)+2) #Electric field at cell centers
ECF = np.zeros(len(xCF)) #Electric field at cell faces

neCChalf = np.zeros(len(xCC)+2) #Electron density	
npCChalf = np.zeros(len(xCC)+2) #Positive ion density
ECChalf = np.zeros(len(xCC)+2) #Electric field at cell centers
ECFhalf = np.zeros(len(xCF)) #Electric field at cell faces
#Initial conditions
neCC[1:-1] = 0.01*np.exp(-(xCC - xb)**2)
npCC[1:-1] = 0.01*np.exp(-(xCC - xb)**2)
ECF = -Eb*np.ones(len(xCF))
ECC[1:-1] = 0.5*(ECF[:-1]+ECF[1:]) 

#Initialize ghost cells
neCC[0] = -neCC[1]
neCC[-1] = neCC[-2]
npCC[0] = -npCC[1]
npCC[-1] = npCC[-2]
#What about the electric field?


t = 0.0
iiter =1
while t < Tfinal:
    #Computing the advective flux. East flux- flux on the right cell face & West flux- flux on the left cell face
    eastFlux = np.multiply(neCC[1:-1],ECC[1:-1]) 
    westFlux = np.multiply(neCC[:-2], ECC[:-2])
    #Ghalf- flux terms
    Ghalf = (eastFlux - westFlux)/dx + (D/dx**2)*(neCC[2:] - 2.0*neCC[1:-1] + neCC[:-2])
    #Shalf- source term
    Shalf = np.multiply(np.multiply(neCC[1:-1], np.abs(ECC[1:-1])), np.exp(-np.divide(np.ones(len(ECC[1:-1])),np.abs(ECC[1:-1]))))
    #Computing the densities at half time step (t_n/2)
    neCChalf[1:-1] = neCC[1:-1] + dt*(Ghalf + Shalf)
    npCChalf[1:-1] = npCC[1:-1] + dt*(Shalf)
    ECFhalf = np.append(ECF[1:] + dx*(neCChalf[1:-1] - npCChalf[1:-1]),-Eb) #This is a simple numerical integration- usually we have to solve the Poisson's equation
    ECChalf[1:-1] = 0.5*(ECFhalf[:-1]+ECFhalf[1:])
    #Updating the ghost cells
    neCChalf[0] = -neCChalf[1]
    neCChalf[-1] = neCChalf[-2]
    npCChalf[0] = -npCChalf[1]
    npCChalf[-1] = npCChalf[-2]
    
    
    #Evaluating all the fluxes and source terms using the solution from half time step
    eastFlux = np.multiply(neCChalf[1:-1], ECChalf[1:-1])
    westFlux = np.multiply(neCChalf[:-2], ECChalf[:-2])
    G = (eastFlux - westFlux)/dx + (D/dx**2)*(neCChalf[2:] - 2.0*neCChalf[1:-1] + neCChalf[:-2])
    S = np.multiply(np.multiply(neCChalf[1:-1], np.abs(ECChalf[1:-1])), np.exp(-np.divide(np.ones(len(ECChalf[1:-1])),np.abs(ECChalf[1:-1]))))
    
    #Computing the solution at the new time step (t_n+1)
    neCC[1:-1] = neCC[1:-1] + 0.5*dt*(G + S + Ghalf + Shalf)
    npCC[1:-1] = npCC[1:-1] + 0.5*dt*(S+ Shalf)
    ECF = np.append(ECFhalf[1:] + dx*(neCC[1:-1] - npCC[1:-1]),-Eb)
    ECC[1:-1] = 0.5*(ECF[:-1]+ECF[1:])
    neCC[0] = -neCC[1]
    neCC[-1] = neCC[-2]
    npCC[0] = -npCC[1]
    npCC[-1] = npCC[-2]
    
    #Saving the solutions every 10 time step
    if (np.mod(iiter,10)) == 0:
        print('Saving')
        np.savetxt(str(iiter)+('.txt'),np.asarray([xCC, neCC[1:-1],npCC[1:-1], neCC[1:-1]-npCC[1:-1], ECC[1:-1]]).T,delimiter=',')
    
    #Uncomment below to plot stuff real time
    #plt.plot(xCC,ECC[1:-1])
    #plt.show()
    
    #This checks if the solution is diverging and stops the script
    if np.isnan(neCC).any():
        print('Solution diverges\n')
        break

    #Incrementing the time step and the iteration
    t = t + dt
    iiter = iiter + 1

print('Simulation done!')
    