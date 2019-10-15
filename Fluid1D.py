#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 16:45:56 2019


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
dx = 0.25
nCells = int((xL-x0)/dx) #Number of finite volumes (cells)
xCF = np.linspace(x0,xL,nCells+1) #Coordinates of cell faces
xCC = xCF[:-1]+dx*0.5 #Coordinates of cell centers
#parameters
D = 0.1
dt = 0.1
Eb = 1
xb = 31 #Position of the initial seed
Tfinal = 200

#Initialize fields
neCC = np.zeros(len(xCC)) #Electron density	
npCC = np.zeros(len(xCC)) #Positive ion density
ECC = np.zeros(len(xCC)) #Electric field at cell centers
ECF = np.zeros(len(xCF)) #Electric field at cell faces
#Initial conditions
neCC = 0.01*np.exp(-(xCC - xb)**2)
npCC = 0.01*np.exp(-(xCC - xb)**2)
ECF = -Eb*np.ones(len(xCF))
ECC = 0.5*(ECF[:-1]+ECF[1:]) 

t = 0.0
iiter =1
while t < Tfinal:
    #Computing the advective flux. East flux- flux on the right cell face & West flux- flux on the left cell face
    eastFlux = np.multiply(neCC,ECC) 
    westFlux = np.multiply(np.append(0.0, neCC[:-1]), np.append(0.0, ECC[:-1]))
    #Constructing the diffusion matrix to compute the diffusion flux
    diffusionMatrix = np.diag(np.ones(len(xCC)-1),1) - 2.0*np.diag(np.ones(len(xCC))) + np.diag(np.ones(len(xCC)-1),-1) #Gives a tridiagonal matrix
    diffusionMatrix[-1,-1] = -1.0 #Adjusting to take care of boundary conditions
    #Ghalf- flux terms
    Ghalf = (eastFlux - westFlux)/dx + (D/dx**2)*(diffusionMatrix.dot(neCC))
    #Shalf- source term
    Shalf = np.multiply(np.multiply(neCC, np.abs(ECC)), np.exp(-np.divide(np.ones(len(ECC)),np.abs(ECC))))
    #Computing the densities at half time step (t_n/2)
    neCChalf = neCC + dt*(Ghalf + Shalf)
    npCChalf = npCC + dt*(Shalf)
    ECFhalf = np.append(ECF[1:] + dx*(neCChalf - npCChalf),-Eb) #This is a simple numerical integration- usually we have to solve the Poisson's equation
    ECChalf = 0.5*(ECFhalf[:-1]+ECFhalf[1:])
    
    
    #Evaluating all the fluxes and source terms using the solution from half time step
    eastFlux = np.multiply(neCChalf, ECChalf)
    westFlux = np.multiply(np.append(0.0, neCChalf[:-1]), np.append(0.0, ECChalf[:-1]))
    G = (eastFlux - westFlux)/dx + (D/dx**2)*(diffusionMatrix.dot(neCChalf))
    S = np.multiply(np.multiply(neCChalf, np.abs(ECChalf)), np.exp(-np.divide(np.ones(len(ECChalf)),np.abs(ECChalf))))
    
    #Computing the solution at the new time step (t_n+1)
    neCC = neCC + 0.5*dt*(G + S + Ghalf + Shalf)
    npCC = npCC + 0.5*dt*(S+ Shalf)
    ECF = np.append(ECFhalf[1:] + dx*(neCC - npCC),-Eb)
    ECC = 0.5*(ECF[:-1]+ECF[1:])
    
    #Saving the solutions every 10 time step
    if (np.mod(iiter,10)) == 0:
        print('Saving')
        np.savetxt(str(iiter)+('.txt'),np.asarray([xCC, neCC,npCC, neCC-npCC, ECC]).T,delimiter=',')
    
    #Uncomment below to plot stuff real time
    #plt.plot(xCC,ECC)
    #plt.show()
    
    #This checks if the solution is diverging and stops the script
    if np.isnan(neCC).any():
        print('Solution diverges\n')
        break

    #Incrementing the time step and the iteration
    t = t + dt
    iiter = iiter + 1

print('Simulation done!')
    
