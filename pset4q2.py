#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 22:35:29 2022

@author: taranjeetbrar
"""
"""
Phys 432: Problem Set 
@author: Tara Brar
March 24th 2022
"""
import numpy as np
import matplotlib.pyplot as pl

# Set up the grid, time and grid spacing, and the sound speed squared

Ngrid = 100
Nsteps = 5000
dt = 0.1
dx = 2.0
cs2 = 1.0 # sound speed squared


x = np.arange(Ngrid) * dx # grid
density = np.ones(Ngrid) # rho
mtmDensity= np.zeros(Ngrid) # rho x u
ergDensity= np.zeros(Ngrid) # rho x E

pressure = np.zeros(Ngrid) 
sndSpeed = np.zeros(Ngrid) 
M = np.zeros(Ngrid) 

u = np.zeros(Ngrid+1) # advective velocity (keep the 1st and last element zero)

# Initial conditions
Amp, sigma = 1, Ngrid/10 #perturbaton amplitude + width appropriate for strong shock (?)
ergDensity = ergDensity + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)

def advection(f, u, dt, dx):
    # calculating flux terms
    J = np.zeros(len(f)+1) # keeping the first and the last term zero
    J[1:-1] = np.where(u[1:-1] > 0, f[:-1] * u[1:-1], f[1:] * u[1:-1])
    f = f - (dt / dx) * (J[1:] - J[:-1]) #update

    return f

# plotting
pl.ion()
fig, (ax1,ax2) = pl.subplots(2,1)

x1, = ax.plot(x, f1, 'ro')

#ax.set_xlim([0, dx*Ngrid+1])
#ax.set_ylim([0.8, 2.1])

ax1.set_xlabel('')
ax1.set_ylabel('Density')


ax2.set_xlabel('$M$')
ax2.set_ylabel('Density')

for ct in range(Nsteps):
    # advection velocity at the cell interface
    u[1:-1] = 0.5 * ((mtmDensity[:-1] / density[:-1]) + (f2[1:] / f1[1:]))

    # update density and momentum
    density = advection(density, u, dt, dx)
    mtmDensity = advection(mtmDensity, u, dt, dx)

    # add the source term
    mtmDensity[1:-1] = mtmDensity[1:-1] - 0.5 * (dt / dx) * cs2 * (density[2:] - f1[:-2])

    # correct for source term at the boundary (reflective)
    mtmDensity[0] = mtmDensity[0] - 0.5 * (dt / dx) * cs2 * (density[1] - density[0])
    mtmDensity[-1] = mtmDensity[-1] - 0.5 * (dt / dx) * cs2 * (density[-1] - density[-2])

    # recalculating advection velocities
    u[1:-1] = 0.5 * ((mtmDensity[:-1] / density[:-1]) + (mtmDensity[1:] / density[1:]))

    # updating energy  
    ergDensity = advection(ergDensity, u, dt, dx)
    
    
    # pressure
    
    
    # update the plot
    x1.set_ydata(f1)
    fig.canvas.draw()
    pl.pause(0.001)