#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHYS 432 W2022
Problem Set 3, Q3
Velocity field of lava flowing down an inclined slope
@author: Tara Brar
March 10, 2020
"""
import numpy as np
import matplotlib.pyplot as pl

# Set up the grid and advection and diffusion parameters
Ngrid = 50
Nsteps = 5000
dt = 0.001
dx = 1

# advection
v = 10 # [cm/s]
alpha = v*dt/2/dx

# diffusion
nu = 10**3 # [cm^2/s]
beta = nu*dt/dx**2

# gravity
angle = (np.pi)/6
g = 980 
gravTerm = g*np.sin(angle)*dt
H = 5

x = np.arange(0, Ngrid*1., dx) / Ngrid 
x = x*H

f1 = np.copy(x)

# no-slip condition, u(0) = 0
f1[0] = 0

# Set up plot
pl.ion()
fig, ax = pl.subplots(1,1)
ax.set_title('D=%3.1f')

# Plotting steady state in the background for reference
def steadyState(y, g, alpha, H, nu):
    return -g/nu * np.sin(alpha) * (1/2 * y**2 - H*y) 

ref = steadyState(x, g, angle, H, nu)
ax.plot(x, ref, 'k-')

plt1, = ax.plot(x, f1, 'ro')

# Setting the axis limits for visualization
ax.set_xlim([0,5])
ax.set_ylim([0,7])

# this draws the objects on the plot
fig.canvas.draw()

for ct in range(Nsteps):

    ## Calculate diffusion first
    # Setting up matrices for diffusion operator
    A1 = np.eye(Ngrid) * (1.0 + 2.0 * beta) + np.eye(Ngrid, k=1) * -beta + np.eye(Ngrid, k=-1) * -beta

    ## Boundary conditions to keep the first element fixed
    # This ensures f in the first cell stays fixed at all times under diffusion
    A1[0][0] = 1.0
    A1[0][1] = 0

    # Stress-free boundary condition at the lava-air interface
    A1[Ngrid - 1][Ngrid - 1] = 1.0 + beta
    
    # Solving for the next timestep
    f1 = np.linalg.solve(A1, f1)

    ## Calculate advection (Godunov)
    f1[1:Ngrid-1] = np.where(v > 0, f1[1:Ngrid-1] - 2*alpha*(f1[1:Ngrid-1]-f1[:Ngrid-2]), f1[1:Ngrid-1] - 2*alpha*(f1[2:]-f1[1:Ngrid-1]))

    # Accounting for gravity (but not correctly)
    #f1 += gravTerm
    
    # update the plot
    plt1.set_ydata(f1)

    fig.canvas.draw()
    pl.pause(0.001)

    ## Calculate diffusion first
    # Setting up matrices for diffusion operator
    A1 = np.eye(Ngrid) * (1.0 + 2.0 * beta) + np.eye(Ngrid, k=1) * -beta + np.eye(Ngrid, k=-1) * -beta

    ## Boundary conditions to keep the first element fixed
    # This ensures f in the first cell stays fixed at all times under diffusion
    A1[0][0] = 1.0
    A1[0][1] = 0

    # Stress-free boundary condition at the lava-air interface
    A1[Ngrid - 1][Ngrid - 1] = 1.0 + beta
    
    # Solving for the next timestep
    f1 = np.linalg.solve(A1, f1)

    ## Calculate advection (Godunov)
    f1[1:Ngrid-1] = np.where(v > 0, f1[1:Ngrid-1] - 2*alpha*(f1[1:Ngrid-1]-f1[:Ngrid-2]), f1[1:Ngrid-1] - 2*alpha*(f1[2:]-f1[1:Ngrid-1]))

    # Accounting for gravity (but not correctly)
    f1 += gravTerm
    
    # update the plot
    plt1.set_ydata(f1)

    fig.canvas.draw()
    pl.pause(0.00