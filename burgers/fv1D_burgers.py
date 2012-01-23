# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 15:10:33 2012

@author: pietro
"""
import burgers

# Define the parameters of the PDE of the mesh
a = 1            # advection speed
xL = -5          # position of the left boundary
xR = 5           # position of the right boundary
tSTART = 0       # initial time
tEND = 4         # final time
x0 = -0
iMAX = 100       # number of finite difference points

C = 0.9          # Courant number

# Boundary conditions B.C.
qL = 0.
qR = 1.

# equidistant distribution of iMAX points between xL and xR
xvect, dx = np.linspace(xL, xR, iMAX, retstep = True)

# Initial condition
q = burgers.shock(xvect, x0, qL, qR)

# equidistant distribution of points between each spatial step
xvect_b, dx_b = np.linspace(xL + dx / 2., xR - dx / 2., iMAX, retstep = True)
tvect, dt = np.linspace(tSTART, tEND, tMAX, retstep = True)
dx = (xR-xL)/(IMAX-1);  % mesh spacing
x = linspace(xL,xR,IMAX);  
xb = linspace(xL+dx/2,xR-dx/2,IMAX-1);

plot(xb,q,'o')
