# -*- coding: utf-8 -*-
"""
Configuration file of the PDE
"""
from burgers import plot_fv1D, plot_burgers1Dexact

#
# BURGERS
#
# Define the parameters of the PDE of the mesh
a = 1            # advection speed

# space domain
xL = -5          # position of the left boundary
xR = 5           # position of the right boundary

# time domain
tSTART = 0       # initial time
tEND = 2.5         # final time

x0 = -3

iMAX = 100      # number of finite difference points in space
tMAX = 100       # number of finite difference points in time

maxiter = 10000  # number of max iteration

courant = 0.9    # Courant number

# Boundary conditions B.C.
qL = 1.
qR = 3.

# OPTIONS
# Pause between each frame
pause=0.0005
# Numerical methods: godSca, lfSca, lwSca, forceSca, roeSca, oshermodSca
nummtd = 'godSca'
# Order of the methods: first, second
order  = 'second'
# find solution at the interface: simple, muscl
solintmtd = 'muscl'         # (only for the second order method)


#===========================================================================
#
# EXECUTE
#

# plot the exact solution
#plot_burgers1Dexact(xL, xR, iMAX, x0, qL, qR, tEND, tMAX, pause=0.005)

# compare the exact solution and the result of different numerical methods
plot_fv1D(xL, xR, x0, iMAX, qL, qR, tSTART, tEND, courant,
          nummtd, order, solintmtd, maxiter=10000, pause=0.005)