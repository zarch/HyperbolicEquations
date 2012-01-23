# -*- coding: utf-8 -*-
"""
Configuration file of the PDE
"""
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

iMAX = 1000      # number of finite difference points in space
tMAX = 100       # number of finite difference points in time

maxiter = 10000  # number of max iteration

courant = 0.9    # Courant number

# Boundary conditions B.C.
qL = 3.
qR = 1.

# Pause between each frame
pause=0.0005