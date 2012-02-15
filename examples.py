# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 17:04:58 2012

@author: pietro
"""

#
# Burgers
#
import num

# Numerical methods: godSca, lfSca, lwSca, forceSca, roeSca, oshermodSca
nummtd = 'godSca'
# Order of the methods: first, second
order  = 'second'
# find solution at the interface: simple, muscl
solintmtd = 'muscl'         # (only for the second order method)


burgers = num.Burgers(advection_speed = 1, x0 = -3, xL = -5, xR = 5,
        tSTART = 0, tEND = 0.8, iMAX = 200, nMAX = 1000, courant_numb=0.9,
        qL = 1., qR = 3., pause = 0.1)
burgers.compute_numerical('Lax Wendrof')
burgers.compute_exact()
burgers.show()
burgers.compute_numerical('Lax Friedrich')
burgers.show()
burgers.compute_numerical('force')
burgers.show()
burgers.compute_numerical('Godunov')
burgers.show()
burgers.compute_numerical('Roe')
burgers.show()
burgers.compute_numerical('Osher')
burgers.show()
#import pdb; pdb.set_trace()

#
# Systems
#



#
# SWE1D
#



#
# SWE2D
#

