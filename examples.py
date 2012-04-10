# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 17:04:58 2012

@author: pietro
"""
import num
import math


#
# Burgers
#


def flux(q):
    """
    function y = flux(q)
        y = 0.5*q.^2;   % componentwise application
        %y = q^4/4;
    """
    return 0.5*q*q #q*q*q*q/4. # Example 21/23, pag. 31/32


def eigenvalue(q):
    """
    function y = eigenvalue(q)
        y = q;
        %y = q.^3;
    """
    return q #q*q*q # Example 21/23, pag. 31/32

#burgers = num.Burgers(advection_speed = 1, x0 = 0, xL = -1, xR = 5,
#        tSTART = 0, tEND = 0.8, iMAX = 200, nMAX = 1000, courant_numb=0.45,
#        qL = 3., qR = 1., pause = 0.01)
#
#burgers.eigenvalue = eigenvalue
#burgers.flux = flux
#
#burgers.compute_numerical('Lax Wendrof')
#burgers.compute_exact()
#burgers.show()
#burgers.compute_numerical('Lax Friedrich')
#burgers.show()
#burgers.compute_numerical('force')
#burgers.show()
#burgers.compute_numerical('Godunov')
#burgers.show()
#burgers.compute_numerical('Roe')
#burgers.show()
#burgers.compute_numerical('Osher')
#burgers.show()

#burgers.compute_numerical('Lax Wendrof', order='second')
#burgers.compute_exact()
#burgers.show()
#burgers.compute_numerical('Lax Friedrich', order='second')
#burgers.show()
#burgers.compute_numerical('force', order='second')
#burgers.compute_exact()
#burgers.show()
#burgers.compute_numerical('Godunov', order='second')
#burgers.show()
#burgers.compute_numerical('Roe', order='second')
#burgers.show()
#burgers.compute_numerical('Osher', order='second')
#burgers.show()

#burgers.compute_numerical('Lax Wendrof', order='second', type='muscl')
#burgers.compute_exact()
#burgers.show()
#burgers.compute_numerical('Lax Friedrich', order='second', type='muscl')
#burgers.compute_exact()
#burgers.show()
#burgers.compute_numerical('force', order='second', type='muscl')
#burgers.compute_exact()
#burgers.show()
#burgers.compute_numerical('Godunov', order='second', type='muscl')
#burgers.compute_exact()
#burgers.show()
#burgers.compute_numerical('Roe', order='second', type='muscl')
#burgers.compute_exact()
#burgers.show()
#burgers.compute_numerical('Osher', order='second', type='muscl')
#burgers.compute_exact()
#burgers.show()


#
# Systems
#
#A0 = [[0, 2],[ 1, 0]]
#R0 = [[-math.sqrt(2), math.sqrt(2)],[1., 1.]]
#D0 = [[-math.sqrt(2), 0],[0., math.sqrt(2)]]
#
#problem0 = num.Solver1Dsys(advection_speed = A0, R = R0, D = D0,
#                           qL = [2,1], qR = [1,0],
#                           iMAX=200, courant_numb=0.9)
#
#
#A1 = [[0, 1],[ 2, 0]]
#R1 = [[-math.sqrt(2), math.sqrt(2)], [1, 1]]
#D1 = [[-math.sqrt(2),0], [0, math.sqrt(2)]]
#
#problem1 = num.Solver1Dsys(advection_speed = A1, R = R1, D = D1,
#                           qL = [2,1], qR = [1,0],
#                           iMAX=200, courant_numb=0.9)
#
#
#problem0.compute_numerical('Lax Wendrof')
#problem0.show2()
#problem1.compute_numerical('Lax Friedrich')
#problem1.show2()
#problem.compute_numerical('Warming Beam')
#problem.show2()
#problem.compute_numerical('Fromm')
#problem.show2()
#problem.compute_numerical('WENO')
#problem.show2()
#problem.compute_numerical('DGRK')
#problem.show2()

dgrk = num.DGRK(iMAX = 20, accuracy = 2,
            M = [(1, 0, 0),( 0, 1./3., 0),(0, 0, 1./5.)],
            K = [(0, 0, 0),(2, 0, 0),(0, 2, 0)])
dgrk.compute_numerical('DGRK')
dgrk.show2()

import swe
#
# SWE1D
#
problem = swe.Swe1D(iMAX=200, nMAX=1000, courant_numb=0.95, )
problem.compute_numerical('Godunov')
problem.show2()


#
# SWE2D
#
#problem = Swe2D(xL = 0, xR = 1, yL = 0, yR = 1, tEND=0.2,
#                iMAX = 20, jMAX = 20, nMAX=10000,)
#                #tri = tri4, neigh = triNeigh4)
#problem.compute_numerical('Godunov')
#problem.exp2vtk('provaVTKmatlab')
