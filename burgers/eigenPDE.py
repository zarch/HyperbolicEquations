# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 09:41:56 2012

@author: pietro
"""

from sympy import pprint
import sys
sys.displayhook = pprint
from sympy.matrices import *    # FIXME:
from sympy import *             # FIXME:
import matplotlib.pyplot as plt

x, y, z, t = symbols('x y z t')
a, b, l = symbols('a b l')
x1, x2 = symbols('x1 x2')

#A = Matrix([[0,a],[b,0]])
#(eval1, multipl1, basis1), (eval2, multipl2, basis2) = A.eigenvects()
#evect1,evect2 = basis1, basis2
#diagI = diag(1,1)
#P = A-l*diagI
#polyP = P.det()
#evals = solve(polyP, l)
#A.eigenvects()
#l1, l2 = A.eigenvals()
#C1 = A - diagI*l1
#C2 = A - diagI*l2
#unkn = Matrix([x1,x2])
#e1,e2 = C*unkn
#solve([e1, e2], set([x1, x2]))
#solve([e1, e2], [x1, x2])

def right_landa_left(matrix):
    eigenvals = []
    eigenvects = []
    multiplicities = []
    for val, mul, vect in matrix.eigenvects():
        eigenvals.append(val)
        multiplicities.append(mul)
        eigenvects.append(vect)
    landa = diag(*eigenvals)
    right = Matrix([[evect[0][i] for evect in eigenvects] for i in range(len(eigenvects))])
    left = right.inv()
    return right, landa, left
    
#A = Matrix([[0,a],[b,0]])
#R, V, L = right_landa_left(A)
#pprint(A)
#pprint(R*V*L)
#
#B = Matrix([[1,2,-2],[-1,1,1],[-1,2,0]])
#R, V, L = right_landa_left(B)
#pprint(B)
#pprint(R*V*L)

#A.berkowitz_eigenvals()
