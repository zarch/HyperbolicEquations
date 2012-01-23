# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 10:21:09 2012

@author: pietro
"""
import __future__
import fluxes as flx
import numpy as np

def g(q, xi):
    """
    function y = g(q,xi)
        y = eigenvalue(q)-xi;
    """
    return flx.eigenvalue(q)-xi


def dg(q, xi, eta=1e-7):
    """
    function y = dg(q,xi)
        eta=1e-7;
        y = (g(q+eta,xi)-g(q-eta,xi))/(2*eta);
    """
    return ( g(q + eta, xi) - g(q - eta, xi) ) / ( 2 * eta )


def newton(q0, xi, tol = 1e-12, iMAX = 100):
    """
    function q=Newton(q0,xi)
        tol=1e-12;
        q=q0;
        for i=1:100
            if(abs(g(q,xi))<tol)
                break
            end
            dq = g(q,xi)/dg(q,xi);
            q=q-dq;
        end
    """
    if isinstance(xi, np.ndarray):
        print('Error: "xi" must be a number')
        raise
    q = q0
    for i in range(iMAX):
        if(abs( g(q, xi) ) < tol):
            break
        dq = g(q, xi) / dg(q, xi)
        q = q-dq
    #if i == iMAX and abs( g(q, xi) ) < tol: print('Not converged using Newton')
    #print('Newton number of iteration to converge: {0:d}'.format(i))
    return q

