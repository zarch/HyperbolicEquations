# -*- coding: utf-8 -*-
"""
Define Flux functions

@author: pietro
"""
import numpy as np
import __future__

def toarray(xvect):
    if isinstance(xvect, np.ndarray):
        return xvect
    else:
        return np.array((xvect,))


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

from solve import newton

def nlers(qL, qR, xi):
    """
    function q = NLERS(qL,qR,xi)
        aL = eigenvalue(qL);
        aR = eigenvalue(qR);
        eta=0.001;
        if(aL > aR)
            % Shock wave
            s = (flux(qR)-flux(qL))/(qR-qL);
            if(xi<=s)
                q = qL;
            else
                q = qR;
            end
        else
            % Rarefaction
            if(xi<aL)
                q = qL;
            elseif(xi>aR)
                q = qR;
            else
                q=0.5*(qL+qR);
                q = Newton(q,xi);
            end
        end
    """
    aL = eigenvalue(qL);
    aR = eigenvalue(qR);
    if aL > aR:
        # is a Shock wave
        s = ( flux(qR) - flux(qL) ) / ( qR - qL );
        if (xi <= s):
            q = qL;
        else:
            q = qR;
    else:
        # is a Rarefaction
        if xi < aL:
            q = qL
        elif xi > aR:
            q = qR
        else:
            #print( 'is a Rarefaction')
            #import pdb; pdb.set_trace()
            q = 0.5 * (qL + qR)
            q = newton(q, xi)
    return q


def shock(xvect, x0, qL, qR):
    """
    >>> shock(range(10), 4, 2, 8)
    [[ 2.  2.  2.  2.  2.  8.  8.  8.  8.  8.]]
    """
    xvect = toarray(xvect)
    solution = np.empty(xvect.shape)
    solution[np.nonzero(xvect <= x0)] = qL
    solution[np.nonzero(xvect > x0)] = qR
    return solution


def rarefaction(xvect, xL, xR, qL, qR):
    xvect = toarray(xvect)
    solution = np.empty(xvect.shape)
    solution[np.nonzero(xvect < xL)] = qL
    solution[np.nonzero(xvect > xR)] = qR
    # find where: xL < xvect < xR
    rarefactionlist, = np.nonzero( (xvect > xL) * (xvect < xR) )
    return solution, rarefactionlist


def nlers2(qL, qR, xvect):
    aL = eigenvalue(qL)
    aR = eigenvalue(qR)
    xvect = toarray(xvect)

    if(aL > aR):
        # is a Shock wave
        s = ( flux(qR) - flux(qL) ) / ( qR - qL )
        solution = shock(xvect, s, qL, qR)
    else:
        # is a Rarefaction
        solution, rarefactionlist = rarefaction(xvect, aL, aR, qL, qR)
        for rarefact in rarefactionlist:
            solution[rarefact] = newton( 0.5 * ( qL + qR ), xvect[rarefact] )
    return solution


def godsca(qL, qR, dx, dt):
    """Godunov flux, pag 55, formula 4.14

    function flu = godSca(qL,qR)
        flu = flux(NLERS(qL,qR,0));
    """
    return flux( nlers(qL, qR, 0) ) # FIXME: perché 0? perché calcolo il flusso sul bordo tra due celle

def lfsca(qL, qR, dx, dt):
    """Lax-Friedrichs flux, pag 54, formula 4.8

    function flu = lfSca(qL,qR,dx,dt)
        flu = 1/2*(flux(qR)+flux(qL))-1/2*dx/dt*(qR-qL);
    """
    return 0.5 * ( flux(qR) + flux(qL) ) - 0.5 * dx/dt * ( qR - qL )

def lwsca(qL, qR, dx, dt):
    """Lax-Wendrof flux, pag 54, formula 4.9

    function flu = lwSca(qL,qR,dx,dt)
        a = 0.5 * ( eigenvalue(qL) + eigenvalue(qR) );
        flu = 1/2 * ( flux(qR)+flux(qL) ) - 1/2 * dt/dx * a^2 * ( qR - qL );
    """
    a = 0.5 * ( eigenvalue(qL) + eigenvalue(qR) )
    return 0.5 * ( flux(qR) + flux(qL) ) - 0.5 * dt/dx * a*a * ( qR - qL )

def forcesca(qL, qR, dx, dt):
    """FORCE flux, pag 54, formula 4.13

    function flu = forceSca(qL,qR,dx,dt)
        fLW = lwSca(qL,qR,dx,dt);
        fLF = lfSca(qL,qR,dx,dt);
        flu = 0.5*(fLW+fLF);
    """
    fLW = lwsca(qL, qR, dx, dt)
    fLF = lfsca(qL, qR, dx, dt)
    return 0.5 * ( fLW + fLF )

def roesca(qL, qR, dx, dt):
    """Roe flux, pag 63, formula 4.71

    function flu = forceSca(qL,qR,dx,dt)
        % Only for burgers
        a = 0.5*(qR+qL);
        flu = 0.5*(flux(qL)+flux(qR))-0.5*abs(a)*(qR-qL);
    """
    a = 0.5 * ( qR + qL )
    return 0.5 * ( flux(qL) + flux(qR) ) - 0.5 * abs(a) * ( qR - qL )

def oshermodsca(qL, qR, dx, dt): # :
    """Osher flux, pag 65, formula 4.83

    function flu = oshermodSca(qL,qR,dx,dt)
        aL = eigenvalue(qL);
        aR = eigenvalue(qR);
        flu = 0.5*(flux(qL)+flux(qR))-0.5*((abs(aL)+abs(aR))/2)*(qR-qL);
    """
    aL, aR = eigenvalue(qL), eigenvalue(qR);
    return 0.5 * ( flux(qL) + flux(qR) ) \
           - 0.5 * ( ( abs(aL) + abs(aR) ) * 0.5) * ( qR - qL )

#godSca, lfSca, lwSca, forceSca, roeSca, oshermodSca
NUM_MTD = {
    'godSca'       : godsca,
    'lfSca'        : lfsca,
    'lwSca'        : lwsca,
    'forceSca'     : forcesca,
    'roeSca'       : roesca,
    'oshermodSca'  : oshermodsca
}
