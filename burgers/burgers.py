# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 16:54:14 2012

@author: Pietro Zambelli

How to use
===========

The program is divided in subfiles:

    * `confPDE` to set the PDE mesh varaibles;
    * `fluxes` contain all the function needed to compute fluxes using
      different methods;
    * `solve` contain Newton method to converge to the solution;
    * `burgers` contain method to run and resolve the PDE and plotting the
      result.

"""
import numpy as np
import matplotlib.pyplot as plt
import fluxes as flx
import __future__

def prettyp(lst):
    return ' | '.join([ '{0:+3.2f}'.format(l) for l in lst])

def minmod(a,b):
    """Minmod, formula. 4.104, pag. 67

    function slope = minmod(a,b)
        if(a*b<0)
            slope=0;
        elseif(abs(a)<=abs(b))
            slope=a;
        elseif(abs(a)>abs(b))
            slope=b;
        end
    """
    if a * b < 0:
        slope = 0
    elif abs(a) <= abs(b):
        slope = a
    elif abs(a) > abs(b):
        slope = b
    return slope

def get_sol_interfaceMUSCL(sol, dx, dt):
    """
    for i=2:IMAX-2
        slope = minmod(q(i)-q(i-1),q(i+1)-q(i));
        qp0 = q(i)-0.5*slope;
        qm0 = q(i)+0.5*slope;
        qp(i) = qp0 + 1/2 * dt/dx * ( flux(qp0) - flux(qm0) );
        qm(i) = qm0 + 1/2 * dt/dx * ( flux(qp0) - flux(qm0) );
    end
    qm(1)=q(1);
    qp(IMAX)=q(IMAX-1);
    """
    slope, sol_p, sol_m = [],[],[]
    for i in range(2, len(sol) - 2):
        slope = minmod(sol[i] - sol[i-1], sol[i+1] - sol[i])
        sol0_p = sol[i] - 0.5 * slope[i]
        sol0_m = sol[i] + 0.5 * slope[i]
        sol_p.append( sol0_p + 0.5* dt/dx * (flx.flux(sol0_p) - flx.flux(sol0_m) ) )
        sol_m.append( sol0_m + 0.5* dt/dx * (flx.flux(sol0_p) - flx.flux(sol0_m) ) )
    # set values at the boundaries
    sol_p[-1] = sol[-1]
    sol_m[ 0] = sol[ 0]
    return sol_p, sol_m


def get_sol_interface(sol, dx, dt):
    """
    for i=2:IMAX-2
        slope = minmod(q(i)-q(i-1),q(i+1)-q(i));
        qp(i) = q(i)-0.5*slope;
        qm(i) = q(i)+0.5*slope;
    end
    qm(1)=q(1);
    qp(IMAX)=q(IMAX-1)
    """
    slope, sol_p, sol_m = np.empty(sol.shape), np.empty(sol.shape), np.empty(sol.shape)
    for i in range(1, len(sol) - 1):
        #import pdb; pdb.set_trace()
        #print i, sol[i-1], sol[i], sol[i+1]
        slope[i] = minmod(sol[i] - sol[i-1], sol[i+1] - sol[i])
        sol_p[i] = sol[i] - 0.5 * slope[i]
        sol_m[i] = sol[i] + 0.5 * slope[i]
    # set values at the boundaries
    sol_p[-1] = sol[-1]
    sol_m[ 0] = sol[ 0]
    print('slope: {}'.format(prettyp(slope[:10])))
    print('sol_p: {}'.format(prettyp(sol_p[:10])))
    print('sol_m: {}'.format(prettyp(sol_m[:10])))
    return sol_p, sol_m



def checkCFL(solution, dx, courant, time, tEND,):
    """Check CFL stability and final time matching

    dt = c * dx / max( abs( eigenvalue(q) ) );
    if(time + dt > t1)
        dt = t1-time;
    end
    if(time >= t1)
        break
    end"""
    maxflux = np.max( np.abs( flx.eigenvalue(solution) ) )
    dt = courant * dx / maxflux
    if dt == np.nan: return -9999
    if time + dt > tEND:
        return tEND - time
    if time >= tEND:
        return -9999
    strformat = 'time={t:6.5f}; dt={dt:6.5f}; dx={dx:6.5f}; c={c:3.2f}; maxflux={fl:6.5f};'
    print(strformat.format(t=time, dt=dt, dx=dx, c=courant, fl=maxflux) )
    return dt


#===========================================================================
def burgers1Dexact(xvect, tvect, x0, qL, qR, iMAX):
    for ti in tvect:
#        yield flx.nlers2(qL, qR, (xvect-x0)/ti)
        xi = (xvect - x0) / ti
        solution = []
        for i in range(iMAX):
            solution.append(flx.nlers(qL,qR,xi[i]))
        yield solution

#===========================================================================


# Define Flux functions

#DIR = {
#'r' : np.array((0, 1)),
#'l' : np.array((-1, 0))
#}

from fluxes import NUM_MTD

def fv1D(solution, dx, qL, qR, tSTART, tEND, courant,
         nummtd, maxiter=10000):
    print nummtd
    time = tSTART
    for _ in range(maxiter):
        dt = checkCFL(solution, dx, courant, time, tEND,)
        if dt == -9999: break
        next_solution = np.empty(solution.shape)
        # UPWIND method #FIXME: perché secondo Lucas questo metodo è upwind non lo posso sapere a priori, o si?
        for i in range(2, len(solution) - 2):
#            fR = godSca(q(i),q(i+1));
#            fL = godSca(q(i-1),q(i));
            fluxL = NUM_MTD[nummtd](solution[i-1], solution[i], dx, dt)
            fluxR = NUM_MTD[nummtd](solution[i], solution[i+1], dx, dt)
            next_solution[i] = solution[i] - dt/dx * (fluxR - fluxL)
            #next_solution[i] = solution[i] + dt/dx * (fluxR - fluxL)
        # set the values for the boundary
        next_solution[0] = qL
        next_solution[-1] = qR
        solution = next_solution
        time += dt
        yield time, solution



def fv1DseconOrder(sol, dx, qL, qR, tSTART, tEND, courant,
                   nummtd, maxiter=10000):
    """

    """
    print nummtd
    time = tSTART
    for _ in range(maxiter):
        dt = checkCFL(sol, dx, courant, time, tEND,)
        if dt == -9999: break
        sol_p, sol_m = get_sol_interface(sol, dx, dt)
        # UPWIND method #FIXME: perché secondo Lucas questo metodo è upwind non lo posso sapere a priori, o si?
        for i in range(2, len(sol) - 2):
#            fR = godSca(q(i),q(i+1));
#            fL = godSca(q(i-1),q(i));
            fluxL = NUM_MTD[nummtd](sol_m[i-1], sol_p[i], dx, dt)
            fluxR = NUM_MTD[nummtd](sol_m[i], sol_p[i+1], dx, dt)
            sol[i] -= dt/dx * (fluxR - fluxL)
            # CHANGE + in -
        # set the values for the boundary
        sol[ 0] = qL
        sol[-1] = qR
        time += dt
        yield time, sol

#
# Define plot functions
#

def plot_burgers1Dexact(xL, xR, iMAX, x0, qL, qR, tEND, tMAX,
                        pause=0.005, ylimext = 0.15, **kargs):
    # compute the domine in space and time
    xvect, dx = np.linspace(xL, xR, iMAX, retstep = True)
    tvect, dt = np.linspace(0, tEND, tMAX, retstep = True)
    print('the space step dx is: {0}'.format(dx))
    print('the time step dt is: {0}'.format(dt))
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    qL_limit, qR_limit = (qL * (1. - ylimext ), qR * ( 1. + ylimext ) )
    # compute initial state t = 0
    line1, = ax.plot(xvect, flx.shock(xvect, x0, qL, qR), 'r-')
    plt.ylim( ( qL_limit, qR_limit ) )
    plt.pause(pause)
    # start the cycle for each step
    for solution in burgers1Dexact(xvect, tvect[1:], x0, qL, qR, iMAX):
        ax.clear()
        line1, = ax.plot(xvect, solution, 'r-')
        plt.ylim( ( qL_limit, qR_limit ) )
        plt.pause(pause)

def get_q_limit(qL, qR, ylimext):
    lst = [qL, qR]
    lst.sort()
    qL, qR = lst
    limiter = qR * ylimext
    return qL - limiter, qR + limiter




def plot_fv1D(xL, xR, x0, iMAX, qL, qR, tSTART, tEND, courant,
              nummtd, maxiter=10000, pause=0.005,
              ylimext = 0.15):
    # equidistant distribution of iMAX points between xL and xR
    xvect, dx = np.linspace(xL, xR, iMAX, retstep = True)
    # equidistant distribution of points between each spatial step
    xvect_b, dx_b = np.linspace(xL + dx/2., xR - dx/2., iMAX, retstep = True)
    # Initial condition
    solution = flx.shock(xvect_b, x0, qL, qR)
    # plot the first solution for t = o
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    qL_limit, qR_limit = get_q_limit(qL, qR, ylimext)
    #import pdb; pdb.set_trace()
    #print xvect_b
    print prettyp(solution[:10])
    line1, = ax.plot(xvect_b, solution, 'o')
    #plt.ylim( ( qL_limit, qR_limit ) )
    for time, solution in fv1DseconOrder(solution, dx, qL, qR, tSTART, tEND, courant, nummtd, maxiter=10000):
        plt.pause(pause)
        ax.clear()
        #print xvect_b
        print prettyp(solution[:10])
        plt.title('Current time t = {:6.5f}'.format(time))
        line1, = ax.plot(xvect_b, solution, 'o')
        xi = ( xvect - x0 ) / time
        line2, = ax.plot(xvect_b, [flx.nlers(qL, qR, xq) for xq in xi], 'r-')
        #line2, = ax.plot(xvect_b, flx.nlers2(qL, qR, x1), 'r-')
        #plt.ylim( ( qL_limit, qR_limit ) )




def main():
    from confPDE import *
    #plot_burgers1Dexact(xL, xR, iMAX, x0, qL, qR, tEND, tMAX, pause=0.005)

    # nummtd = godSca, lfSca, lwSca, forceSca, roeSca, oshermodSca
    plot_fv1D(xL, xR, x0, iMAX, qL, qR, tSTART, tEND, courant,
              'lwSca', maxiter=10000, pause=0.005)

main()

#>>> config.sections()
#[Burgers 1D exact]
#import pdb; pdb.set_trace()

