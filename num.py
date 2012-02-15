# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 15:06:53 2012

@author: pietro



How to use Solver1D
======================


>>> problem = Solver1D(advection_speed = -1)
>>> problem.compute('Lax Wendrof')
>>> problem.show()



"""
import numpy as np
import scipy as sc
import __future__
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import datetime

SOLVER = """Solver1D
Domain
    space : [{xL}, {xR}]
    time  : [{tSTART}, {tEND}]
Problem
    left value : {qL}
    right value : {qR}
    position discontinuity : {x0}
    advection speed : {advec}
Numerical Methods
    available : {mthds}
    used : {used}
    number of finite difference points : {iMAX}
    number of max iteration in time : {nMAX}
    courant number : {c}
    spatial resolution = {dx}
    time resolution = {dt}

"""



class Solver(object):
    def __init__(self, x0 = 0, xL = -1, xR = 1,
                 tSTART = 0, tEND = 0.2, iMAX = 200, nMAX = 10000):
        """
        """
        self.x0 = x0
        # Define the parameters of the mesh and time
        self.xL = xL           # position of the left boundary
        self.xR =  xR          # position of the right boundary
        self.tSTART =  tSTART  # initial time
        self.tEND =  tEND      # final time
        self.iMAX = iMAX       # number of finite difference points
        self.nMAX = nMAX       # number of max iteration
        self.xvect, self.dx = np.linspace(self.xL, self.xR, self.iMAX,
                                          retstep = True)
        # set a dictionary with all the numerical methods available
        self.mthds = {}
        self.used = None
        DTYPESOL = np.dtype([('time',np.float),
                             ('sol', np.float, (self.iMAX,))])
        self.numeric = np.empty((self.nMAX,), dtype = DTYPESOL)
        self.exact   = np.empty((self.nMAX,), dtype = DTYPESOL)
        self.ylimext = 0.15
        self.pause = 0.01


    def compute(self, funct):
        # set intial condition
        print("Start compute using: {mth}".format(mth=funct))
        self.used = funct
        start = datetime.datetime.now()
        self.dat[0] = self.tSTART, self.shock(self.x0)
        for itime in range(0, self.nMAX):
            time = self.dat['time'][itime]
            if  time + self.dt > self.tEND:
                dt = self.tEND - time
            else:
                dt = self.dt
            self.dat[itime+1] = time + dt, self.mthds[funct](itime, dt)
            if time + dt >= self.tEND: break
            #TODO: add the flush bar time+dt/self.tEND*100
        end = datetime.datetime.now()
        print("Time necessary to compute: {0}".format(end - start))

    def show(self):
        # FIXME:tart a process or a thread ?
        #plt.ion()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        first = True
        for ntime, nsol, etime, esol in zip(self.numeric, self.exact):
            if ntime!=0 or first:
                #import pdb; pdb.set_trace()
                first = False
                #print(time)
                plt.pause(self.pause)
                ax.clear()
                # plot the numerical solution
                line_numeric, = ax.plot(self.xvect, nsol,  'o')
                plt.title('{m} - Current time t = {t:6.5f}'.format(m=self.used,
                                                                   t=ntime))
                line_exact,   = ax.plot(self.xvect, esol, 'r-')
            else:
                break
        plt.show()

    def shock(self, x0):
        """
        >>> shock(range(10), 4, 2, 8)
        [[ 2.  2.  2.  2.  2.  8.  8.  8.  8.  8.]]
        """
        solution = np.empty(self.xvect.shape)
        solution[np.nonzero(self.xvect <= x0)] = self.qL
        solution[np.nonzero(self.xvect >  x0)] = self.qR
        return solution


class Solver1D(object):
    def __init__(self, advection_speed = 0, x0 = 0, xL = -1, xR = 1,
        tSTART = 0, tEND = 0.2, iMAX = 200, nMAX = 10000, courant_numb=0.5,
        qL = 1, qR = 0):
        """
        """
        self.advec = advection_speed # advection speed
        self.x0 = x0
        # Define the parameters of the mesh and time
        self.xL = xL           # position of the left boundary
        self.xR =  xR          # position of the right boundary
        self.tSTART =  tSTART  # initial time
        self.tEND =  tEND      # final time
        self.iMAX = iMAX       # number of finite difference points
        self.nMAX = nMAX       # number of max iteration
        self.xvect, self.dx = np.linspace(self.xL, self.xR, self.iMAX,
                                          retstep = True)
        self.courant = courant_numb  # Courant number
        # Boundary conditions
        self.qL = qL
        self.qR = qR
        # set a dictionary with all the numerical methods available
        self.mthds = {'Lax Wendrof'   : self.lw,
                      'Lax Friedrich' : self.lf,
                      'Warming Beam'  : self.wb,
                      'Fromm'         : self.fromm,}
        self.used = None
        self.solution = []
        DTYPESOL = np.dtype([('time', np.float64),
                             ('sol',  np.float64, (self.iMAX,))])
        self.dat = np.empty((self.nMAX,), dtype = DTYPESOL)
        self.ylimext = 0.15
        self.pause = 0.01
        dt = self.courant * self.dx / np.max( np.abs(self.advec) )
        self.dt = dt
        xexact = np.linspace(xL, xR, iMAX*5)
        self.xexact = xexact
        self.exact = np.zeros((len(self.xexact), ))
        self.exact_sol = []
        self.DTYPEEXACT = np.dtype([
                    ('time', np.float64),
                    ('sol',  np.float64, (len(self.xexact),))
                    ])
        self.last = None
        self.boundary = (np.array([ 0, -1]),)

    def __str__(self):
        return SOLVER.format(xL = self.xL, xR = self.xR,
                             tSTART = self.tSTART, tEND = self.tEND,
                             qL = self.qL, qR = self.qR, x0 = self.x0,
                             advec = self.advec,
                             mthds = ','.join(self.mthds.keys(),),
                             used = self.used,
                             iMAX = self.iMAX, nMAX = self.nMAX,
                             c = self.courant,
                             dx = self.dx, dt = self.dt)

#    def CFLcondition(time, self): #FIXME: tenerla o rimuoverla
#        """Return stop, and deltatime
#        stop i boolean and dt is float"""
#        #dt = self.maxflux * dx
#        if time + self.dt > self.tEND:
#            return False, self.tEND - time
#        if time >= self.tEND:
#            return True, None
#        return False, self.dt

    def compute(self, funct):
        # set intial condition
        print("Start compute using: {mth}".format(mth=funct))
        self.used = funct
        start = datetime.datetime.now()
        self.dat[0] = self.tSTART, self.shock(self.x0)
        for itime in range(0, self.nMAX):
            time = self.dat['time'][itime]
            if  time + self.dt > self.tEND:
                dt = self.tEND - time
            else:
                dt = self.dt
            self.dat[itime+1] = time + dt, self.mthds[funct](itime, dt)
            if time + dt >= self.tEND: break
        end = datetime.datetime.now()
        print("Time necessary to compute: {0}".format(end - start))
            #TODO: add the flush bar time+dt/self.tEND*100


    def show(self):
        # FIXME:tart a process or a thread ?
        #plt.ion()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        first = True
        for time, solution in self.dat:
            if time!=0 or first:
                #import pdb; pdb.set_trace()
                first = False
                #print(time)
                plt.pause(self.pause)
                ax.clear()
                # plot the numerical solution
                line_solution, = ax.plot(self.xvect, solution, 'o')
                plt.title('{m} - Current time t = {t:6.5f}'.format(m=self.used,
                                                                   t=time))
                # plot the exact solution #TODO: add a process to compute the exact
                x0 = self.x0 + self.advec*time
                line_exact, = ax.plot(self.xvect, self.shock(x0),  'r-')
            else:
                break
        plt.show()

    def shock(self, x0):
        """
        >>> shock(range(10), 4, 2, 8)
        [[ 2.  2.  2.  2.  2.  8.  8.  8.  8.  8.]]
        """
        solution = np.empty(self.xvect.shape)
        solution[np.nonzero(self.xvect <= x0)] = self.qL
        solution[np.nonzero(self.xvect >  x0)] = self.qR
        return solution

#======================

    def lw(self, itime, dt):
        """Lax Wendrof

        for i=2:IMAX-1
            q1(i) = q(i) - a*dt/dx/2       * (q(i+1)-       q(i-1)) ...
                         + a^2*dt^2/2/dx^2 * (q(i+1)-2*q(i)+q(i-1));
        end
        q1(1) = qL;
        q1(IMAX)=qR;
        """
        q0 = self.dat[itime]['sol']
        q1 = self.dat[itime+1]['sol']
        c =  self.advec * dt / self.dx
        for i in range(1, self.iMAX-1):
            # TODO: use numerical expression?
            q1[i] = q0[i] - c * 0.5     * ( q0[i+1] -           q0[i-1]) \
                          + c * c * 0.5 * ( q0[i+1] - 2*q0[i] + q0[i-1]) \
            # set boundary condition
        q1[(np.array([0, -1]),)] = self.qL, self.qR
        self.dat[itime+1]['time'] = self.dat[itime]['time'] + dt
        return q1 #FIXME: decidere se far restituire qualcosa dal metodo o no



    def lf(self, itime, dt):
        """Lax Friedrich
        for i=2:IMAX-1
            ap = 0.5*(a+abs(a));
            am = 0.5*(a-abs(a));
            q1(i) = 0.5 *           ( q(i+1) + q(i-1) ) ...
                   - ap * dt/dx/2 * ( q(i+1) - q(i-1) )
                   - am * dt/dx/2 * ( q(i+1) - q(i-1) );
        end
        q1(1) = qL;
        q1(IMAX)=qR;
        """
        q0 = self.dat[itime]['sol']
        q1 = self.dat[itime+1]['sol']
        ap = 0.5 * ( self.advec + abs( self.advec ) )
        am = 0.5 * ( self.advec - abs( self.advec ) )
        cp = ap * dt/self.dx/2.
        cm = am * dt/self.dx/2
        for i in range(1, len(q0)-1):
            q1[i] =  0.5 * ( q0[i+1] + q0[i-1]) \
                    - cp * ( q0[i+1] - q0[i-1]) \
                    - cm * ( q0[i+1] - q0[i-1]) \
            # set boundary condition
        q1[(np.array([0, -1]),)] = self.qL, self.qR
        return q1

    def wb(self, itime, dt):
        """Warming Beam

        for i=3:IMAX-2
            ap = 0.5*(a+abs(a));
            am = 0.5*(a-abs(a));
            q1(i) = q(i) - ap*dt/dx/2       * (q(i-2)-4*q(i-1)+3*q(i))...
                         + ap^2*dt^2/dx^2/2 * (q(i-2)-2*q(i-1)+  q(i))...
                         + am*dt/dx/2       * (q(i+2)-4*q(i+1)+3*q(i))...
                         + am^2*dt^2/dx^2/2 * (q(i+2)-2*q(i+1)+  q(i));
        end
        q1(1) = qL;
        q1(2) = qL;
        q1(IMAX)=qR;
        q1(IMAX-1)=qR;
        q=q1;
        """
        q0 = self.dat[itime]['sol']
        q1 = self.dat[itime+1]['sol']
        ap = 0.5 * ( self.advec + abs( self.advec ) )
        am = 0.5 * ( self.advec - abs( self.advec ) )
        cp = ap * dt/self.dx/2.
        cm = am * dt/self.dx/2
        cp2 = (cp*2)**2/2
        cm2 = (cm*2)**2/2
        for i in range(2, len(q0)-2):
            q1[i] = q0[i] - cp  * ( q0[i-2] - 4 * q0[i-1] + 3 * q0[i] ) \
                          + cp2 * ( q0[i-2] - 2 * q0[i-1] +     q0[i] ) \
                          + cm  * ( q0[i+2] - 4 * q0[i+1] + 3 * q0[i] ) \
                          + cm2 * ( q0[i+2] - 2 * q0[i+1] +     q0[i] )
        # set boundary condition
        q1[(np.array([0, 1, -1, -2]),)] = q0[(np.array([0, 1, -1, -2]),)]
        return q1

    def fromm(self, itime, dt):
        """
        for i=3:IMAX-2
            ap = 0.5*(a+abs(a));
            am = 0.5*(a-abs(a));
            q1(i) = q(i) - ap*dt/dx/4      *(q(i-2)-5*q(i-1)+3*q(i)+q(i+1))...
                         + ap^2*dt^2/dx^2/4*(q(i-2)-  q(i-1)-  q(i)+q(i+1))...
                         + am*dt/dx/4      *(q(i+2)-5*q(i+1)+3*q(i)+q(i-1))...
                         + am^2*dt^2/dx^2/4*(q(i+2)-  q(i+1)-  q(i)+q(i-1));
        end
        q1(1) = qL;
        q1(2) = qL;
        q1(IMAX)=qR;
        q1(IMAX-1)=qR;
        """
        q0 = self.dat[itime]['sol']
        q1 = self.dat[itime+1]['sol']
        ap = 0.5 * ( self.advec + abs( self.advec ) )
        am = 0.5 * ( self.advec - abs( self.advec ) )
        cp = ap * dt/self.dx/4.
        cm = am * dt/self.dx/4.
        cp2 = (cp*4)**2./4.
        cm2 = (cm*4)**2./4.
        for i in range(2, len(q0)-2):
            q1[i] = q0[i]-cp  * (q0[i-2] - 5 * q0[i-1] + 3 * q0[i] + q0[i+1]) \
                         +cp2 * (q0[i-2] -     q0[i-1] -     q0[i] + q0[i+1]) \
                         +cm  * (q0[i+2] - 5 * q0[i+1] + 3 * q0[i] + q0[i-1]) \
                         +cm2 * (q0[i+2] -     q0[i+1] -     q0[i] + q0[i-1])
        # set boundary condition
        q1[(np.array([0, 1, -1, -2]),)] = q0[(np.array([0, 1, -1, -2]),)]
        return q1


############################################################################

class Burgers(Solver1D):
    def __init__(self, advection_speed = 1, x0 = -3, xL = -5, xR = 5,
        tSTART = 0, tEND = 0.2, iMAX = 200, nMAX = 10000, courant_numb=0.5,
        qL = 1, qR = 3, pause = 0.01):
        Solver1D.__init__(self, advection_speed = advection_speed,
                 x0 = x0, xL = xL, xR = xR,
                 tSTART = tSTART, tEND = tEND,
                 iMAX = iMAX, nMAX = nMAX, courant_numb=courant_numb,
                 qL = qL, qR = qR)
        self.mthds = {'Lax Wendrof'   : self.lw,
                      'Lax Friedrich' : self.lf,
                      'force'         : self.force,
                      'Godunov'       : self.godunov,
                      'Roe'           : self.roe,
                      'Osher'         : self.osher}
        self.pause = pause

    def get_dt(self, solution):
        maxflux = np.max( np.abs( self.eigenvalue(solution) ) )
        return self.courant * self.dx / maxflux

    def flux(self, q):
        """
        function y = flux(q)
            y = 0.5*q.^2;   % componentwise application
            %y = q^4/4;
        """
        return 0.5*q*q #q*q*q*q/4. # Example 21/23, pag. 31/32


    def eigenvalue(self, q):
        """
        function y = eigenvalue(q)
            y = q;
            %y = q.^3;
        """
        return q #q*q*q # Example 21/23, pag. 31/32

    def shock(self, x0, qL, qR, xvect=None):
        """
        >>> shock(range(10), 4, 2, 8)
        [[ 2.  2.  2.  2.  2.  8.  8.  8.  8.  8.]]
        """
        if xvect==None: xvect = self.xvect
        solution = np.empty(xvect.shape)
        solution[np.nonzero(xvect <= x0)] = qL
        solution[np.nonzero(xvect > x0)] = qR
        return solution

    def nlers(self, qL, qR, xi):
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
        aL = self.eigenvalue(qL);
        aR = self.eigenvalue(qR);
        if aL > aR:
            # is a Shock wave
            s = ( self.flux(qR) - self.flux(qL) ) / ( qR - qL );
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
                q = self.newton(q, xi)
        return q

    def g(self, q, xi):
        """
        function y = g(q,xi)
            y = eigenvalue(q)-xi;
        """
        return self.eigenvalue(q)-xi


    def dg(self, q, xi, eta=1e-7):
        """
        function y = dg(q,xi)
            eta=1e-7;
            y = (g(q+eta,xi)-g(q-eta,xi))/(2*eta);
        """
        return ( self.g(q + eta, xi) - self.g(q - eta, xi) ) / ( 2 * eta )


    def newton(self, q0, xi, tol = 1e-12, iMAX = 100):
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
            if(abs( self.g(q, xi) ) < tol):
                break
            dq = self.g(q, xi) / self.dg(q, xi)
            q = q-dq
        #if i == iMAX and abs( g(q, xi) ) < tol: print('Not converged using Newton')
        #print('Newton number of iteration to converge: {0:d}'.format(i))
        return q

    def lf(self, qL, qR, dx, dt):
        return 0.5 * (self.flux(qL) + self.flux(qR)) - 0.5 * dx/dt * (qR - qL)

    def lw(self, qL, qR, dx, dt):
        """Lax-Wendrof flux, pag 54, formula 4.9

        function flu = lwSca(qL,qR,dx,dt)
            a = 0.5 * ( eigenvalue(qL) + eigenvalue(qR) );
            flu = 1/2 * ( flux(qR)+flux(qL) ) - 1/2 * dt/dx * a^2 * ( qR - qL );
        """
        #import pdb; pdb.set_trace()
        a = 0.5 * ( self.eigenvalue(qL) + self.eigenvalue(qR) )
        return   0.5 * ( self.flux(qR) + self.flux(qL) ) \
               - 0.5 * dt/dx * a*a * ( qR - qL )

    def force(self, qL, qR, dx, dt):
        """FORCE flux, pag 54, formula 4.13

        function flu = forceSca(qL,qR,dx,dt)
            fLW = lwSca(qL,qR,dx,dt);
            fLF = lfSca(qL,qR,dx,dt);
            flu = 0.5*(fLW+fLF);
        """
        fLW = self.lw(qL, qR, dx, dt)
        fLF = self.lf(qL, qR, dx, dt)
        return 0.5 * ( fLW + fLF )

    def godunov(self, qL, qR, dx, dt):
        """Godunov flux, pag 55, formula 4.14

        function flu = godSca(qL,qR)
            flu = flux(NLERS(qL,qR,0));
        """
        return self.flux( self.nlers(qL, qR, 0) )

    def roe(self, qL, qR, dx, dt):
        """Roe flux, pag 63, formula 4.71

        function flu = forceSca(qL,qR,dx,dt)
            % Only for burgers
            a = 0.5*(qR+qL);
            flu = 0.5*(flux(qL)+flux(qR))-0.5*abs(a)*(qR-qL);
        """
        a = 0.5 * ( qR + qL )
        return   0.5 * ( self.flux(qL) + self.flux(qR) ) \
               - 0.5 * abs(a) * ( qR - qL )

    def osher(self, qL, qR, dx, dt): # :
        """Osher flux, pag 65, formula 4.83

        function flu = oshermodSca(qL,qR,dx,dt)
            aL = eigenvalue(qL);
            aR = eigenvalue(qR);
            flu = 0.5*(flux(qL)+flux(qR))-0.5*((abs(aL)+abs(aR))/2)*(qR-qL);
        """
        aL, aR = self.eigenvalue(qL), self.eigenvalue(qR);
        return   0.5 * ( self.flux(qL) + self.flux(qR) ) \
               - 0.5 * ( ( abs(aL) + abs(aR) ) * 0.5) * ( qR - qL )


    def get_solution(self, itime, dt, funct):
        q0 = self.dat[itime]['sol']
        q1 = self.dat[itime+1]['sol']
        #import pdb; pdb.set_trace()
        for i in range(1, len(q0)-1):
            #if q0[i-1] != q0[i+1]: import pdb; pdb.set_trace()
            fluxL = self.mthds[funct](q0[i-1], q0[i], self.dx, dt)
            fluxR = self.mthds[funct](q0[i], q0[i+1], self.dx, dt)
            q1[i] =  q0[i] - dt/self.dx * (fluxR - fluxL)
            # set boundary condition
        q1[self.boundary] = q0[self.boundary]
        return q1

    def compute_numerical(self, funct):
        # set intial condition
        print("Start compute using: {mth}".format(mth=funct))
        self.used = funct
        start = datetime.datetime.now()
        self.dat[0] = self.tSTART, self.shock(self.x0, self.qL, self.qR)
        for itime in range(0, self.nMAX):
            time = self.dat['time'][itime]
            dt = self.get_dt(self.dat['sol'][itime])
            if  time + dt > self.tEND: dt = self.tEND - time
            self.dat[itime+1] = time + dt, self.get_solution(itime, dt, funct)
            #print self.dat[itime+1]
            if time + dt >= self.tEND: break
        end = datetime.datetime.now()
        print("Time necessary to compute: {0}".format(end - start))
            #TODO: add the flush bar time+dt/self.tEND*100

#    def burgers1Dexact(xvect, tvect, x0, qL, qR, iMAX):
#        for ti in tvect:
#    #        yield flx.nlers2(qL, qR, (xvect-x0)/ti)
#            xi = (xvect - x0) / ti
#            solution = []
#            for i in range(iMAX):
#                solution.append(self.nlers(qL,qR,xi[i]))
#            yield solution

    def compute_exact(self):
        self.last = np.nonzero(self.dat['time'])[0][-1]+1
        #import pdb; pdb.set_trace()
        self.exact_sol = np.zeros((self.last,), dtype=self.DTYPEEXACT )
        print("Start compute exact solution")
        start = datetime.datetime.now()
        for itime, time in enumerate(self.dat['time'][:self.last]):
            if itime == 0:
                self.exact_sol[itime] = time, self.shock(self.x0,
                                                         self.qL, self.qR,
                                                         self.xexact)
            else:
                for i, xt in enumerate((self.xexact - self.x0)/time):
                    #import pdb; pdb.set_trace()
                    self.exact[i] = self.nlers(self.qL, self.qR, xt)
                #print(itime)
                self.exact_sol[itime] = time, self.exact
        end = datetime.datetime.now()
        print("Time necessary to compute: {0}".format(end - start))

    def show(self):
        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1)
        first = True
        for (time_d, sol_d), (time_e, sol_e) in zip(self.dat, self.exact_sol):
            if time_d!=0 or first:
                first = False
                plt.pause(self.pause)
                ax1.clear()
                plt.title('{m} - Current time t = {t:6.5f}'.format(m=self.used,
                                                                   t=time_d))
                line1e, = ax1.plot(self.xexact, sol_e, 'r-')
                line1d, = ax1.plot(self.xvect,  sol_d, 'o')
                #TODO: add a process to compute the exact
            else:
                break
        plt.show()

############################################################################




class Solver1Dsys(object):
    def __init__(self, advection_speed = np.array([[0., 2.], [ 1., 0.]]),
                 R = np.array([[-np.sqrt(2), np.sqrt(2)],[1., 1.]]),
                 D = np.array([[-np.sqrt(2), 0],[0., np.sqrt(2)]]),
                 x0 = 0, xL = -1, xR = 1,
                 tSTART = 0, tEND = 0.2,
                 iMAX = 200, nMAX = 10000, courant_numb=0.5,
                 qL = np.array([2., 1.]), qR = np.array([1., 0.])):
        """
        """

        self.x0 = x0
        # Define the parameters of the mesh and time
        self.xL = xL           # position of the left boundary
        self.xR =  xR          # position of the right boundary
        self.tSTART =  tSTART  # initial time
        self.tEND =  tEND      # final time
        self.iMAX = iMAX       # number of finite difference points
        self.nMAX = nMAX       # number of max iteration
        self.xvect, self.dx = np.linspace(self.xL, self.xR, self.iMAX,
                                          retstep = True)
        self.courant = courant_numb  # Courant number

        # Boundary conditions
        self.qL = np.array(qL)
        self.qR = np.array(qR)
        self.A = np.array(advection_speed) # advection speed
        self.R0 = np.array(R)
        self.D0 = np.array(D)
        B = np.dot(self.R0, np.dot(np.abs( self.D0 ), np.linalg.inv( self.R0 )))
        self.B = B
        lambd, R, D, Id, cL, cR = self.get_exactparameter()
        self.lambd, self.R, self.D, self.Id, self.cL, self.cR = lambd, R, D, Id, cL, cR
        Aplus =  0.5 * ( self.A + self.B )
        self.Aplus = Aplus
        Aminus = 0.5 * ( self.A - self.B )
        self.Aminus = Aminus

        # set a dictionary with all the numerical methods available
        self.mthds = {'Lax Wendrof'   : self.lw,
                      'Lax Friedrich' : self.lf,
                      'Warming Beam'  : self.wb,
                      'Fromm'         : self.fromm,
                      'WENO'          : self.weno,}
                      #'DGRK'          : self.dgrk}
        self.used = None
        self.solution = []
        self.DTYPESOL = np.dtype([
                             ('time', np.float64),
                             ('sol',  np.float64, (self.iMAX, len(self.qL)))])
        self.dat = np.zeros((self.nMAX,), dtype = self.DTYPESOL)

        self.ylimext = 0.15
        self.pause = 0.01
        # dt = c * dx / max( max( abs( lambda ) ) );
        dt = self.courant * self.dx / np.max( np.abs(self.lambd) )
        self.dt = dt
        self.boundary = (np.array([ 0,  1, -1, -2]),)
        xexact = np.linspace(xL, xR, iMAX*5)
        self.xexact = xexact
        self.exact = np.zeros((len(self.xexact), 2))
        self.exact_sol = []
        self.DTYPEEXACT = np.dtype([
                    ('time', np.float64),
                    ('sol',  np.float64, (len(self.xexact), len(self.qL)))
                    ])
        self.last = None

        #=============== used only for DGRK



    def __str__(self):
        return SOLVER.format(xL = self.xL, xR = self.xR,
                             tSTART = self.tSTART, tEND = self.tEND,
                             qL = self.qL, qR = self.qR, x0 = self.x0,
                             advec = self.A,
                             mthds = ','.join(self.mthds.keys(),),
                             used = self.used,
                             iMAX = self.iMAX, nMAX = self.nMAX,
                             c = self.courant,
                             dx = self.dx, dt = self.dt)


    def get_exactparameter(self, qL=None, qR=None):
        if (qL, qR) == (None, None): qL, qR = self.qL, self.qR
        d, R = np.linalg.eig(self.A)
        D = np.diag(d)
        iR = np.linalg.inv(R)
        cL = np.dot(iR, qL)
        cR = np.dot(iR, qR)
        Id = sc.eye(len(qL))
        return d, R, D, Id, cL, cR

    def compute_numerical(self, funct):
        # set intial condition
        print("Start compute using: {mth}".format(mth=funct))
        self.used = funct
        start = datetime.datetime.now()
        self.dat[0] = self.tSTART, self.shock(self.x0)
        for itime in range(0, self.nMAX):
            time = self.dat['time'][itime]
            if  time + self.dt > self.tEND:
                dt = self.tEND - time
            else:
                dt = self.dt
            self.dat[itime+1] = time + dt, self.mthds[funct](itime, dt)
            if time + dt >= self.tEND: break
        end = datetime.datetime.now()
        print("Time necessary to compute: {0}".format(end - start))
            #TODO: add the flush bar time+dt/self.tEND*100

    def ers(self, qL, qR, A, xt):
        d, R = np.linalg.eig(A)
        D = np.diag(d)
        iR = np.linalg.inv(R)
        cL = np.dot(iR, qL) #qL
        cR = np.dot(iR, qR) #qR
        Id = np.identity(2)
        return  np.dot(0.5 * R, np.dot(Id + np.sign( D - xt * Id ),cL)) \
                  + np.dot(0.5 * R, np.dot(Id - np.sign( D - xt * Id ),cR))


    def compute_exact(self, time):
        if time == 0:
            self.exact = self.shock(self.x0, self.xexact)
            return self.exact
        else:
            for i, xt in enumerate(self.xexact/time):
                 self.exact[i,:] = self.ers(self.qL, self.qR, self.A, xt)
            return self.exact

    def compute_exact2(self):
        self.last = np.nonzero(self.dat['time'])[0][-1]+1
        #import pdb; pdb.set_trace()
        self.exact_sol = np.zeros((self.last,), dtype=self.DTYPEEXACT )
        print("Start compute exact solution")
        start = datetime.datetime.now()
        for itime, time in enumerate(self.dat['time'][:self.last]):
            if itime == 0:
                self.exact_sol[itime] = time, self.shock(self.x0, self.xexact)
            else:
                for i, xt in enumerate(self.xexact/time):
                    #import pdb; pdb.set_trace()
                    self.exact[i,:] = self.ers(self.qL, self.qR, self.A, xt)
                #print(itime)
                self.exact_sol[itime] = time, self.exact
        end = datetime.datetime.now()
        print("Time necessary to compute: {0}".format(end - start))


    def show(self):
        # FIXME:tart a process or a thread ?
        #plt.ion()
        fig = plt.figure()
        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,1,2)
        first = True
        for time, solution in self.dat:
            if time!=0 or first:
                #import pdb; pdb.set_trace()
                first = False
                qexact = self.compute_exact(time)
                plt.pause(self.pause)
                ax1.clear()
                ax2.clear()
                qexact = self.compute_exact(time)
                #import pdb; pdb.set_trace()
                line1e, = ax1.plot(self.xexact, qexact[:,0], 'r-')
                line1, = ax1.plot(self.xvect, solution[:,0], 'o')
                plt.title('Current time t = {:6.5f}'.format(time))
                line2e, = ax2.plot(self.xexact, qexact[:,1], 'r-')
                line2, = ax2.plot(self.xvect, solution[:,1], 'o')
                # plot the numerical solution
                plt.title('{m} - Current time t = {t:6.5f}'.format(m=self.used,
                                                                   t=time))
                # plot the exact solution #TODO: add a process to compute the exact
            else:
                break
        plt.show()


    def show2(self):
        # FIXME:tart a process or a thread ?
        self.compute_exact2()
        fig = plt.figure()
        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,1,2)
        first = True
        for (time_d, sol_d), (time_e, sol_e) in zip(self.dat, self.exact_sol):
            if time_d!=0 or first:
                first = False
                plt.pause(self.pause)
                ax1.clear()
                ax2.clear()
                plt.title('{m} - Current time t = {t:6.5f}'.format(m=self.used,
                                                                   t=time_d))
                line1e, = ax1.plot(self.xexact, sol_e[:,0], 'r-')
                line1d, = ax1.plot(self.xvect,  sol_d[:,0], 'o')
                line2e, = ax2.plot(self.xexact, sol_e[:,1], 'r-')
                line2d, = ax2.plot(self.xvect,  sol_d[:,1], 'o')
                #TODO: add a process to compute the exact
            else:
                break
        plt.show()

    def show3(self):
        # FIXME:tart a process or a thread ?
        self.compute_exact2()
        fig = plt.figure()
        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,1,2)
        axcolor = 'lightgoldenrodyellow'
        axtime = plt.axes([0.25, 0.25, 0.65, 0.03], axisbg=axcolor)
        time_min, time_max = self.dat['time'][0], self.dat['time'][self.last-1]
        #ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f', closedmin=True, closedmax=True, slidermin=None, slidermax=None, dragging=True, **kwargs
        first = True
        for (time_d, sol_d), (time_e, sol_e) in zip(self.dat, self.exact_sol):
            if time_d!=0 or first:
                first = False
                plt.pause(self.pause)
                ax1.clear()
                ax2.clear()
                plt.title('{m} - Current time t = {t:6.5f}'.format(m=self.used,
                                                                   t=time_d))
                line1e, = ax1.plot(self.xexact, sol_e[:,0], 'r-')
                line1d, = ax1.plot(self.xvect,  sol_d[:,0], 'o')
                line2e, = ax2.plot(self.xexact, sol_e[:,1], 'r-')
                line2d, = ax2.plot(self.xvect,  sol_d[:,1], 'o')
                slider_time = Slider(axtime, 'time', time_min,
                                     time_max, valinit=time_d)
                #TODO: add a process to compute the exact
            else:
                break

        # Update slider
        t_data = self.dat['time'][:self.last]
        def update(val):
            time = slider_time.val
            difftime = t_data - time
            itime = difftime.argmin()
            plt.title('{m} - Current time t = {t:6.5f}'.format(m=self.used,
                      t=self.dat[itime]['time']))
            line1e.set_ydata(self.exact_sol[itime]['sol'][:,0])
            line1d.set_ydata(      self.dat[itime]['sol'][:,0])
            line2e.set_ydata(self.exact_sol[itime]['sol'][:,1])
            line2d.set_ydata(      self.dat[itime]['sol'][:,1])
            plt.draw()
        slider_time.on_changed(update)
        plt.show()


    def shock(self, x0, xvect = None):
        """
        %% Set and plot initial conditions
        for i=1:IMAX
            if(x(i)<=0.)
                q(:,i) = qL;
            else
                q(:,i) = qR;
            end
        end
        """
        if xvect==None: xvect = self.xvect
        solution = np.empty((len(xvect),len(self.qL)))
        solution[np.nonzero(xvect <= x0), :] = self.qL
        solution[np.nonzero(xvect >  x0), :] = self.qR
        return solution


    def lw(self, itime, dt):
        """
        q1(:,i) = q(:,i) - A * dt/dx/2     * (q(:,i+1) -        q(:,i-1) )...
                         + A*A*dt^2/dx^2/2 * (q(:,i+1)-2*q(:,i)+q(:,i-1) );
        q1(:,i) = q(:,i) - A *     dt/dx/2 *    ( q(:,i+1) - q(:,i-1) ) ...
                         + A * A * dt^2/2/dx^2 *( q(:,i+1) -2*q(:,i) + q(:,i-1));

        q1(:,i) = q(:,i)-A*dt/dx/2      *(q(:,i+1)         -q(:,i-1)) ...
                        +A*A*dt^2/2/dx^2*(q(:,i+1)-2*q(:,i)+q(:,i-1));

        """
        q0 = self.dat[itime]['sol']
        q1 = self.dat[itime+1]['sol']
        Cp  = self.A * dt / self.dx * 0.5
        Cp2 = np.dot(self.A, self.A) * dt**2 / self.dx**2 * 0.5
        for i in range(1, len(q0)-2):
            q1[i,:] = q0[i,:] - np.dot(Cp,  q0[i+1,:] -             q0[i-1,:])\
                              + np.dot(Cp2, q0[i+1,:] - 2*q0[i,:] + q0[i-1,:])
        # set boundary condition
        q1[self.boundary,:] = q0[self.boundary,:]
        return q1

    def lf(self, itime, dt):
        """
        q1(:,i) =             0.5 * ( q(:,i+1) + q(:,i-1) ) ...
                 - Aplus*dt/dx/2  * ( q(:,i+1) - q(:,i-1) ) ...
                 - Aminus*dt/dx/2 * ( q(:,i+1) - q(:,i-1) );

        """
        q0 = self.dat[itime]['sol']
        q1 = self.dat[itime+1]['sol']
        Cp = self.Aplus * dt / self.dx * 0.5
        Cm = self.Aminus * dt / self.dx * 0.5
        for i in range(1, len(q0)-2):
            q1[i,:] =   0.5 * (     q0[i+1,:] + q0[i-1,:] ) \
                      - np.dot( Cp, q0[i+1,:] - q0[i-1,:]) \
                      - np.dot( Cm, q0[i+1,:] - q0[i-1,:])
        # set boundary condition
        q1[self.boundary,:] = q0[self.boundary,:]
        return q1

    def wb(self, itime, dt):
        """Warming Beam
        for i=3:IMAX-2
            q1(:,i) = q(:,i) - Aplus*dt/dx/2             * (q(:,i-2)-4*q(:,i-1)+3*q(:,i))...
                             + Aplus*Aplus*dt^2/dx^2/2   * (q(:,i-2)-2*q(:,i-1)+  q(:,i))...
                             + Aminus*dt/dx/2            * (q(:,i+2)-4*q(:,i+1)+3*q(:,i))...
                             + Aminus*Aminus*dt^2/dx^2/2 * (q(:,i+2)-2*q(:,i+1)+  q(:,i))    ;
        end
        """
        q0 = self.dat[itime]['sol']
        q1 = self.dat[itime+1]['sol']
        Ap = self.Aplus * dt / self.dx
        Ap2 = np.dot(Ap,Ap)
        Am = self.Aminus * dt / self.dx
        Am2 = np.dot(Am,Am)
        for i in range(2, len(q0)-2):
            q1[i,:] = q0[i,:]-np.dot(Ap * 0.5, q0[i-2,:]-4*q0[i-1,:]+3*q0[i,:])\
                             +np.dot(Ap2* 0.5, q0[i-2,:]-2*q0[i-1,:]+  q0[i,:])\
                             +np.dot(Am * 0.5, q0[i+2,:]-4*q0[i+1,:]+3*q0[i,:])\
                             +np.dot(Am2* 0.5, q0[i+2,:]-2*q0[i+1,:]+  q0[i,:])
        # set boundary condition
        q1[self.boundary,:] = q0[self.boundary,:]
        return q1

    def fromm(self, itime, dt):
        """
        q1(:,i) = q(:,i)-Aplus*dt/dx/4       *(q(:,i-2)-5*q(:,i-1)+3*q(:,i)+q(:,i+1)) ...
                        +Aplus^2*dt^2/dx^2/4 *(q(:,i-2)-  q(:,i-1)-  q(:,i)+q(:,i+1)) ...
                        +Aminus*dt/dx/4      *(q(:,i+2)-5*q(:,i+1)+3*q(:,i)+q(:,i-1)) ...
                        +Aminus^2*dt^2/dx^2/4*(q(:,i+2)-  q(:,i+1)-  q(:,i)+q(:,i-1));
        """
        q0 = self.dat[itime]['sol']
        q1 = self.dat[itime+1]['sol']
        Ap = self.Aplus * dt / self.dx
        Ap2 = np.dot(Ap,Ap)
        Am = self.Aminus * dt / self.dx
        Am2 = np.dot(Am,Am)
        Cp, Cm, Cp2, Cm2 = Ap * 0.25, Am * 0.25, Ap2 * 0.25,  Am2 * 0.25
        for i in range(2, len(q0)-2):
            q1[i,:] = q0[i,:]-np.dot(Cp, q0[i-2,:]-5*q0[i-1,:]+3*q0[i,:]+q0[i+1,:]) \
                             +np.dot(Cp2,q0[i-2,:]-  q0[i-1,:]-  q0[i,:]+q0[i+1,:])  \
                             +np.dot(Cm, q0[i+2,:]-5*q0[i+1,:]+3*q0[i,:]+q0[i-1,:]) \
                             +np.dot(Cm2,q0[i+2,:]-  q0[i+1,:]-  q0[i,:]+q0[i-1,:])
        # set boundary condition
        q1[self.boundary,:] = q0[self.boundary,:]
        return q1

    def weno_qm_qp(self, itime, dt):
        """
        for i=3:length(q)-2
            eta=1e-14;
            r=8;
            % Compute polynomial coefficients for the left stencil
            l0 = q(:,i)';
            l1 = 1/4*q(:,i-2)'-q(:,i-1)'+3/4*q(:,i)';
            l2 = 1/12*q(:,i-2)'-1/6*q(:,i-1)'+1/12*q(:,i)';
            % Compute polynomial coefficients for the central stencil
            c0 = q(:,i)';
            c1 = -1/4*q(:,i-1)'+1/4*q(:,i+1)';
            c2 = 1/12*q(:,i-1)'-1/6*q(:,i)'+1/12*q(:,i+1)';
            % Compute polynomial coefficients for the right stencil
            r0 = q(:,i)';
            r1 = -3/4*q(:,i)'+q(:,i+1)'-1/4*q(:,i+2)';
            r2 = 1/12*q(:,i)'-1/6*q(:,i+1)'+1/12*q(:,i+2)';
            % Compute the Oscillation Indicator (OI)
            sigmaL = 156*l2.^2+4*l1.^2;
            sigmaR = 156*r2.^2+4*r1.^2;
            sigmaC = 156*c2.^2+4*c1.^2;
            % Compute the non-linearized non-linear weights
            whatL = 1./(sigmaL+eta).^r;
            whatR = 1./(sigmaR+eta).^r;
            whatC = 100000./(sigmaC+eta).^r;
            sumwhat = whatL+whatC+whatR;
            % Normalize non-linear weights
            wL = whatL./sumwhat;
            wR = whatR./sumwhat;
            wC = whatC./sumwhat;
            % Compute weights for final polynomial
            omega0 = wL.*l0+wC.*c0+wR.*r0;
            omega1 = wL.*l1+wC.*c1+wR.*r1;
            omega2 = wL.*l2+wC.*c2+wR.*r2;
            % Compute x_i+1/2 and x_i-1/2 values, spatial and temporal
            % derivatives
            q_right     = omega0+omega1+omega2;
            dq_x_right  = (2*omega1+6*omega2)/dx;
            dq_xx_right = omega2*12/dx^2;
            dq_t_right   = (-A*dq_x_right')';
            dq_tt_right  = (A^2*dq_xx_right')';
            q_left      = omega0-omega1+omega2;
            dq_x_left   = (omega1*2-omega2*6)/dx;
            dq_xx_left  = omega2*12/dx^2;
            dq_t_left    = (-A*dq_x_left')';
            dq_tt_left   = (A^2*dq_xx_left')';
            qm(:,i) = q_right+dt/2*dq_t_right+dt^2/6*dq_tt_right;  % right side of cell
            qp(:,i) = q_left+dt/2*dq_t_left+dt^2/6*dq_tt_left;  % left side of cell
        end
        qm(:,1)=q(:,1);qm(:,2)=q(:,1);qm(:,3)=q(:,1);
        qp(:,1)=q(:,1);qp(:,2)=q(:,1);qp(:,3)=q(:,1);
        qm(:,IMAX-2)=q(:,IMAX-1);qm(:,IMAX-1)=q(:,IMAX-1);qm(:,IMAX)=q(:,IMAX-1);
        qp(:,IMAX-2)=q(:,IMAX-1);qp(:,IMAX-1)=q(:,IMAX-1);qp(:,IMAX)=q(:,IMAX-1);
        % Compute fluxes and update the numerical solution
        for i=3:IMAX-2
            fp = A*ERS(qm(:,i),qp(:,i+1),A,0);
            fm = A*ERS(qm(:,i-1),qp(:,i),A,0);
            q1(:,i) = q(:,i)-dt/dx*(fp-fm);
        end
        q1(:,1) = qL;q1(:,2) = qL;q1(:,3) = qL;
        q1(:,IMAX)=qR;q1(:,IMAX-1)=qR;q1(:,IMAX-2)=qR;
        q=q1;
        """
        #import pdb; pdb.set_trace()
        r=8
        eta = 1e-14
        q0 = self.dat[itime]['sol']
        A2 = np.dot(self.A, self.A)
        qm, qp = np.empty(q0.shape), np.empty(q0.shape)
        for i in range(2, len(q0)-2):
            # Compute polynomial coefficients for the left stencil
            l0 = q0[i,:]
            l1 = 1./4.  * q0[i-2,:] -         q0[i-1,:] + 3./4.  * q0[i,:]
            l2 = 1./12. * q0[i-2,:] - 1./6. * q0[i-1,:] + 1./12. * q0[i,:]
            # Compute polynomial coefficients for the central stencil
            c0 =  q0[i,:]
            c1 = - 1./4. * q0[i-1,:]                    + 1./4.  * q0[i+1,:]
            c2 =  1./12. * q0[i-1,:] - 1./6. * q0[i,:]  + 1./12. * q0[i+1,:]
            # Compute polynomial coefficients for the right stencil
            r0 =           q0[i,:]
            r1 = - 3./4. * q0[i,:] +         q0[i+1,:] - 1./4.  * q0[i+2,:]
            r2 =  1./12. * q0[i,:] - 1./6. * q0[i+1,:] + 1./12. * q0[i+2,:]
            # Compute the Oscillation Indicator (OI)
            sigmaL = 156. * l2**2. + 4. * l1**2.
            sigmaR = 156. * r2**2. + 4. * r1**2.
            sigmaC = 156. * c2**2. + 4. * c1**2.
            # Compute the non-linearized non-linear weights
            whatL = 1./(sigmaL+eta)**r
            whatR = 1./(sigmaR+eta)**r
            whatC = 100000./(sigmaC+eta)**r
            sumwhat = whatL + whatC + whatR
            # Normalize non-linear weights
            wL = whatL/sumwhat
            wR = whatR/sumwhat
            wC = whatC/sumwhat
            # Compute weights for final polynomial
            omega0 = wL*l0 + wC*c0 + wR*r0
            omega1 = wL*l1 + wC*c1 + wR*r1
            omega2 = wL*l2 + wC*c2 + wR*r2
            # Compute x_i+1/2 and x_i-1/2 values, spatial and temporal
            # derivatives
            q_right     = omega0 + omega1 + omega2
            dq_x_right  = ( 2. * omega1 + 6. * omega2 )/self.dx
            dq_xx_right = omega2 * 12./self.dx**2.
            dq_t_right   = np.dot(-self.A, dq_x_right)
            dq_tt_right  = np.dot(A2, dq_xx_right)
            q_left      =  omega0 - omega1 + omega2
            dq_x_left   = ( omega1 * 2. - omega2 * 6. )/self.dx
            dq_xx_left  = omega2 * 12./self.dx**2.
            dq_t_left    = np.dot(-self.A, dq_x_left)
            dq_tt_left   = np.dot(A2, dq_xx_left)
            # right side of cell
            qm[i,:] = q_right + dt/2. * dq_t_right + dt**2./6. * dq_tt_right
            # left side of cell
            qp[i,:] = q_left + dt/2. * dq_t_left + dt**2./6. * dq_tt_left
            #if i==10: import pdb; pdb.set_trace()
        boundary = (np.array([ 0,  1, 2, -1,-2,-3]),)
        qm[boundary,:] = q0[boundary,:]
        qp[boundary,:] = q0[boundary,:]
        return qm, qp


    def weno(self, itime, dt):
        """WENO
        % WENO 3 order
        for i=3:IMAX-2
            fp = A*ERS(qm(:,i),qp(:,i+1),A,0);
            fm = A*ERS(qm(:,i-1),qp(:,i),A,0);
            q1(:,i) = q(:,i)-dt/dx*(fp-fm);
        end
        q1(:,1) = qL;q1(:,2) = qL;q1(:,3) = qL;
        q1(:,IMAX)=qR;q1(:,IMAX-1)=qR;q1(:,IMAX-2)=qR;
        q=q1;
        """
        #import pdb; pdb.set_trace()
        q0 = self.dat[itime]['sol']
        q1 = self.dat[itime+1]['sol']
        qm, qp = self.weno_qm_qp(itime, dt)
        for i in range(2, len(q0)-2):
            fp = np.dot(self.A, self.ers(qm[i,:],  qp[i+1,:],self.A,0))
            fm = np.dot(self.A, self.ers(qm[i-1,:],qp[i,:],self.A,0))
            q1[i,:] = q0[i,:] - dt/self.dx * ( fp - fm )
        # set boundary conditions
        q1[self.boundary,:] = q0[self.boundary,:]
        return q1







class DGRK(Solver1Dsys):
    def __init__(self, advection_speed = np.array([[0., 2.], [ 1., 0.]]),
                 R = np.array([[-np.sqrt(2), np.sqrt(2)],[1., 1.]]),
                 D = np.array([[-np.sqrt(2), 0],[0., np.sqrt(2)]]),
                 x0 = 0, xL = -1, xR = 1,
                 tSTART = 0, tEND = 0.2,
                 iMAX = 200, nMAX = 10000, courant_numb=0.5,
                 qL = np.array([2., 1.]), qR = np.array([1., 0.]),
                 accuracy=2,
                 M=[(1, 0, 0),( 0, 1./3., 0),(0, 0, 1./5.)],
                 K=[(0, 0, 0),(2, 0, 0),(0, 2, 0)] ):
        #super(Solver1D, self).__init__() #return AttributeError: 'DGRK' object has no attribute 'dx'
        Solver1Dsys.__init__(self, advection_speed = advection_speed,
                 R = R, D = D, x0 = x0, xL = xL, xR = xR,
                 tSTART = tSTART, tEND = tEND,
                 iMAX = iMAX, nMAX = nMAX, courant_numb=courant_numb,
                 qL = qL, qR = qR)
        # add the new methods
        self.mthds['DGRK'] = self.dgrk
        self.mthds['DGRK2'] = self.dgrk2
        self.nVar,self.temp = self.A.shape
        #return TypeError: object.__init__() takes no parameters
        self.accuracy = int(accuracy)         # polynomial degree of the test and basis functions
        self.nDOFs = self.accuracy+1;    # number of spatial degrees of freedom
        self.nDOF  = (self.accuracy+1)^2 # number of space-time degrees of freedom
        self.courant = 1./(2.*self.accuracy+1.);  # Courant number
        dt = self.courant * self.dx / np.max( np.abs(self.lambd) )
        self.dt = dt
        M = self.dx * np.array(M, dtype = np.float64)
        self.M = M
        iM = np.linalg.inv(self.M)
        self.iM = iM
        # Define mass matrix for int(diff(phi_k,x)*phi_l)
        self.K = np.array(K, dtype = np.float64)
        qhat0 = np.zeros((self.nVar, self.nDOFs, self.iMAX))
        self.qhat0 = qhat0
        self.qhat0[:,0,:] = self.shock(self.x0).T
        self.qhat0[:,1:2,:] = 0
        qhat1 = np.zeros((self.nVar,self.nDOFs,self.iMAX))
        self.qhat1 = qhat1
        Lh0 = np.zeros((self.nVar,self.nDOFs,self.iMAX))
        self.Lh0 = Lh0
        Lh1 = np.zeros((self.nVar,self.nDOFs,self.iMAX))
        self.Lh1 = Lh1
        phiL = self.phis(0)
        phiR = self.phis(1)
        self.phiL = phiL
        self.phiR = phiR
        FL, FR = np.empty(self.qhat0.shape), np.empty(self.qhat0.shape)
        self.FL = FL
        self.FR = FR

    def phis(self, xi):
        """
        function y=phis(xi)
            y(1) = 1;
            y(2) = 2*xi-1;
            y(3) = 1-6*xi+6*xi^2;"""
        return np.array([1.,
                         2. * xi - 1.,
                         1. - 6. * xi + 6. * xi**2.])


    def get_qlqr(self, phiL, phiR, qhat):
        """
        for i=1:IMAX
            % Compute left and right boundary cell values
            phiL = phis(0);
            phiR = phis(1);
            for j=1:2
                ql(j,i) = phiL(1)*qhat(j,1,i)+phiL(2)*qhat(j,2,i)+phiL(3)*qhat(j,3,i);
                qr(j,i) = phiR(1)*qhat(j,1,i)+phiR(2)*qhat(j,2,i)+phiR(3)*qhat(j,3,i);
            end
        end
        """
        ql, qr = np.empty((self.iMAX, 2)), np.empty((self.iMAX, 2))
        # non si può calcolare direttamente?
        for i in range(self.iMAX):
            # Compute left and right boundary cell values
            for j in (0,1):
                #TODO: use numpy casting
                # compute left
                ql[i,j] =  phiL[0] * qhat[j,0,i]\
                          +phiL[1] * qhat[j,1,i]\
                          +phiL[2] * qhat[j,2,i]
                # compute right
                qr[i,j] =  phiR[0] * qhat[j,0,i]\
                          +phiR[1] * qhat[j,1,i]\
                          +phiR[2] * qhat[j,2,i]
            #import pdb; pdb.set_trace()
        return ql, qr

    def get_qlqr2(self, phiL, phiR, qhat):
        """
        for i=1:IMAX
            % Compute left and right boundary cell values
            phiL = phis(0);
            phiR = phis(1);
            for j=1:2
                ql(j,i) = phiL(1)*qhat(j,1,i)+phiL(2)*qhat(j,2,i)+phiL(3)*qhat(j,3,i);
                qr(j,i) = phiR(1)*qhat(j,1,i)+phiR(2)*qhat(j,2,i)+phiR(3)*qhat(j,3,i);
            end
        end
        """
        #import pdb; pdb.set_trace()
        muL = qhat*np.array((phiL[:, np.newaxis],phiL[:,np.newaxis]))
        muR = qhat*np.array((phiR[:, np.newaxis],phiR[:,np.newaxis]))
        return np.sum(muL, axis=1).T, np.sum(muR, axis=1).T

    def get_FLFR(self, phiL, phiR, ql, qr):
        """
        for i=1:IMAX
            phiL = phis(0);
            phiR = phis(1);
            % Left and right values for base functions
            if(i==1)
                fluL = A*ql(:,i);
                fluR = A*ERS(qr(:,i),ql(:,i+1),A,0);
            elseif(i==IMAX)
                fluR = A*qr(:,i-1);
                fluL = A*ERS(qr(:,i-1),ql(:,i),A,0);
            else
                fluR = A*ERS(qr(:,i),ql(:,i+1),A,0);
                fluL = A*ERS(qr(:,i-1),ql(:,i),A,0);
            end
            for j=1:3
                FR(:,j,i) = phiR(j)*fluR;
                FM(:,j,i) = phiL(j)*fluL;
            end
        end
        """
        FL, FR = np.empty(self.qhat0.shape), np.empty(self.qhat0.shape)
        for i in range(self.iMAX):
            # Left and right values for base functions
            if i==0: # the first
                fluL = np.dot(self.A, ql[i, :])
                fluR = np.dot(self.A, self.ers(qr[i  ,:], ql[i+1,:], self.A, 0))
            elif i==self.iMAX-1: # the last
                fluL = np.dot(self.A, self.ers(qr[i-1,:], ql[i  ,:], self.A, 0))
                fluR = np.dot(self.A, qr[i-1, :])
            else:
                fluL = np.dot(self.A, self.ers(qr[i-1,:], ql[i  ,:], self.A, 0))
                fluR = np.dot(self.A, self.ers(qr[i  ,:], ql[i+1,:], self.A, 0))

            for j in (0,1,2): #TODO: use numpy vectors:

                FL[:, j, i] = phiL[j]*fluL
                FR[:, j, i] = phiR[j]*fluR
            #import pdb; pdb.set_trace()
            #phiL * fluL[:,np.newaxis]
        return FL, FR

    def get_FLFR2(self, phiL, phiR, ql, qr):
        # i==0
        fluL = np.dot(self.A, ql[0, :])
        fluR = np.dot(self.A, self.ers(qr[0  ,:], ql[1,:], self.A, 0))
        self.FL[:, :, 0] = phiL*fluL[:,np.newaxis]
        self.FR[:, :, 0] = phiR*fluR[:,np.newaxis]
        for i in range(1,self.iMAX-1):
            fluL = np.dot(self.A, self.ers(qr[i-1,:], ql[i  ,:], self.A, 0))
            fluR = np.dot(self.A, self.ers(qr[i  ,:], ql[i+1,:], self.A, 0))
            self.FL[:, :, i] = phiL*fluL[:,np.newaxis]
            self.FR[:, :, i] = phiR*fluR[:,np.newaxis]
        # i==-1
        fluL = np.dot(self.A, self.ers(qr[-2,:], ql[-1  ,:], self.A, 0))
        fluR = np.dot(self.A, qr[-2, :])
        self.FL[:, :, -1] = phiL*fluL[:,np.newaxis]
        self.FR[:, :, -1] = phiR*fluR[:,np.newaxis]
        return self.FL, self.FR


    def lhgen(self, qhat):
        """"
        function Lh = Lhgen(IMAX,qhat,iM,K,A)
        [=> get_qlqr]
        [=> get_FLFR]
        %
        for i=1:IMAX
            for p=1:3
                 flus(:,p) = A*qhat(:,p,i);
            end
            for j=1:2
                Lh(j,:,i) = iM*((K*flus(j,:)')'-FR(j,:,i)+FM(j,:,i))';
            end
        end
        """
        ql, qr = self.get_qlqr(self.phiL, self.phiR, qhat)
        FL, FR = self.get_FLFR(self.phiL, self.phiR, ql, qr)
        flus = np.empty((2,3))
        Lh = np.empty(qhat.shape)
        for i in range(self.iMAX):
            for p in (0, 1, 2):
                flus[:,p] = np.dot(self.A, qhat[:,p,i])
            for j in (0, 1):
                Lh[j, :, i] = np.dot(self.iM,
                                    (np.dot(self.K, flus[j,:].T)\
                                    - FR[j,:,i] + FL[j, :, i]).T)
        return Lh

    def lhgen2(self, qhat):
        """"
        function Lh = Lhgen(IMAX,qhat,iM,K,A)
        [=> get_qlqr]
        [=> get_FLFR]
        %
        for i=1:IMAX
            for p=1:3
                 flus(:,p) = A*qhat(:,p,i);
            end
            for j=1:2
                Lh(j,:,i) = iM*((K*flus(j,:)')'-FR(j,:,i)+FM(j,:,i))';
            end
        end
        """
        ql, qr = self.get_qlqr2(self.phiL, self.phiR, qhat)
        #import pdb; pdb.set_trace()
        FL, FR = self.get_FLFR2(self.phiL, self.phiR, ql, qr)
        flus = np.empty((2,3))
        Lh = np.empty(qhat.shape)
        for i in range(self.iMAX):
            for p in (0, 1, 2):
                flus[:,p] = np.dot(self.A, qhat[:,p,i])
            for j in (0, 1):
                Lh[j, :, i] = np.dot(self.iM,
                                    (np.dot(self.K, flus[j,:].T)\
                                    - FR[j,:,i] + FL[j, :, i]).T)
        return Lh


    def dgrk(self, itime, dt):
        # TVD Runge-Kutta (Jiang and Shu)
        Lh0 = self.lhgen(self.qhat0)
        #import pdb; pdb.set_trace()
        qhat1 = self.qhat0 + self.dt*Lh0
        Lh1 = self.lhgen(qhat1)
        qhat2 = 3./4. * self.qhat0 + 1./4. * qhat1 + 1./4. * self.dt*Lh1
        Lh2 = self.lhgen(qhat2)
        self.qhat0 = 1./3. * self.qhat0 + 2./3. * qhat2 + 2./3. * self.dt*Lh2;
        return np.array([self.qhat0[0,0,:], self.qhat0[1,0,:]]).T

    def dgrk2(self, itime, dt):
        # TVD Runge-Kutta (Jiang and Shu)
        Lh0 = self.lhgen2(self.qhat0)
        #import pdb; pdb.set_trace()
        qhat1 = self.qhat0 + self.dt*Lh0
        Lh1 = self.lhgen2(qhat1)
        qhat2 = 3./4. * self.qhat0 + 1./4. * qhat1 + 1./4. * self.dt*Lh1
        Lh2 = self.lhgen2(qhat2)
        self.qhat0 = 1./3. * self.qhat0 + 2./3. * qhat2 + 2./3. * self.dt*Lh2;
        return np.array([self.qhat0[0,0,:], self.qhat0[1,0,:]]).T








#
# TEST
#

## Solver1D
#problem = Solver1D(advection_speed = -1)
#problem.compute('Lax Wendrof')
#problem.show()
#problem.compute('Lax Friedrich')
#problem.show()
#problem.compute('Warming Beam')
#problem.show()
#problem.compute('Fromm')
#problem.show()

A0 = [[0, 2],[ 1, 0]]
R0 = [[-np.sqrt(2), np.sqrt(2)],[1., 1.]]
D0 = [[-np.sqrt(2), 0],[0., np.sqrt(2)]]

A1 = [[0, 1],[ 2, 0]]
R1 = [[-np.sqrt(2), np.sqrt(2)], [1, 1]]
D1 = [[-np.sqrt(2),0], [0, np.sqrt(2)]]

# Solver1Dsys
#problem = Solver1Dsys(iMAX=200, courant_numb=0.9)
#problem.compute_numerical('Lax Wendrof')
#problem.show2()
#problem = Solver1Dsys(advection_speed = [[0, 1],[ 2, 0]],
#                      R = [[-np.sqrt(2), np.sqrt(2)], [1, 1]],
#                      D = [[-np.sqrt(2),0], [0, np.sqrt(2)]],
#                      qL = [1,0], qR = [0,1],
#                      iMAX=200, courant_numb=0.9)
#problem.compute_numerical('Lax Friedrich')
#problem.show2()
#problem = Solver1Dsys(advection_speed = A0,
#                      R = R0,
#                      D = D0,
#                      qL = [2,1], qR = [1,0],
#                      iMAX=200, courant_numb=0.9)
#problem.compute_numerical('Warming Beam')
#problem.show2()
#problem.compute_numerical('Fromm')
#problem.show2()
#problem.compute_numerical('WENO')
#problem.show2()
#problem.compute_numerical('DGRK')
#problem.show2()
# RuntimeError: underlying C/C++ object has been deleted
#dgrk = DGRK(iMAX = 20, accuracy = 2,
#            M = [(1, 0, 0),( 0, 1./3., 0),(0, 0, 1./5.)],
#            K = [(0, 0, 0),(2, 0, 0),(0, 2, 0)])
#dgrk.compute_numerical('DGRK')
#dgrk.show2()
#dgrk = DGRK(iMAX = 20, accuracy = 2,
#            M = [(1, 0, 0),( 0, 1./3., 0),(0, 0, 1./5.)],
#            K = [(0, 0, 0),(2, 0, 0),(0, 2, 0)])
#dgrk.compute_numerical('DGRK2')
#dgrk.show2()
#import pdb; pdb.set_trace()

