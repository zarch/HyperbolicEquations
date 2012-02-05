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
                      'Weno'          : self.weno,}
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
        #import pdb; pdb.set_trace()


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
problem = Solver1Dsys(iMAX=200, courant_numb=0.9)
problem.compute_numerical('Lax Wendrof')
problem.show2()
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
#problem.compute_numerical('Weno')
#problem.show2()
# RuntimeError: underlying C/C++ object has been deleted
