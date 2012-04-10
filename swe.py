# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 09:34:27 2012

@author: pietro
"""

import numpy as np
import math
from num import Solver1Dsys
import matplotlib.delaunay as delaunay
#import matplotlib.pyplot as plt
#from matplotlib import cm
#from mpl_toolkits.mplot3d.axes3d import get_test_data
from vtk import exportVTK
import os

import datetime


class Swe1D(Solver1Dsys):
    def __init__(self, advection_speed = np.array([[0., 2.], [ 1., 0.]]),
                 R = np.array([[-np.sqrt(2), np.sqrt(2)],[1., 1.]]),
                 D = np.array([[-np.sqrt(2), 0],[0., np.sqrt(2)]]),
                 x0 = 0, xL = -1, xR = 1,
                 tSTART = 0, tEND = 0.2,
                 iMAX = 3, nMAX = 3, courant_numb=0.5,
                 qL = np.array([2., 1.]), qR = np.array([1., 0.]),
                 n_iter = 500, tol = 1e-6):
        Solver1Dsys.__init__(self, advection_speed = advection_speed,
                 R = R, D = D, x0 = x0, xL = xL, xR = xR,
                 tSTART = tSTART, tEND = tEND,
                 iMAX = iMAX, nMAX = nMAX, courant_numb=courant_numb,
                 qL = qL, qR = qR)
        # add the new methods
        self.mthds['Godunov'] = self.godunov
        self.dt = 2 #FIXME: capire meglio
        self.n_iter = n_iter
        self.tol = tol

    def get_dt(self, q):
        """
        dt = 2;
        for i=1:IMAX
            dt = min(dt,c*dx/max(max(abs(eigSWE(q(:,i))))));
        end
        """
        #dt1 = self.dt
        dt = self.courant * self.dx / np.abs(self.eig(q)).max()
        return dt

    def eig(self, q):
        """
        function lambda = eigSWE(q)

        c = sqrt(9.81*q(1));

        lambda(1)=q(2)/q(1)-c;
        lambda(2)=q(2)/q(1)+c;
        """
        q1, q2 = q.T
        c = np.sqrt(9.81*q1)
        lambda1 = q2/q1 - c
        lambda2 = q2/q1 + c
        return np.array([lambda1, lambda2])


    def wavetype(self, hS, uS, cS, hX, uX, cX, xt, sign):
        """
        # Sample left wave
        if(hS >= hL)
            % Left shock
            qL = sqrt((hS + hL)*hS/(2.0*hL*hL));
            sL = uL - cL*qL;
            if(s <= sL)
                % Sample point lies to the left of the shock
                h = hL;
                u = uL;
            else
                % Sample point lies to the right of the shock
                h = hS;
                u = uS;
            end
         else
            % Left rarefaction
            shL = uL - cL;
            if(s <= shL)
                % Sample point lies to the right of the rarefaction
                h = hL;
                u = uL;
            else
                stL = uS - cS;
                if(s <= stL)
                    % Sample point lies inside the rarefaction
                    u = (uL + 2.0*cL + 2.0*s)/3.0;
                    c = (uL + 2.0*cL - s)/3.0;
                    h = c*c/gravit;
                else
                    % Sample point lies in the STAR region
                    h = hS ;
                    u = uS;
                end
            end
        end
        """
        if hS >= hX:
            # Shock
            qX = math.sqrt( (hS + hX) * hS / (2. * hX * hX) )
            sX = uX + sign * cX * qX
            if xt <= sX:
                # Sample point lies to the left of the shock
                if sign==-1: return hX, uX
                else: return hS, uS
            else:
                # Sample point lies to the right of the shock
                if sign==-1: return hS, uS
                else: return hX, uX
        else:
            # Rarefaction
            shX =  uX + sign * cX
            if xt <= shX:
                # Sample point lies to the right of the rarefaction
                return hX, uX
            else:
                stX = uS + sign * cS
                if xt <= stX:
                    # Sample point lies inside the rarefaction
                    u = (uX + 2.*cX + 2.*xt)/3.
                    c = (uX - sign * 2. * cX + sign * xt)/3.
                    h = c * c / 9.81
                    return h, u
                else:
                    # Sample point lies in the STAR region
                    return hS, uS




    def samwet(self, hL, hR, hS, uL, uR, uS, cL, cR, cS, xt):
        """
        function [h,u] = samwet(hL,hR,hS,uL,uR,uS,cL,cR,cS,s)

        gravit = 9.81;
        if(s <= uS)
            % Sample left wave
            if(hS >= hL)
                % Left shock
                qL = sqrt((hS + hL)*hS/(2.0*hL*hL));
                sL = uL - cL*qL;
                if(s <= sL)
                    % Sample point lies to the left of the shock
                    h = hL;
                    u = uL;
                else
                    % Sample point lies to the right of the shock
                    h = hS;
                    u = uS;
                end
             else
                % Left rarefaction
                shL = uL - cL;
                if(s <= shL)
                    % Sample point lies to the right of the rarefaction
                    h = hL;
                    u = uL;
                else
                    stL = uS - cS;
                    if(s <= stL)
                        % Sample point lies inside the rarefaction
                        u = (uL + 2.0*cL + 2.0*s)/3.0;
                        c = (uL + 2.0*cL - s)/3.0;
                        h = c*c/gravit;
                    else
                        % Sample point lies in the STAR region
                        h = hS ;
                        u = uS;
                    end
                end
            end
        else
            % Sample right wave
            if(hS >= hR)
                % Right shock
                qR = sqrt((hS + hR)*hS/(2.0*hR*hR));
                sR = uR + cR*qR;
                if(s >= sR)
                    % Sample point lies to the right of the shock
                    h = hR;
                    u = uR;
                else
                    % Sample point lies to the left of the shock
                    h = hS;
                    u = uS;
                end
            else
                % Right rarefaction
                shR = uR + cR;
                if(s >= shR)
                    % Sample point lies to the right of the rarefaction
                    h = hR;
                    u = uR;
                else
                    stR = uS + cS;
                    if(s >= stR)
                        % Sample point lies inside the rarefaction
                        u = (uR  - 2.0*cR + 2.0*s)/3.0;
                        c = (-uR + 2.0*cR + s)/3.0;
                        h = c*c/gravit;
                    else
                        % Sample point lies in the STAR region
                        h = hS;
                        u = uS;
                    end
                end
            end
        end
        """
        if xt <= uS:
            # Sample left wave
            return self.wavetype(hS, uS, cS, hL, uL, cL, xt, -1.)
        else:
            # Sample right wave
            return self.wavetype(hS, uS, cS, hR, uR, cR, xt, 1.)


    def ers(self, qL, qR, A, xt):
        """
        function qE = ERSSWE(qL,qR,s)

        gravit = 9.81;
        niter = 1000;
        tol = 1e-6;
        % Compute left and right conserved variables, celerities and fluxes
        hL = qL(1); uL = qL(2)/qL(1);
        hR = qR(1); uR = qR(2)/qR(1);
        cL = sqrt(gravit*qL(1));
        cR = sqrt(gravit*qR(1));
        fL = fluxSWE(qL);
        fR = fluxSWE(qR);
        % Find starting value for Newton-Rhapson iterative method
        hS = starte(hL,hR,cL,cR,uL,uR);
        h0 = hS;
        % Start iteration
        for it=1:niter
            [fL,fld] = geofun(hS,hL,cL);
            [fR,frd] = geofun(hS,hR,cR);
            hS  = hS - (fL + fR + uR-uL)/(fld + frd);
            cha = abs(hS-h0)/(0.5*(hS+h0));
            if(cha <= tol)
                break
            end
            if(hS < 0.)
                hS = tol;
            end
            h0 = hS;
        end
        % Converged solution for depth hS in Star Region.
        % Compute velocity uS in Star Region
        uS = 0.5*(uL + uR) + 0.5*(fR - fL);
        cS = sqrt(gravit*hS);
        % Sample solution
        [h,u] = samwet(hL,hR,hS,uL,uR,uS,cL,cR,cS,s);
        qE(1)=h;
        qE(2)=h*u;
        """
        hL, hR, uL, uR, cL, cR = self.get_param(qL, qR)
        fL = self.flux(qL)
        fR = self.flux(qR)
        # Find starting value for Newton-Rhapson iterative method
        hS = self.start_exact(hL, hR, cL, cR, uL, uR)
        h0 = hS
        # Start iteration
        for _ in range(self.n_iter):
            fL, fld = self.geofun(hS, hL, cL)
            fR, frd = self.geofun(hS, hR, cR)
            hS  = hS - (fL + fR + uR - uL)/(fld + frd)
            if abs(hS - h0) / (0.5 * (hS + h0) ) <= self.tol: break
            if hS < 0.: hS = self.tol
            h0 = hS
        # Check if converged or not
        if _ == self.n_iter and abs(hS - h0) / (0.5 * (hS + h0) ) > self.tol:
            raise 'Solution not converged, try to increase n_iter value'
        else:
            # Converged solution for depth hS in Star Region.
            # Compute velocity uS in Star Region
            uS = 0.5*(uL + uR) + 0.5*(fR - fL);
            cS = math.sqrt(9.81*hS);
            # Sample solution
            h, u = self.samwet(hL, hR, hS, uL, uR, uS, cL, cR, cS, xt)
            return h, h*u

    def start_exact(self, hL, hR, cL, cR, uL, uR):
        """function hS = starte(hL,hR,cL,cR,uL,uR)

        gravit = 9.81;
        dMin = min(hL,hR);
        % Use Two-Rarefaction (TRRS) solution as starting value
        hS = (1/gravit)*(0.5*(cL+cR)-0.25*(uR-uL))^2;
        if(hS > dMin)
            % Use two-shock (TSRS) solution as starting value
            % with DS as computed from TRRS as estimate
            gel = sqrt(0.5*gravit*(hS+hL)/(hS*hL));
            ger = sqrt(0.5*gravit*(hS+hR)/(hS*hR));
            hS = (gel*hL+ger*hR-(uR-uL))/(gel+ger);
        end
        """
        dMin = min(hL, hR)
        # Use Two-Rarefaction (TRRS) solution as starting value
        hS = (1/9.81)*(0.5*(cL + cR) - 0.25 * ( uR - uL ) )**2
        if hS > dMin:
            # Use two-shock (TSRS) solution as starting value
            # with DS as computed from TRRS as estimate
            gel = math.sqrt(0.5*9.81 * (hS + hL)/(hS * hL))
            ger = math.sqrt(0.5*9.81 * (hS + hR)/(hS * hR))
            return (gel*hL+ger*hR-(uR-uL))/(gel+ger);
        return hS


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


    def geofun(self, h, hK, cK):
        """
        function [f,fD] = geofun(h,hK,cK)

        gravit = 9.81;
        if(h < hK)
            % Wave is rarefaction wave (or depression)
            c = sqrt(gravit*h);
            f = 2*(c-cK);
            fD = gravit/c;
        else
            % Wave is shock wave (or bore)
            ges = sqrt(0.5*gravit*(h+hK)/(h*hK));
            f = (h-hK)*ges;
            fD = ges - 0.25*gravit*(h-hK)/(ges*h*h);
        end
        """
        if (h < hK):
            # Wave is rarefaction wave (or depression)
            c = math.sqrt(9.81 * h) # we use math because is faster
            return 2 * (c-cK), 9.81/c
        else:
            # Wave is shock wave (or bore)
            ges = math.sqrt(0.5*9.81 * ( h + hK )/( h * hK ) )
            return ( h - hK ) * ges, ges - 0.25*9.81 * (h-hK)/(ges * h*h)


    def flux(self, q):
        """
        function f = fluxSWE(q)

        h = q(1);
        if(q(1) > 0)
            u = q(2)/q(1);
        else
            u = 0;
        end


        f(1) = h*u;
        f(2) = h*u^2+1/2*9.81*h^2;
        """
        q1, q2 = q
        if q1 > 0:
            u = q2/q1
        else:
            u = 0
        return np.array([q1*u, q1 * u*u + 0.5*9.81 * q1*q1])

    def get_param(self, qL, qR):
        """Compute left and right conserved variables, celerities and fluxes
        hL = qL(1); uL = qL(2)/qL(1);
        hR = qR(1); uR = qR(2)/qR(1);
        cL = sqrt(gravit*qL(1));
        cR = sqrt(gravit*qR(1));
        fL = fluxSWE(qL);
        fR = fluxSWE(qR);

        """
        hL, hR = qL[0],    qR[0]
        uL, uR = qL[1]/hL, qR[1]/hR
        cL, cR = np.sqrt(9.81*hL), np.sqrt(9.81*hR)
        return hL, hR, uL, uR, cL, cR,

    def fhll(self, qL, qR):
        """
        function f = FHLL(qL,qR)

        hL = qL(1);
        hR = qR(1);
        uL = qL(2)/qL(1);
        uR = qR(2)/qR(1);
        cL = sqrt(9.81*hL);
        cR = sqrt(9.81*hR);

        fL = [hL*uL; hL*uL^2+0.5*9.81*hL^2];
        fR = [hR*uR; hR*uR^2+0.5*9.81*hR^2];

        uS = 0.5*(uL+uR) + cL - cR;
        cS = 0.25*(uL-uR) + 0.5*(cL+cR);

        sL = min(uL-cL,uS-cS);
        sR = max(uR+cR,uS+cS);

        fS = (sR*fL-sL*fR+sL*sR*(qR-qL))/(sR-sL);

        if(sL>=0)
            f = fL;
        elseif(sL<0 && sR>0)
            f = fS;
        else
            f = fR;
        end
        """
        hL, hR, uL, uR, cL, cR = self.get_param(qL, qR)
        fL = np.array([hL*uL, hL*uL*uL+0.5*9.81*hL*hL])
        fR = np.array([hR*uR, hR*uR*uR+0.5*9.81*hR*hR])
        uS = 0.5*(uL+uR) + cL - cR
        cS = 0.25*(uL-uR) + 0.5*(cL+cR)
        sL = min(uL-cL,uS-cS);
        sR = max(uR+cR,uS+cS)
        fS = (sR*fL-sL*fR+sL*sR*(qR-qL))/(sR-sL)
        if sL>=0:
            return fL
        elif sL<0 and sR>0:
            return fS
        else:
            return fR


    def weno(self, itime, dt):
        """
        for i=2:IMAX-1
            fp = fluxSWE(ERSSWE(qm(:,i),qp(:,i+1),0))';
            fm = fluxSWE(ERSSWE(qm(:,i-1),qp(:,i),0))';
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
            fp = self.flux(self.ers(qm[i,:],  qp[i+1,:],self.A,0))
            fm = self.flux(self.ers(qm[i-1,:],qp[i,:],self.A,0))
            q1[i,:] = q0[i,:] - dt/self.dx * ( fp - fm )
        # set boundary conditions
        q1[self.boundary,:] = q0[self.boundary,:]
        return q1

    def godunov(self, itime, dt):
        """
        % Godunov 1st order
        for i=2:IMAX-1
    %         fp = fluxSWE(ERSSWE(q(:,i),q(:,i+1),0))';
            fp = FHLL(q(:,i),q(:,i+1));
    %         fm = fluxSWE(ERSSWE(q(:,i-1),q(:,i),0))';
            fm = FHLL(q(:,i-1),q(:,i));
            q1(:,i) = q(:,i)-dt/dx*(fp-fm);
        end
        q1(:,1) = qL;    q1(:,2) = qL;       q1(:,3) = qL;
        q1(:,IMAX)= qR;  q1(:,IMAX-1) = qR;  q1(:,IMAX-2) = qR;
        q = q1;
        """
        q0 = self.dat[itime]['sol']
        #try:
        q1 = self.dat[itime+1]['sol']
        #except IndexError:
        #    import pdb; pdb.set_trace()
        for i in range(2, len(q0)-2):
            fp = self.fhll(q0[i  ,:], q0[i+1,:])
            fm = self.fhll(q0[i-1,:], q0[i  ,:])
            q1[i,:] = q0[i,:] - dt/self.dx * ( fp - fm )
        # set boundary conditions
        q1[self.boundary,:] = q0[self.boundary,:]
        return q1

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

    def set_initial_condition(self):
        return self.shock(self.x0)

    def compute_numerical(self, funct):
        # set intial condition
        print("Start compute using: {mth}".format(mth=funct))
        self.used = funct
        start = datetime.datetime.now()
        self.dat[0] = self.tSTART, self.set_initial_condition()
        #import pdb; pdb.set_trace()
        for itime in range(0, self.nMAX):
            time = self.dat['time'][itime]
            dt = self.get_dt(self.dat[itime]['sol']) #TODO: add this row to the Solver1Dsys method
            if  time + dt > self.tEND:
                dt = self.tEND - time
            self.dat[itime+1] = time + dt, self.mthds[funct](itime, dt)
            if time + dt >= self.tEND: break
        end = datetime.datetime.now()
        self.last_time = itime
        print("Time necessary to compute: {0}".format(end - start))


def grid(xvect, yvect):
    points = []
    for row in xvect:
        for col in yvect:
            points.append((row, col))
    return points



class Swe2D(Swe1D):
    def __init__(self, advection_speed = np.array([[0., 2.], [ 1., 0.]]),
                 R = np.array([[-np.sqrt(2), np.sqrt(2)],[1., 1.]]),
                 D = np.array([[-np.sqrt(2), 0],[0., np.sqrt(2)]]),
                 x0 = 0, xL = 0, xR = 1, yL = 0, yR = 1,
                 tSTART = 0, tEND = 0.2,
                 iMAX = 3, nMAX = 3, courant_numb=0.45,
                 qL = np.array([2., 1.]), qR = np.array([1., 0.]),
                 n_iter = 500, tol = 1e-6,
                 jMAX = 3, tri = None, neigh = None, side=None): #FIX use **kargs

        Swe1D.__init__(self, advection_speed = advection_speed,
                 R = R, D = D, x0 = x0, xL = xL, xR = xR,
                 tSTART = tSTART, tEND = tEND,
                 iMAX = iMAX, nMAX = nMAX,
                 courant_numb=courant_numb,
                 qL = qL, qR = qR,
                 n_iter = n_iter, tol = tol)
        # Definition of the computational domain
        self.yL, self.yR = yL, yR
        self.tSTART, self.tEND = tSTART, tEND
        # Define the number of iteration
        self.courant = courant_numb
        self.nMAX = nMAX
        # Define the mesh size
        self.iMAX, self.jMAX = iMAX, jMAX
        xvect, dx = np.linspace(self.xL, self.xR, self.iMAX,
                                          retstep = True)
        yvect, dy = np.linspace(self.yL, self.yR, self.jMAX,
                                          retstep = True)
        self.xvect, self.yvect = xvect, yvect
        self.dx, self.dy = dx, dy
        # generates the points of a corresponding cartesian mesh
        xx, yy = np.meshgrid(self.xvect, self.yvect)
        self.xx, self.yy = xx, yy
        points = grid(self.xvect, self.yvect)
        self.points = np.array(points)
        #import pdb; pdb.set_trace()
        xpt, ypt = self.points.T
        self.xpt, self.ypt = xpt, ypt # FIXME define if it is possible tosimply assign directly
        #import pdb; pdb.set_trace()
        if (tri, neigh) == (None, None):
            circumctrs, edges, tri_pnts, tri_neighbs = delaunay.delaunay(xpt, ypt)
        else:
            edges, tri_pnts, tri_neighbs = None, tri, neigh
        self.edges, self.tri_pnts, self.tri_neighbs = edges, tri_pnts, tri_neighbs
        #import pdb; pdb.set_trace()
        side_length = np.zeros(tri_pnts.shape)
        self.side_length = side_length
        V = np.zeros(len(tri_pnts))
        lmin = np.zeros(len(tri_pnts))
        self.V, self.lmin = V, lmin
        self.DTYPESOL = np.dtype([
                             ('time', np.float64),
                             ('sol',  np.float64, (len(self.tri_pnts), 3))])
        dat = np.zeros((self.nMAX,), dtype = self.DTYPESOL)
        self.dat = dat
        barycenters = np.array(self.compute_barycenters())
        self.bary = barycenters
        normvect = self.compute_normalvect()
        self.normvect = normvect
        xpt_b, ypt_b  = self.bary.T
        circumctrs_b, edges_b, tri_pnts_b, tri_neighbs_b = delaunay.delaunay(xpt_b, ypt_b)
        self.edges_b, self.tri_pnts_b, self.tri_neighbs_b = edges_b, tri_pnts_b, tri_neighbs_b
        self.last_time = 0


    def get_dt(self, q):
        """
        amax = 0;
        for i=1:nElem
            lambda = PDEEigenvalues(q(:,i),g);
            amax = max(amax,max(abs(lambda)));
        end
        dt = cfl*min(lmin)/amax;
        """
        #import pdb; pdb.set_trace()
        #dt =  self.courant * self.dx / np.abs(self.pde_eigv(q)).max()
        return self.courant * min(self.lmin) / np.abs(self.pde_eigv(q)).max()

    def compute_barycenters(self):
        """Return the triangle barycenter area and inscribed circle diameter"""
        barycenters = []
        for i,tri in enumerate(self.tri_pnts):
            #import pdb; pdb.set_trace()
            barycenters.append(1./3.*self.points[tri].sum(axis=0))
            tri0, tri1, tri2 = tri
            a = np.linalg.norm(self.points[tri1]-self.points[tri0])
            b = np.linalg.norm(self.points[tri2]-self.points[tri1])
            c = np.linalg.norm(self.points[tri0]-self.points[tri2])
            #import pdb; pdb.set_trace()
            self.side_length[i,:] = a,b,c
            s = 0.5 * (a+b+c)
            self.V[i] = math.sqrt(s*(s-a)*(s-b)*(s-c));
            self.lmin[i] = 2*self.V[i]/s
        return np.array(barycenters)

    def compute_normalvect(self):
        """Return the triangle barycenter area and inscribed circle diameter"""
        normvect = np.zeros((len(self.tri_pnts),3,3))
        zvec = np.array([0, 0, 1])
        for itri, tri in enumerate(self.tri_pnts):
            #import pdb; pdb.set_trace()
            tri0, tri1, tri2 = tri
            x1,y1 = self.points[tri1]-self.points[tri0]
            v1 = np.array([x1,y1,0])
            x2,y2 = self.points[tri2]-self.points[tri1]
            v2 = np.array([x2,y2,0])
            x3,y3 = self.points[tri0]-self.points[tri2]
            v3 = np.array([x3,y3,0])
            v1 = v1/np.linalg.norm(v1)
            v2 = v2/np.linalg.norm(v2)
            v3 = v3/np.linalg.norm(v3)
            #import pdb; pdb.set_trace()
            normvect[itri,:,:] = np.cross(v1,zvec), np.cross(v2,zvec), np.cross(v3,zvec)
        #import pdb; pdb.set_trace()
        return normvect

    def set_initial_condition(self):
        """nVar = 3;
        X0 = [0.5;0.5];
        q = zeros(nVar,nElem);
        for i=1:nElem
            q(1,i) = 1+1*exp(-0.5*((XB(i,1)-X0(1))^2+(XB(i,2)-X0(2))^2)/0.1^2);
        end
        """
        X0 = np.array([0.5, 0.5])
        XB = self.bary
        q0 = 1 + np.exp(-0.5*(np.sum((XB-X0[np.newaxis])**2., axis=1))/0.1**2)
        q1 = np.zeros(q0.shape)
        #import pdb; pdb.set_trace()
        return np.array([q0, q1, q1]).T


    def get_avg_points(self):
        """%% Define averaged values over nodes for plotting purposes only
        triB = delaunay(XB(:,1),XB(:,2));
        trisurf(triB,XB(:,1),XB(:,2),q(1,:));
        """
        pass

    def pde_eigv(self, u):
        """
        %
        % This function computes the eigenvalues of the
        % Jacobian of the system
        %
        function lambda = PDEEigenvalues(u,g)
          %

          %
          c(1) = sqrt(g*u(1));
          %
          vel(1) = sqrt((u(2)/u(1))^2 + (u(3)/u(1))^2);
          %
          lambda(1,1) = vel(1) - c(1);   % Left-going acoustic wave
          lambda(2,1) = vel(1);            % Contact / entropy wave
          lambda(3,1) = vel(1) + c(1);   % Right-going acoustic wave
        """
        u0, u1, u2 = u.T
        c = np.sqrt(9.81*u0)
        vel = np.sqrt((u1/u0)**2 + (u2/u0)**2)
        return np.array([vel-c, vel, vel+c])

    def flux(self, u):
        """
        %
        % Compute the nonlinear flux as a function
        % of the vector of conserved variables u
        %
        function [fx,fy]=PDEFlux(u,g)
           %
           % Euler equations
           %
           % u(1,:) = rho
           % u(2,:) = rho*u
           % u(3,:) = rho*E
           %

           %
           fx(1,1) = u(2);
           fx(2,1) = u(1)*(u(2)/u(1))^2+0.5*g*(u(1))^2;
           fx(3,1) = u(2).*u(3)/u(1);
           %
           fy(1,1) = u(3);
           fy(2,1) = u(3)*u(2)/u(1);
           fy(3,1) = u(1)*(u(3)/u(1))^2+0.5*g*(u(1))^2;

           %
        """
        flu = np.zeros((3,2), dtype=np.float64)
        flu[0,0] = u[1]
        flu[1,0] = u[0] * (u[1]/u[0])**2 + 0.5 * 9.81*u[0]**2
        flu[2,0] = u[1] * u[2]/u[0] #FIXME attenzione che c'è il punto controllare se sono scalari o vettori'
        flu[0,1] = u[2]
        flu[1,1] = u[2] * u[1]/u[0]
        flu[2,1] = u[0] * (u[2]/u[0])**2 + 0.5 * 9.81*u[0]**2
        return flu

    def godunov(self, itime, dt):
        """
        q1 = q;
        for i=1:nElem
            qL = q1(:,i);
            [fL,gL] = PDEFlux(qL,g);
            for j=1:3
                iNeighbor = triNeigh(i,j);
                nx = triNvec(i,j,1);
                ny = triNvec(i,j,2);
                if(iNeighbor<=0)
                    qR = qL;
                    qR(2:3) = -qL(2:3);
                else
                    qR = q(:,iNeighbor);
                end
                [fR,gR] = PDEFlux(qR,g);
                sL = max(max(abs(PDEEigenvalues(qL,g))));
                sR = max(max(abs(PDEEigenvalues(qR,g))));
                smax = max(sL,sR);
                fij = 0.5*((fL+fR)*nx+(gL+gR)*ny) - 0.5*smax*(qR-qL);
                q1(:,i) = q1(:,i)-dt/V(i)*sideLength(i,j)*fij;
            end
        end
        q=q1
        """
        #if itime==2: import pdb; pdb.set_trace()
        #print '='*30
        #print itime, dt, self.dat[itime]['time']+dt
        q0 = self.dat[itime]['sol']
        q1 = np.copy(q0[:])
        for i in range(len(self.tri_pnts)):
            #import pdb; pdb.set_trace()
            qL = np.copy(q0[i,:])
            fL, gL = self.flux(qL).T
            for ji in (0, 1, 2):
                #if i == 2 and ji == 1: import pdb; pdb.set_trace()
                nx, ny, nz = self.normvect[i,ji]
                ineighbor = self.tri_neighbs[i,ji]
                if ineighbor == -1:
                    # if not neighbours mirroring the numbers
                    qR = np.copy(qL)
                    qR[1:] = -qL[1:]
                else:
                    qR = np.copy(q0[ineighbor,:])
                fR, gR = self.flux(qR).T
                sL = np.abs(self.pde_eigv(qL)).max()
                sR = np.abs(self.pde_eigv(qR)).max()
                smax = max(sL,sR)
                fij = 0.5 * (np.dot((fL + fR ), nx) + np.dot((gL + gR), ny))- 0.5 * smax * (qR - qL)
                #import pdb; pdb.set_trace()
                q1[i,:] = q1[i,:] - dt/self.V[i] * np.dot(self.side_length[i,ji], fij)
        #print q1
        return q1

    def show2(self):
        """
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib import cm
        from mpl_toolkits.mplot3d.axes3d import get_test_data

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        X = np.arange(-5, 5, 0.25)
        Y = np.arange(-5, 5, 0.25)
        X, Y = np.meshgrid(X, Y)
        R = np.sqrt(X**2 + Y**2)
        Z = np.sin(R)

        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
                linewidth=0, antialiased=False)
        fig.colorbar(surf, shrink=0.5, aspect=10)
        plt.show()
        """
        #zfactor = 4
        xb, yb = self.bary.T
        sol0 = self.dat[0]['sol'][:,0]
        triangles = self.tri_pnts_b
        import mayavi.mlab as mlab
        fig = mlab.figure(bgcolor = (0.1, 0.1, 0.1),
                          size = (1280, 800))
        @mlab.animate()
        def showdat():
            """Example from:
            http://github.enthought.com/mayavi/mayavi/tips.html#animating-a-visualization
            """
            # triangular_mesh see:
            # http://github.enthought.com/mayavi/mayavi/auto/mlab_helper_functions.html?highlight=triangular_mesh#mayavi.mlab.triangular_mesh
            img = mlab.triangular_mesh(xb, yb, sol0, triangles, scalars=sol0)
            #fig = mlab.gcf()
            ms = img.mlab_source
            for t, s in self.dat:
                # see: http://github.enthought.com/mayavi/mayavi/mlab_animating.html?highlight=animating
                ms.set(scalars=s[:,0])
                yield
        a = showdat()

    def exp2vtk(self, dirname):
        try:
            os.mkdir(dirname)
        except OSError:
            pass
        name = 0
        xpt_b, ypt_b  = self.bary.T
        triangles = self.tri_pnts_b
        pointsdata = {'prova' : [1,] * len(self.bary.T)}
        print "Start to export in vtk"
        for time, data in self.dat[:self.last_time+1]:
            #import pdb; pdb.set_trace()
            pointsdata = {'prova' : data[:,0]}
            points = np.array([xpt_b,ypt_b, data[:,0]]).T
            fname = os.path.join(dirname, 'swe{0:06d}.vtk'.format(name))
            exportVTK(fname, points, triangles, pointsdata,
                      description = str(time))
            name += 1



def get_sides(triangle):
    tri0, tri1, tri2 = triangle
    for side in (tri1, tri2), (tri2, tri0), (tri0, tri1):
        yield side

def check_triangles_mesh(triangles, points):
    xpt, ypt = points.T
    neighbors = []
    for iC, triC in enumerate(triangles):
        tri0, tri1, tri2 = triC
        #import pdb; pdb.set_trace()
        v1 = np.array([xpt[tri1] - xpt[tri0], ypt[tri1] - ypt[tri0], 0])
        v2 = np.array([xpt[tri2] - xpt[tri0], ypt[tri1] - ypt[tri0], 0])
        cp = np.cross(v1, v2)
        if cp[2]<0: print('Element flipped', triC); break
        neigh = np.array([-1, -1, -1])
        print('-'*30)
        print(iC)
        count = 0
        # Identify the neighbours
        for sideC in get_sides(triC): # for each side of the triangle
            for iD, triD in enumerate(triangles): #for all the triangles
                if np.all(triC == triD): continue     # exclude the triangle himself
                for sideD in get_sides(triD): # for each side of each triangle
                    if sideC[0] == sideD[1] and sideC[1] == sideD[0]:
                        neigh[count] = iD
                        count+=1;
                        #print(repr(sideC), repr(sideD), iD)
        neighbors.append(neigh)
    return neighbors






tri3 = np.array([
     [3 ,    2 ,    6],
     [2 ,    1  ,   4],
     [2 ,    5  ,   6],
     [8 ,    9  ,   6],
     [5 ,    4  ,   8],
     [6 ,    5  ,   8],
     [2 ,    4  ,   5],
     [4 ,    7  ,   8],]) -1

triNeigh3 = np.array([
     [0   ,  3  ,   0],
    [ 0  ,   0  ,   7],
    [ 7  ,   6  ,   1],
    [ 0  ,   0  ,   6],
     [7  ,   8  ,   6],
     [3  ,   5  ,   4],
    [ 2   ,  5  ,   3],
    [ 0   ,  0  ,   5],]) -1

sideLength3 = np.array([
  [ 5.0000e-01,   7.0711e-01,   5.0000e-01],
  [ 5.0000e-01,   5.0000e-01,   7.0711e-01],
  [ 5.0000e-01,   5.0000e-01,   7.0711e-01],
   [5.0000e-01,   5.0000e-01,   7.0711e-01],
   [5.0000e-01,   7.0711e-01,   5.0000e-01],
   [5.0000e-01,   5.0000e-01,   7.0711e-01],
   [7.0711e-01,   5.0000e-01,   5.0000e-01],
   [5.0000e-01,   5.0000e-01,   7.0711e-01],])


tri4 = np.array([
     [8,     4,     7],
    [ 3 ,    2 ,    6],
   [  2  ,   1  ,   5],
   [ 10   ,  5   ,  9],
   [  7    , 6    ,10],
   [  2     ,5     ,6],
   [  7,     3,     6],
   [ 12 ,    8 ,   11],
   [  4  ,   3  ,   7],
   [ 11   ,  8   ,  7],
   [ 11    ,10    ,15],
   [  6     ,5    ,10],
   [ 16,    11,    15],
   [  7,    10 ,   11],
   [ 16 ,   12  ,  11],
   [ 10  ,   9   , 13],
   [ 10   , 14    ,15],
   [ 10    ,13    ,14],]) -1

triNeigh4 = np.array([
     [0  ,   9,    10],
    [ 0   ,  6 ,    7],
    [ 0    , 0  ,   6],
   [ 12,     0   , 16],
    [ 7 ,   12    ,14],
    [ 3  ,  12,     2],
    [ 9   ,  2 ,    5],
    [ 0,    10  ,  15],
   [  0 ,    7   ,  1],
   [  8  ,   1    ,14],
   [ 14   , 17    ,13],
   [  6    , 4,     5],
   [ 15,    11 ,    0],
   [  5 ,   11  ,  10],
   [  0  ,   8   , 13],
    [ 4   ,  0    ,18],
    [18    , 0,    11],
   [ 16     ,0 ,   17],]) -1

#problem = Swe2D(xL = 0, xR = 1, yL = 0, yR = 1, tEND=0.2,
#                iMAX = 20, jMAX = 20, nMAX=10000,)
#                #tri = tri4, neigh = triNeigh4)
#problem.compute_numerical('Godunov')
#problem.exp2vtk('provaVTKmatlab')

#bary = problem.compute_barycenters()








