#!/usr/bin/python

'''
* This file is part of IDHOA software.
*
* Copyright (C) 2013 Barcelona Media - www.barcelonamedia.org
* (Written by Davide Scaini <davide.scaini@barcelonamedia.org> for Barcelona Media)
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>,
* or write to the Free Software Foundation, Inc.,
* 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
'''



import sys
import numpy as np
from scipy import special, linalg
import math as mh

from __builtin__ import *


## checks
def controls(NSPK,DEG):
    if NSPK < (DEG+1)**2:
        sys.exit("\n##!!!!!!!!!!!!!!!!!!!!!!!!!!##\nYou don't have enough speakers for this decoding order.")




def gauss_legendre(n):
    # http://www.scientificpython.net/1/post/2012/04/gausslegendre1.html
    k=np.arange(1.0,n)       
    a_band = np.zeros((2,n)) 
    a_band[1,0:n-1] = k/np.sqrt(4*k*k-1) 
    x,V = linalg.eig_banded(a_band,lower=True) 
    w=2*np.real(np.power(V[0,:],2)) 
    return x, w

def maxReG1(order):
    x,w = gauss_legendre(order+1)
    return max(x)

def maxReCoeffs(order):
    coeffs = [special.legendre(i)(maxReG1(order)) for i in range(order+1)]
    return coeffs




def ambisoniC(phi,theta,DEC,DEG,inversion):
    NUM = len(phi)         # number of speakers
    
    if nD == '3D':
        #things to be initializated
        g0 = [None]*(DEG+1)
        g1 = [None]*(DEG+1)
        W  = [None]*(len(phi))
        if DEG>=1:
            X  = [None]*(len(phi))
            Y  = [None]*(len(phi))
            Z  = [None]*(len(phi))
        if DEG>=2:
            V  = [None]*(len(phi))
            T  = [None]*(len(phi))
            R  = [None]*(len(phi))
            S  = [None]*(len(phi))
            U  = [None]*(len(phi))
        if DEG>=3:
            Q  = [None]*(len(phi))
            O  = [None]*(len(phi))
            M  = [None]*(len(phi))
            K  = [None]*(len(phi))
            L  = [None]*(len(phi))
            N  = [None]*(len(phi))
            P  = [None]*(len(phi))
        if DEG>=4:
            c16  = [None]*(len(phi))
            c17  = [None]*(len(phi))
            c18  = [None]*(len(phi))
            c19  = [None]*(len(phi))
            c20  = [None]*(len(phi))
            c21  = [None]*(len(phi))
            c22  = [None]*(len(phi))
            c23  = [None]*(len(phi))
            c24  = [None]*(len(phi))
        if DEG>=5:
            c25  = [None]*(len(phi))
            c26  = [None]*(len(phi))
            c27  = [None]*(len(phi))
            c28  = [None]*(len(phi))
            c29  = [None]*(len(phi))
            c30  = [None]*(len(phi))
            c31  = [None]*(len(phi))
            c32  = [None]*(len(phi))
            c33  = [None]*(len(phi))
            c34  = [None]*(len(phi))
            c35  = [None]*(len(phi))
        if DEG>=6:
            sys.exit("DEG =",DEG," is not implemented yet\n")
    
    
        #####################################################
        ## Calculating the decoding dependent coefficients ##
        #####################################################
        if DEC == 'basic':
            for i in range(0,DEG+1):
                g0[i] = 1.
                g1[i] = 1.
    
        elif DEC == 'maxRe':
            g1p = maxReCoeffs(DEG) 
            Egm = sum([(2.*i+1.)*g1p[i]**2 for i in range(DEG+1)])
            g0 = [np.sqrt(NUM/Egm)]*4
            g1 = [g1p[i]*g0[0] for i in range(len(g1p))]
    
    
        elif DEC == 'phase':
            g0d = np.sqrt(3. *NUM/4.)            # from Dani 1st order Ambisonics
            g1d = g0d*1./3.                    # from Dani 1st order Ambisonics
            g1p = [None]*(DEG+1)
    
            for i in range(0,DEG+1):
                g0[i] = np.sqrt(NUM*(2*DEG+1)) / (DEG+1)
                g1[i] = g0[i]*(mh.factorial(DEG)*mh.factorial(DEG+1)) / float(mh.factorial(DEG+i+1)*mh.factorial(DEG-i))  # note that g1(0) is equal to g0
                g1p[i]= (mh.factorial(DEG)*mh.factorial(DEG+1)) / float(mh.factorial(DEG+i+1)*mh.factorial(DEG-i))          # just for debugging purposes
                
    
        else:
            sys.exit('Decoding scheme unknow: ', DEC, ' Possible ones are: basic, maxRe, phase')
    
    
    
        ##########################################
        ## Calculating the "naive" coefficients ##
        ##########################################
        for i in range(0,len(phi)):
            if DEG>=0:
                W[i] = 1.0/NUM 
                # W has this sqrt(2) factor to take into account that W is not directly the pressure
    
            if DEG>=1: 
                # 1st order
                X[i] = np.sqrt(3)* np.cos(phi[i])*np.cos(theta[i])/NUM;    # factor sqrt(3) comes from sqrt(2m+1) - see tab.3.3 p.156 Jerome Daniel
                Y[i] = np.sqrt(3)* np.sin(phi[i])*np.cos(theta[i])/NUM;   # The decoding matrix is given for SN3D (aka 'semi-normalized') so to go to N3D
                Z[i] = np.sqrt(3)* np.sin(theta[i])/NUM;                # you have to use the alphas. 
    
            if DEG>=2:
                # 2nd order
                V[i] = np.sqrt(5.)* np.sqrt(3.)/2.*np.sin(2.*phi[i])*np.cos(theta[i])**2/NUM ;
                T[i] = np.sqrt(5.)* np.sqrt(3.)/2.*np.sin(phi[i])*np.sin(2.*theta[i])/NUM ;
                R[i] = np.sqrt(5.)* (3.*np.sin(theta[i])**2.-1.)/2./NUM ;   
                S[i] = np.sqrt(5.)* np.sqrt(3.)/2.*np.cos(phi[i])*np.sin(2.*theta[i])/NUM ;
                U[i] = np.sqrt(5.)* np.sqrt(3.)/2.*np.cos(2.*phi[i])*np.cos(theta[i])**2/NUM ; 
    
            if DEG>=3:
                # 3rd order
                Q[i] = np.sqrt(7.)* np.sqrt(5./8.)*np.sin(3.*phi[i])*np.cos(theta[i])**3/NUM;
                O[i] = np.sqrt(7.)* np.sqrt(15.)/2.*np.sin(2.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**2/NUM; 
                M[i] = np.sqrt(7.)* np.sqrt(3./8.)*np.sin(phi[i])*np.cos(theta[i])*(5.*np.sin(theta[i])**2-1.)/NUM; 
                K[i] = np.sqrt(7.)* np.sin(theta[i])*(5.*np.sin(theta[i])**2-3.)/2./NUM;
                L[i] = np.sqrt(7.)* np.sqrt(3./8.)*np.cos(phi[i])*np.cos(theta[i])*(5.*np.sin(theta[i])**2-1.)/NUM; 
                N[i] = np.sqrt(7.)* np.sqrt(15.)/2.*np.cos(2.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**2/NUM;  
                P[i] = np.sqrt(7.)* np.sqrt(5./8.)*np.cos(3.*phi[i])*np.cos(theta[i])**3/NUM;
    
            if DEG>=4:
                # 4th order
                c16[i] = np.sqrt(9.)* np.sqrt(35./2.)*3./8.* np.sin(4.*phi[i])*np.cos(theta[i])**4  /NUM;  # (4,-4)
                c17[i] = np.sqrt(9.)* np.sqrt(35.)*3./4.* np.sin(3.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**3  /NUM;  # (4,-3)
                c18[i] = np.sqrt(9.)* np.sqrt(5./2.)/4.*( -3.* np.sin(2.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**3 + 21.*np.sin(2.*phi[i])*np.sin(theta[i])**2*np.cos(theta[i])**2 ) /NUM;  # (4,-2)
                c19[i] = np.sqrt(9.)* np.sqrt(5.)/4.* ( 21.*np.sin(phi[i])*np.cos(theta[i])*np.sin(theta[i])**3 - 9.*np.sin(phi[i])*np.cos(theta[i])*np.sin(theta[i]) )  /NUM;  # (4,-1)
                c20[i] = np.sqrt(9.)* 3./64.*( 20.*np.cos(2.*theta[i]) - 35.*np.cos(4.*theta[i]) -9. )  /NUM;  # (4,-0)
                c21[i] = np.sqrt(9.)* np.sqrt(5.)/4.* ( 21.*np.cos(phi[i])*np.cos(theta[i])*np.sin(theta[i])**3 - 9.*np.cos(phi[i])*np.cos(theta[i])*np.sin(theta[i]) )  /NUM;  # (4, 1)
                c22[i] = np.sqrt(9.)* np.sqrt(5./2.)/4.*( -3.* np.cos(2.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**3 + 21.*np.cos(2.*phi[i])*np.sin(theta[i])**2*np.cos(theta[i])**2 )  /NUM;  # (4, 2)
                c23[i] = np.sqrt(9.)* np.sqrt(35.)*3./4.* np.cos(3.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**3  /NUM;  # (4, 3)
                c24[i] = np.sqrt(9.)* np.sqrt(35./2.)*3./8.* np.cos(4.*phi[i])*np.cos(theta[i])**4  /NUM;  # (4, 4)
    
            if DEG>=5:
                # 5th order
                c25[i] = np.sqrt(11.)* 3./16.*np.sqrt(77.)* np.sin(5.*phi[i])*np.cos(theta[i])**5  /NUM;  # (5,-5)
                c26[i] = np.sqrt(11.)* 3./8.*np.sqrt(385./2.)* np.sin(4.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**4  /NUM;  # (5,-4) 
                c27[i] = np.sqrt(11.)* np.sqrt(385.)/16.* ( 9.*np.sin(3.*phi[i])*np.cos(theta[i])**3*np.sin(theta[i])**2 )  /NUM;  # (5, -3)
                c28[i] = np.sqrt(11.)* np.sqrt(1155./2.)/4.* ( 3.*np.sin(2.*phi[i])*np.cos(theta[i])**2*np.sin(theta[i])**3 - np.sin(2.*phi[i])*np.cos(theta[i])**2*np.sin(theta[i]) )  /NUM;  # (5,-2)
                c29[i] = np.sqrt(11.)* np.sqrt(165./2.)/8.* ( np.sin(phi[i])*np.cos(theta[i]) - 7.*np.sin(phi[i])*np.cos(theta[i])*np.sin(theta[i])**2 + 21.*np.sin(phi[i])*np.cos(theta[i])*np.sin(theta[i])**4 )  /NUM;  # (5,-1)
                c30[i] = np.sqrt(11.)* np.sqrt(11.)/8.* ( 15.*np.sin(theta[i]) - 70.*np.sin(theta[i])**3 + 63.*np.sin(theta[i])**5 ) /NUM;  # (5, 0)
                c31[i] = np.sqrt(11.)* np.sqrt(165./2.)/8.* ( np.cos(phi[i])*np.cos(theta[i]) - 7.*np.cos(phi[i])*np.cos(theta[i])*np.sin(theta[i])**2 + 21.*np.cos(phi[i])*np.cos(theta[i])*np.sin(theta[i])**4 )  /NUM;  # (5, 1)
                c32[i] = np.sqrt(11.)* np.sqrt(1155./2.)/4.* ( 3.*np.cos(2.*phi[i])*np.cos(theta[i])**2*np.sin(theta[i])**3 - np.cos(2.*phi[i])*np.cos(theta[i])**2*np.sin(theta[i]) )  /NUM;  # (5, 2)
                c33[i] = np.sqrt(11.)* np.sqrt(385.)/16.* ( 9.*np.cos(3.*phi[i])*np.cos(theta[i])**3*np.sin(theta[i])**2 )  /NUM;  # (5, 3)
                c34[i] = np.sqrt(11.)* 3./8.*np.sqrt(385./2.)* np.cos(4.*phi[i])*np.sin(theta[i])*np.cos(theta[i])**4   /NUM;  # (5, 4)
                c35[i] = np.sqrt(11.)* 3./16.*np.sqrt(77.)* np.cos(5.*phi[i])*np.cos(theta[i])**5  /NUM;  # (5, 5)
    
            if DEG>6:
                print "DEG =",DEG," is not implemented yet\n"
    
     
    
    
    
    
        if DEG==1:
            coeffs = np.array([W, Y, Z, X])
        elif DEG==2:
            coeffs = np.array([W, Y, Z, X, V, T, R, S, U])
        elif DEG==3:
            coeffs = np.array([W, Y, Z, X, V, T, R, S, U, Q, O, M, K, L, N, P])
        elif DEG==4:
            coeffs = np.array([W, Y, Z, X, V, T, R, S, U, Q, O, M, K, L, N, P, c16, c17, c18, c19, c20, c21, c22, c23, c24])
        elif DEG==5:
            coeffs = np.array([W, Y, Z, X, V, T, R, S, U, Q, O, M, K, L, N, P, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, c28, c29, c30, c31, c32, c33, c34, c35])
     
        if inversion ==1:
            coeffs =  np.linalg.pinv(coeffs, rcond=1e-8).T /NUM
            coeffs[abs(coeffs)<1e-8] = 0. # because inversion gives some very small values somewhere
    	
    
        ### MULTIPLYING FOR THE SELECTED DECODING SCHEME
        coeff = np.empty(coeffs.shape,dtype=np.float64)
        g1[0] = g0[0]
        for i in range(DEG+1):
            for jj in range(i**2,(i+1)**2):
                coeff[jj]=g1[i]*coeffs[jj]
        


############################################################################################
    ################################
    ## horizontal Ambisonics only ##
    ################################
    if nD == '2D':
        #things to be initializated
        g0 = [None]*(DEG+1)
        g1 = [None]*(DEG+1)
        W  = [None]*(len(phi))
        if DEG>=1:
            X  = [None]*(len(phi))
            Y  = [None]*(len(phi))
        if DEG>=2:
            V  = [None]*(len(phi))
            U  = [None]*(len(phi))
        if DEG>=3:
            Q  = [None]*(len(phi))
            P  = [None]*(len(phi))
        if DEG>=4:
            c8  = [None]*(len(phi))
            c9  = [None]*(len(phi))
        if DEG>=5:
            c10  = [None]*(len(phi))
            c11  = [None]*(len(phi))
        if DEG>=6:
            sys.exit("DEG =",DEG," is not implemented yet\n")
 
 
        #####################################################
        ## Calculating the decoding dependent coefficients ##
        #####################################################
# TOFIX TOBEDONE FIXME 
        if DEC == 'basic':
            for i in range(0,DEG+1):
                g0[i] = 1.
                g1[i] = 1.
    
        elif DEC == 'maxRe':
            g1p = [np.cos(i*np.pi/(2*DEG+2)) for i in range(0,DEG+1)]
            g0 = [np.sqrt(float(NUM)/(DEG+1))]*2
            g1 = [g1p[i]*g0[0] for i in range(len(g1p))]
            
            #print "g0 ",g0
            #print "gmprime ",g1p
            #print "rE ",g1p[1]
    
    
        elif DEC == 'phase':
            g1p = [( mh.factorial(DEG)**2 / float(mh.factorial(DEG+i)*mh.factorial(DEG-i) ) ) for i in range(0,DEG+1)] # ok
            Eg1p= g1p[0]**2 + 2.0* sum([g1p[i]**2 for i in range(1,DEG+1) ]) 
            g0[0]  = np.sqrt(float(NUM)/Eg1p)
            g1 = [x * g0[0] for x in g1p]

            #print "g0 ",g0
            #print "gmprime ",g1p
            #print "rE ", 2.0*DEG/(2.0*DEG+1)
            #print "E ", NUM/g0[0]**2, " 4th expected ", 2.627
            #print "E ", NUM/g0[0]**2, " 3rd expected ", 2.310
            #print g1
    
        else:
            sys.exit('Decoding scheme unknow: ', DEC, ' Possible ones are: basic, maxRe, phase')
   



        ##########################################
        ## Calculating the "naive" coefficients ##
        ##########################################
        for i in range(0,len(phi)):
            if DEG>=0:
                W[i] = 1.0/NUM 
                # W has this sqrt(2) factor to take into account that W is not directly the pressure
    
            if DEG>=1: 
                # 1st order
                X[i] = np.sqrt(2.)* np.cos(phi[i])/NUM;
                Y[i] = np.sqrt(2.)* np.sin(phi[i])/NUM;
    
            if DEG>=2:
                # 2nd order
                V[i] = np.sqrt(2.)* np.cos(2.*phi[i])/NUM;
                U[i] = np.sqrt(2.)* np.sin(2.*phi[i])/NUM; 
    
            if DEG>=3:
                # 3rd order
                Q[i] = np.sqrt(2.)* np.cos(3.*phi[i])/NUM;
                P[i] = np.sqrt(2.)* np.sin(3.*phi[i])/NUM;
    
            if DEG>=4:
                # 4th order
                c8[i] = np.sqrt(2.)* np.cos(4.*phi[i])/NUM;
                c9[i] = np.sqrt(2.)* np.sin(4.*phi[i])/NUM;
    

            if DEG>=5:
                # 5th order
                c10[i] = np.sqrt(2.)* np.cos(5.*phi[i])/NUM;
                c11[i] = np.sqrt(2.)* np.sin(5.*phi[i])/NUM;
 
            if DEG>=6:
                # 6th order
                c12[i] = np.sqrt(2.)* np.cos(6.*phi[i])/NUM;
                c13[i] = np.sqrt(2.)* np.sin(6.*phi[i])/NUM;
 
            if DEG>=7:
                # 7th order
                c14[i] = np.sqrt(2.)* np.cos(7.*phi[i])/NUM;
                c15[i] = np.sqrt(2.)* np.sin(7.*phi[i])/NUM;
 
            if DEG>=8:
                # 8th order
                c16[i] = np.sqrt(2.)* np.cos(8.*phi[i])/NUM;
                c17[i] = np.sqrt(2.)* np.sin(8.*phi[i])/NUM;
 
            if DEG>=9:
                # 9th order
                c18[i] = np.sqrt(2.)* np.cos(9.*phi[i])/NUM;
                c19[i] = np.sqrt(2.)* np.sin(9.*phi[i])/NUM;


            if DEG>9:
                print "DEG =",DEG," is not implemented yet\n"
    


        if DEG==1:
            coeffs = np.array([W, Y, X])
        elif DEG==2:
            coeffs = np.array([W, Y, X, V, U])
        elif DEG==3:
            coeffs = np.array([W, Y, X, V, U, Q, P])
        elif DEG==4:
            coeffs = np.array([W, Y, X, V, U, Q, P, c9, c8])
        elif DEG==5:
            coeffs = np.array([W, Y, X, V, U, Q, P, c9, c8, c11, c10])
        elif DEG==6:
            coeffs = np.array([W, Y, X, V, U, Q, P, c9, c8, c11, c10, c13, c12])
        elif DEG==7:
            coeffs = np.array([W, Y, X, V, U, Q, P, c9, c8, c11, c10, c13, c12, c15, c14])
        elif DEG==8:
            coeffs = np.array([W, Y, X, V, U, Q, P, c9, c8, c11, c10, c13, c12, c15, c14, c17, c16])
        elif DEG==9:
            coeffs = np.array([W, Y, X, V, U, Q, P, c9, c8, c11, c10, c13, c12, c15, c14, c17, c16, c19, c18])
     
        if inversion ==1:
            coeffs =  np.linalg.pinv(coeffs, rcond=1e-8).T /NUM
            coeffs[abs(coeffs)<1e-8] = 0. # because inversion gives some very small values somewhere
    	
    
        ### MULTIPLYING FOR THE SELECTED DECODING SCHEME
        coeff = np.empty(coeffs.shape,dtype=np.float64)
        g1[0] = g0[0]
        for i in range(DEG+1):
            for jj in range(i*2-1,i*2+1):
                if jj >= 0:
                    coeff[jj]=g1[i]*coeffs[jj]
     
 
        # at the end, please return something ;)
    return coeff





##########################
## supporting functions ##
##########################

def Sij(coeffSpk, coeffDir,NSphPt):
    sij = np.dot(coeffSpk.T, coeffDir*NSphPt) # this will have the dimensions of NSPK*NPOINTS
    if sij.shape != (NSPK,NSphPt): sys.exit("Wrong dimensions in Sij\n")
    return sij



def physOmni(Sij):
    # LOW FREQUENCIES
    # pressure
    pressure = sum(Sij)
    # velocity
    V = velocity(Sij)


    # HIGH FREQUENCIES
    # energy density
    energyD = sum(Sij*Sij)
    # intensity
    J = velocity(Sij*Sij)

    return pressure, V, energyD, J



def velocity(Sij):
    phi = np.asarray(PHI)
    theta = np.asarray(THETA)
    Sx = np.cos(phi) * np.cos(theta)
    Sy = np.sin(phi) * np.cos(theta)
    Sz = np.sin(theta)
    # rewrite these five lines using a call to ambisoniC (first order)
    # and fixing the coefficient of X,Y,Z

    if Sx.shape[0] == Sij.shape[0] and Sy.shape[0] == Sij.shape[0] and Sz.shape[0] == Sij.shape[0]:

        # since Sij and Sx are numpy arrays, the * is the element-wise multiplication
        # http://wiki.scipy.org/NumPy_for_Matlab_Users#head-e9a492daa18afcd86e84e07cd2824a9b1b651935
        Vx = sum((Sij.T * Sx).T) / sum(Sij)
        Vy = sum((Sij.T * Sy).T) / sum(Sij)
        Vz = sum((Sij.T * Sz).T) / sum(Sij)

    else: 
        sys.exit("Something went wrong calculating velocity\n")

    return Vx, Vy, Vz



def Radial(A,B):
    Sum = 0
    if len(A)!=len(B): sys.exit("I'm dying in Radial function. Arrays with different shapes.")
    for i in range(len(A)):
        Sum += A[i]*B[i]

    return Sum


def Tang(A,B):    # aka Cross
    Sum = 0
    if len(A)!=len(B): sys.exit("I'm dying in Tang function. Arrays with different shapes.")
    for i in range(len(A)):
        Sum += (A[i]*B[i-1]-A[i-1]*B[i])**2.0

    return np.sqrt(Sum)


def physDir(Sij,phi,theta):
    phi = np.asarray(phi)
    theta = np.asarray(theta)
    Zx = np.cos(phi) * np.cos(theta)
    Zy = np.sin(phi) * np.cos(theta)
    Zz = np.sin(theta)
    # rewrite these five lines using a call to ambisoniC (first order)
    # and fixing the coefficient of X,Y,Z
    
    Z = Zx, Zy, Zz

    press, V, energyD, J = physOmni(Sij)
    
    Vradial = Radial(V,Z)
    Jradial = Radial(J,Z)

    Vtangential = Tang(V,Z)
    Jtangential = Tang(J,Z)

    return press, V, energyD, J, Vradial, Jradial, Vtangential, Jtangential



def oppGain(Sij):
    # FOR IN-PHASE DECODING
    oppGain = np.zeros((NSPK,NPOINTS),dtype=np.float64)
    oppGain[Sij<0] = Sij[Sij<0]
    if oppGain.shape != (NSPK,NPOINTS): sys.exit("Wrong dimensions in oppGain\n")
    oppGain = sum(oppGain*oppGain)

    return oppGain




def sph2cart(phi,theta):
    """Acoustics convention!"""
    x = np.cos(phi) * np.cos(theta)
    y = np.sin(phi) * np.cos(theta)
    z = np.sin(theta)
    return x, y, z

def Physph2cart(phi,theta):
    """Physics convention!"""
    x = np.cos(phi) * np.sin(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(theta)
    return x, y, z

def PlotOverSphere(phi,theta,rho):
    """Acoustics convention!"""
    x = np.cos(phi) * np.cos(theta)*rho
    y = np.sin(phi) * np.cos(theta)*rho
    z = np.sin(theta)*rho
    return x, y, z


def eval_grad(f, theta):
    theta = algopy.UTPM.init_jacobian(theta)
    return algopy.UTPM.extract_jacobian(f(theta))

def eval_hess(f, theta):
    theta = algopy.UTPM.init_hessian(theta)
    return algopy.UTPM.extract_hessian(len(theta), f(theta))


def vtomat(vector):
    if nD == "2D":
        #tmp = np.reshape(vector,((DEG*2+1),NSPKmatch))
        tmp = vector.reshape((DEG*2+1),NSPKmatch)
    if nD == "3D":
        #tmp = np.reshape(vector,((DEG+1)**2,NSPKmatch))
        tmp = vector.reshape(((DEG+1)**2,NSPKmatch))
    return tmp


def mattov(mat):
    if nD == "2D":
        tmp = mat.reshape((1,(DEG*2+1)*NSPKmatch))[0]
    if nD == "3D":
        tmp = mat.reshape((1,(DEG+1)**2*NSPKmatch))[0]
    return tmp


def zeroing_bounds(nodes,upbound,lowbound,run):
    if nD == "3D":
        for i in range((DEG+1)**2*NSPKmatch):
            if np.asarray(abs(nodes)<3.*10e-4).reshape(1,((DEG+1)**2*NSPKmatch))[0,i]: upbound[i] = 0.; lowbound[i] = 0.; # putting to zero the speakers that are in a node of a SH

    if nD == "2D":
        for i in range((DEG*2+1)*NSPKmatch):
            if np.asarray(abs(nodes)<3.*10e-4).reshape(1,((DEG*2+1)*NSPKmatch))[0,i]: upbound[i] = 0.; lowbound[i] = 0.; # putting to zero the speakers that are in a node of a SH


    return nodes, upbound, lowbound


################################
# The function to be minimized #
################################
def function(VarCoeffSpk,*args):
    if len(args)>1 :
        grad = args[0]
        print "You have to implement the gradient or use another algorithm"

    shapeCoeff = np.shape(VarCoeffSpk)
    # this because the minimization library destroys the shape of VarCoeffSpk
    # in 3D
    if nD == "3D":
        if (shapeCoeff == ((DEG+1)**2*NSPKmatch,)):
            VarCoeffSpk = vtomat(VarCoeffSpk) 
        elif (shapeCoeff == ((DEG+1)**2*NSPK,)):
            sys.exit("Here the dimension is "+str(shapeCoeff)+" while it should be "+str(((DEG+1)**2*NSPKmatch,)) )
        elif (shapeCoeff != ((DEG+1)**2,NSPKmatch)) and  (shapeCoeff != ((DEG+1)**2,NSPK)):
            sys.exit("Strange dimensions of VarCoeffSpk in -function-."+str(shapeCoeff))


        if MATCHSPK: VarCoeffSpk = upscalematrix(VarCoeffSpk,MATCHED)
        if (np.shape(VarCoeffSpk) != ((DEG+1)**2,NSPK)):
            sys.exit("Here the matrix should be ((DEG+1)**2,NSPK)")


    # in 2D
    if nD == "2D":
        if (shapeCoeff == ((DEG*2+1)*NSPKmatch,)):
            VarCoeffSpk = vtomat(VarCoeffSpk) 
        elif (shapeCoeff == ((DEG*2+1)*NSPK,)):
            sys.exit("Here the dimension is "+str(shapeCoeff)+" while it should be "+str(( (DEG*2+1)*NSPKmatch,)) )
        elif (shapeCoeff != ((DEG*2+1),NSPKmatch)) and  (shapeCoeff != ((DEG*2+1),NSPK)):
            sys.exit("Strange dimensions of VarCoeffSpk in -function-."+str(shapeCoeff))


        if MATCHSPK: VarCoeffSpk = upscalematrix(VarCoeffSpk,MATCHED)
        if (np.shape(VarCoeffSpk) != ((DEG*2+1),NSPK)):
            sys.exit("Here the matrix should be ((DEG*2+1),NSPK)")



    """The function to be minimized"""
#    if (CP == None): CP = 400. # arbitrary coefficients
#    if (CE == None): CE = 400.
#    if (CR == None): CR = 100.
#    if (CT == None): CT = 100.
#    if (CPH == None):CPH= 1000.
#    if (CV == None): CV = 25.
    Wj = np.ones(coeffDir.shape[1]) # biasing factor dependent on direction j
    Mj = np.ones(coeffDir.shape[1]) # biasing factor dependent on direction j


    sij = Sij(VarCoeffSpk,coeffDir,NPOINTS)
    pressure, V, energyD, J, Vradial, Jradial, Vtang, Jtang = physDir(sij,phiTest,thetaTest)

    # weighting functions
    if WFRONT: Wj = Wj*WfrontVec # maybe it's better if you calculate Wfront Wplane and Wbinary outside the function, so that you calculate them only once...
    if WPLANE: Wj = Wj*WplaneVec
    if WBIN:   Wj = Wj*WbinVec
    if WAUTOREM: Mj = Mj*WremVec
    if WFRONT or WPLANE or WBIN: Wj = Wj*float(len(Wj))/sum(Wj)


    if DEC=='basic':
        Tpressure = ((1.-pressure)**2*Wj)/NPOINTS
        TVlon = ((1.-Vradial)**2*Wj*Mj)/NPOINTS
        TVtang = ((Vtang)**2*Wj*Mj)/NPOINTS
        Tvar = 0
        if PREFHOM:
            Tvar = np.var(VarCoeffSpk[0])/(np.mean(VarCoeffSpk[0]))**2
        
        target = np.sum(CP*Tpressure +  CR*TVlon + CT* TVtang) + CV*Tvar 
        # one sum instead of three (in Tpressure, TVlon,...)
        # anyway norm is (much) faster...

    elif DEC=='maxRe':
        TenergyD = ((1.-energyD)**2*Wj)/NPOINTS
        #TenergyD = ((1.-energyD)**2*Wj*Mj*2)/NPOINTS
        TJrad = ((1.-Jradial)**2*Wj*Mj)/NPOINTS
        TJtang = ((Jtang)**2*Wj*Mj)/NPOINTS
        Tvar = 0
        TJrad2 = np.ones(coeffDir.shape[1])
        if PREFHOM:
            Tvar = np.var(VarCoeffSpk[0])/(np.mean(VarCoeffSpk[0]))**2
        if (WAUTOREM or WBIN):
            TJrad2 = TJrad2*(~WremVec)*((Jradial)**2*Wj)/NPOINTS
            # ~array : negation of an array of bools
        
        target = np.sum(CE*TenergyD +  CR*TJrad + CT*TJtang + TJrad2) + CV*Tvar 
        
        #TenergyD = np.linalg.norm((1.-energyD)*Wj)**2/NPOINTS
        #TJrad = np.linalg.norm((1.-Jradial)*Wj*Mj)**2/NPOINTS
        #TJtang = np.linalg.norm((Jtang)*Wj*Mj)**2/NPOINTS
        #Tvar = 0
        #if PREFHOM:
        #    Tvar = np.var(VarCoeffSpk[0])/(np.mean(VarCoeffSpk[0]))**2
        #
        #target = CE*TenergyD + CR*TJrad + CT*TJtang + CV*Tvar


    elif DEC=='phase':
        TenergyD = ((1.-energyD)**2* Wj)/NPOINTS
        TJrad = ((1.-Jradial)**2*Wj*Mj)/NPOINTS
        TJtang = ((Jtang)**2*Wj*Mj)/NPOINTS
        Tvar = 0
        if PREFHOM:
            Tvar = np.var(VarCoeffSpk[0])/(np.mean(VarCoeffSpk[0]))**2
        ToppGain = np.linalg.norm(oppGain(sij))**2/NPOINTS

        target = np.sum(CE*TenergyD + CR*TJrad + CT*TJtang) + CPH*ToppGain * CV*Tvar  # missing some extra factors (see dani's paper page 5)

    return target


## with nlopt you can also write linear and non linear constraints... have a look at the reference
## http://ab-initio.mit.edu/wiki/index.php/NLopt_Python_Reference
def eqconstr(result,x,grad):
	if grad.size > 0:
		print "gradient to be implemented"
		
		
	nodes = ambisoniC(PHI,THETA,'basic',DEG,0)	# calculates the coefficients in the basic or naive scheme. 
							# These values are used to figure out which speakers lay in the nodes of some spherical harmonic
	for i in range((DEG+1)**2*NSPK):
		if np.asarray(abs(nodes)<10e-8).reshape(1,((DEG+1)**2*NSPK))[0,i]: result[i] = x[i] # putting to zero the speakers that are in a node of a SH
	return


def inconstr(result,x,grad):
	if grad.size > 0:
		print "gradient to be implemented"
		
		
	nodes = ambisoniC(PHI,THETA,'basic',DEG,0)	# calculates the coefficients in the basic or naive scheme. 
							# These values are used to figure out which speakers lay in the nodes of some spherical harmonic
	for i in range((DEG+1)**2*NSPK):
		if np.asarray(nodes>0.).reshape(1,((DEG+1)**2*NSPK))[0,i]: result[i] > x[i]		# keeping the sign
		if np.asarray(nodes<0.).reshape(1,((DEG+1)**2*NSPK))[0,i]: result[i] < x[i]		# keeping the sign
	return


def downscalematrix(CoeffMat,matched):
    '''
    * here we need two things:
    * - the matrix of coefficients that we want to downscale 
    * - the vector of paired speakers
    *
    * The downscaled (shrinked) matrix will be our output here
    '''
    for idx, spki in enumerate(matched):
        if (spki >= 0): 
            # the new matrix has to have the idx rows only
            if (idx==0): ResShrinked = CoeffMat.T[idx]
            else: ResShrinked = np.vstack((ResShrinked, CoeffMat.T[idx]))
    if MATCHSPK: return ResShrinked.T
    else: return CoeffMat


def upscalematrix(ResShrinked,matched):
    '''
    * here we need three things:
    * - the matrix of coefficients that we want to downscale 
    * - the vector of paired speakers
    * - the vector of SH signs that change with left/right symmetry
    *
    * The upscaled matrix will be our output here
    '''
    
    if (np.shape(ResShrinked) == ((DEG+1)**2,NSPK)) or (np.shape(ResShrinked) == ((DEG*2+1),NSPK)):
        return ResShrinked
        
    else:
        counter = 0 # you don't strictly need a counter, you can use len(shrinked) ... I guess

        for idx, spki in enumerate(matched):
            if (idx == 0) : 
                ResUp = ResShrinked.T[counter]
                counter += 1
                continue

            if (spki >= 0): 
               # if (idx == 0) : 
               #     ResUp = ResShrinked.T[counter]
               # else:
               #     # the new matrix has to have the idx rows only
               #     ResUp = np.vstack((ResUp, ResShrinked.T[counter]))
                
                ResUp = np.vstack((ResUp, ResShrinked.T[counter]))
                counter += 1
        
            elif (spki < 0):
                # search for the matched speaker
                # (which is the one labelled with positive spki)
                #wantedrow = np.argwhere(matched[matched>=0] == abs(matched[idx]))[0][0] # slow version
                wantedrow = MAP[-spki-1] # faster version
                ResUp = np.vstack((ResUp, ResShrinked.T[wantedrow]*SIGNSvec)) # multiply with vetor of signs
            else: sys.exit("Ehm ... this is impossible! [upscalematrix function]")
        
        if MATCHSPK: return ResUp.T
        else: return ResShrinked


