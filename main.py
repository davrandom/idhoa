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



import numpy as np
import math as mh
import sys
import nlopt
# on Mac:
# export PYTHONPATH=$PYTHONPATH:/usr/local/lib/python2.7/site-packages
import time
import os
from argparse import ArgumentParser

import auxiliary as aux

parser = ArgumentParser(description="Generate Higher Order Ambisonics decoding for Irregular Layouts")
parser.add_argument('configfile',metavar='configfile',type=str,help="Give a configuration file to calculate your coefficients. (Look at the example.ini)")
args = parser.parse_args()

configfile = args.configfile
aux.configConst(configfile)

from functions import *
from plotting import *


print "You choose ", DEC, " at order: ",DEG
print "The number of speakers is ", NSPK
print "Number of sampling points ", NPOINTS

############################
## Saving stuff in a file ##
############################
pyfilename = case+"-"+str(DEG)+"-"+str(DEC)+"-rem"+str(AUTOREM)+"-sym"+str(MATCHSPK)

if (DEC=="basic"): pyfilename += "CP"+str(CP)+"CV"+str(CV)+".py" 
if (DEC=="maxRe"): pyfilename += "CR"+str(CR)+"CT"+str(CT)+"CE"+str(CE)+".py"
if (DEC=="phase"): pyfilename += "CR"+str(CR)+"CT"+str(CT)+"CPH"+str(CPH)+".py"

if os.path.exists(pyfilename): 
    message = "Looks like you already have a configuration/output file with the same name. \nRemove or rename "+pyfilename+" and run again idhoa."
    sys.exit(message)
 

# controls(NSPK,DEG) ## checks before running the main program
SpeakersPlotting(PHI,THETA,1)
if WBIN: SpeakersPlotting(phiTest*WbinVec,thetaTest*WbinVec,1)
if AUTOREM: SpeakersPlotting(phiTest*WremVec,thetaTest*WremVec,1)

start = time.time()


## INITIAL PARAMETERS 
coeffDir = ambisoniC(phiTest,thetaTest,'basic',DEG,0) # calculating ambisonics coefficients in test directions
Guess0 = ambisoniC(PHI,THETA,DEC,DEG,0)     # initial guess
GuessPinv = ambisoniC(PHI,THETA,DEC,DEG,1)     # initial guess


f0 = function(Guess0,coeffDir[:])
fPV = function(GuessPinv,coeffDir[:])

#FIXME I don't like this, call the initial Guess like InitGuess and don't mix 0 with PV stuff they are two different things
InitGuess = Guess0
if f0 > fPV and np.asarray([abs(GuessPinv)<1.]).all(): InitGuess = GuessPinv     # This way we choose the starting point closest to the minimum


sij = Sij(InitGuess,coeffDir[:],NPOINTS)
pressure0, V0, energyD0, J0, Vradial0, Jradial0, Vtang0, Jtang0 = physDir(sij,phiTest,thetaTest)

####################################
## Initial PLOTTING
####################################
phiPl,thetaPl = PlSpherePtRotat(SEED)
Npt = len(phiPl)
phi = np.asarray(phiPl)
theta = np.asarray(thetaPl)

coeffPl = ambisoniC(phiPl,thetaPl,'basic',DEG,0)
sij = Sij(Guess0,coeffPl,Npt)
#sij = Sij(GuessPinv,coeffPl,Npt) # if you want to look at the pinv

pressureZ, VZ, energyDZ, JZ, VradialZ, JradialZ, VtangZ, JtangZ = physDir(sij,phiPl,thetaPl)

if DEC!='basic':
    Polar("Naive-Horizontal",phi[theta==0.],energyDZ[theta==0.],JradialZ[theta==0.],JtangZ[theta==0.],('energy','intensity L','intensity T'))
if DEC=='basic':
    Polar("Naive-Horizontal",phi[theta==0.],pressureZ[theta==0.],VradialZ[theta==0.],VtangZ[theta==0.],('pressure','velocity L','velocity T'))


if nD == '3D':
    if DEC!='basic':
        Polar("Naive-Vertical",theta[phi==0],energyDZ[phi==0],JradialZ[phi==0],JtangZ[phi==0],('energy','intensity L','intensity T'))
    if DEC=='basic':
        Polar("Naive-Vertical",theta[phi==0],pressureZ[phi==0],VradialZ[phi==0],VtangZ[phi==0],('pressure','velocity L','velocity T'))




#####################################
## MINIMIZATION
#####################################

## here whe have to shrink the initGuess

InitGuess = downscalematrix(InitGuess,MATCHED)
initvect = mattov(InitGuess)
## Global Optimization
# GN_DIRECT, GN_DIRECT_L, GN_ORIG_DIRECT, GN_ORIG_DIRECT_L, GN_DIRECT_NOSCAL, GN_DIRECT_L_NOSCAL, GN_DIRECT_L_RAND_NOSCAL
# GN_CRS2_LM
# G_MLSL_LDS, G_MLSL
# GD_STOGO, GD_STOGO_RAND
# GN_ISRES
prog = []
tprogress = []
nonzeroprog = []
run = 0
minstart = time.time()

while True:
    ## Local Optimization
    # LN_COBYLA, LN_BOBYQA, LN_NEWUOA, LN_NEWUOA_BOUND, LN_PRAXIS, LN_NELDERMEAD, LN_SBPLX
    if run > 0: 
        ResCoeff = downscalematrix(ResCoeff,MATCHED)
        initvect = mattov(ResCoeff)
    
    opt = nlopt.opt(nlopt.LN_SBPLX,len(initvect)) 
    opt.set_min_objective(function)
    tol = np.asarray([0.1]*len(initvect))
    #opt.add_equality_mconstraint(eqconstr, tol)    
    #opt.add_inequality_mconstraint(inconstr, tol)  
    opt.set_initial_step([0.2]*len(initvect))
    if run > 0: opt.set_initial_step([0.9]*len(initvect))
    #opt.set_initial_step(np.random.rand(len(initvect))); 


    upbound = [1.]*len(initvect)
    lowbound = [-1.]*len(initvect)
    
    nodes = ambisoniC(PHI,THETA,'basic',DEG,0)
    nodes = downscalematrix(nodes,MATCHED)
    # calculates the coefficients in the basic or naive scheme. 
    # These values are used to figure out which speakers lay in the nodes of some spherical harmonic
    nodes, upbound, lowbound = zeroing_bounds(nodes, upbound, lowbound, run)
    if run>0: nodes, upbound, lowbound = zeroing_bounds(res, upbound, lowbound, run)


    # number of non-zero elements 
    print len(np.nonzero(initvect)[0]), " non-zero elements"
    nonzeroprog.append(len(np.nonzero(initvect)[0]))

    opt.set_upper_bounds(upbound)
    opt.set_lower_bounds(lowbound)

    opt.set_xtol_abs(10e-6)
    opt.set_ftol_abs(10e-6)

    if run > 0: 
        opt.set_xtol_abs(10e-8)
        opt.set_ftol_abs(10e-8)

    ## getting the NEW minimization matrix
    res = opt.optimize(initvect)
    result = opt.last_optimize_result() #alternative way to get the result
    
    ResCoeff = vtomat(res)
    ResCoeff = upscalematrix(ResCoeff,MATCHED)

    print "Function value for Naive: ", function(Guess0,coeffDir[:]), " for Pinv: ", function(GuessPinv,coeffDir[:]), " for NL-min: ", function(ResCoeff,coeffDir[:]), " Elapsed time: ", time.time()-minstart
    prog.append(function(ResCoeff,coeffDir[:]))
    tprogress.append(time.time()-minstart)
    minstart = time.time()
    
    #####################
    ## exit condition
    # if (run>0 and (str(prog[run-1])[0:6]==str(prog[run])[0:6] or prog[run]>min(prog)+1)): break # for PRAXIS use this
    if (run>0 and ("{:7.4f}".format(prog[run-1])=="{:7.4f}".format(prog[run]) or prog[run]>min(prog)+1) or not MUTESPKRS ): break # for SBPLX use this
    run+=1




#####################################
## RESULTS 
#####################################
print "Minimization Results:"
if result == nlopt.SUCCESS: print "Successful minimization"
if result == nlopt.STOPVAL_REACHED: print "Optimization stopped because stopval was reached"
if result == nlopt.FTOL_REACHED: print "Optimization stopped because ftol_rel or ftol_abs was reached"
if result == nlopt.XTOL_REACHED: print "Optimization stopped because xtol_rel or xtol_abs was reached"
if result == nlopt.MAXEVAL_REACHED: print "Optimization stopped because maxeval was reached"
if result == nlopt.MAXTIME_REACHED: print "Optimization stopped because maxtime was reached"
if result == nlopt.FAILURE: print "Minimization FAILED" 
if result == nlopt.INVALID_ARGS: print "Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera)"
if result == nlopt.OUT_OF_MEMORY: print "Ran out of memory."
if result == nlopt.ROUNDOFF_LIMITED: print "Halted because roundoff errors limited progress. (In this case, the optimization still typically returns a useful result.)" 
if result == nlopt.FORCED_STOP: print "Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization's nlopt_opt object opt from the user's objective function or constraints." 


ResCoeff = vtomat(res)
#np.reshape(res, ((DEG+1)**2,NSPKmatch))
ResCoeff = upscalematrix(ResCoeff,MATCHED)

print "\nCoefficients \n", ResCoeff.T
if ResCoeff.T[abs(ResCoeff.T>1.)].any(): print "WARNING: You reached a bad minimum."




##############################
## PLOTTING
##############################
## results plots
coeffPl = ambisoniC(phiPl,thetaPl,'basic',DEG,0)
sij = Sij(ResCoeff,coeffPl,Npt)
pressure, V, energyD, J, Vradial, Jradial, Vtang, Jtang = physDir(sij,phiPl,thetaPl)
if DEC!='basic':
    Polar("Horizontal",phi[theta==0.],energyD[theta==0.],Jradial[theta==0.],Jtang[theta==0.],('energy','intensity L','intensity T'))
if DEC=='basic':
    Polar("Horizontal",phi[theta==0.],pressure[theta==0.],Vradial[theta==0.],Vtang[theta==0.],('pressure','velocity L','velocity T'))


if nD == '3D':
    if DEC!='basic':
        Polar("Vertical",theta[phi==0],energyD[phi==0],Jradial[phi==0],Jtang[phi==0],('energy','intensity L','intensity T'))
    if DEC=='basic':
        Polar("Vertical",theta[phi==0],pressure[phi==0],Vradial[phi==0],Vtang[phi==0],('pressure','velocity L','velocity T'))

SpeakersPlotting(phi,theta,Jradial)

print "Elapsed time ", time.time()-start

###############################
# output some files
filename = str(DEG)+"-"+str(DEC)+"-rem"+str(AUTOREM)
np.savetxt(filename+".idhoa",ResCoeff.T,fmt="%f",delimiter="  ")
np.savetxt(filename+"-Gpinv.idhoa",GuessPinv.T,fmt="%f",delimiter="  ")
np.savetxt(filename+"-G0.idhoa",Guess0.T,fmt="%f",delimiter="  ")

##
# Output to a python file for later use (as simple as import the python file)
aux.outp("Layout",case,pyfilename)
aux.outp("DEC",DEC,pyfilename)
aux.outp('PHI',PHI,pyfilename) 
aux.outp('THETA',THETA,pyfilename) 
aux.outp('DEC',DEC,pyfilename)
aux.outp('DEG',DEG,pyfilename)
aux.outp('nD',nD,pyfilename)
aux.outp('NSPK',NSPK,pyfilename)
aux.outp('thetaThreshold',thetaThreshold,pyfilename)
aux.outp('phiTest',phiTest,pyfilename) 
aux.outp('thetaTest',thetaTest,pyfilename) 
aux.outp('NPOINTS',NPOINTS,pyfilename) 
aux.outp('SEED',SEED,pyfilename) 
aux.outp('WfrontVec',WfrontVec,pyfilename) 
aux.outp('WplaneVec',WplaneVec,pyfilename)
aux.outp('WbinVec',WbinVec,pyfilename) 
aux.outp('WremVec',WremVec,pyfilename)
aux.outp("CP ", CP ,pyfilename)  
aux.outp("CV ", CV ,pyfilename)
aux.outp("CE ", CE ,pyfilename)
aux.outp("CR ", CR ,pyfilename)
aux.outp("CT ", CT ,pyfilename)
aux.outp("CPH", CPH,pyfilename)
aux.outp("MATCHSPK", MATCHSPK,pyfilename)
aux.outp("MATCHTOL", MATCHTOL,pyfilename)
aux.outp("ResCoeff",ResCoeff,pyfilename)
aux.outp("Guess0",Guess0,pyfilename)
aux.outp("GuessPinv",GuessPinv,pyfilename)
aux.outp("fprogress",prog,pyfilename)
aux.outp("nonzeroprog",nonzeroprog,pyfilename)
aux.outp("tprogress",tprogress,pyfilename)
aux.outp("time",time.time()-start,pyfilename)

##
# save MATLAB output
# patched by AH for compatibility with adt 
# bitbucket.org/ambidecodertoolbox/adt.git
import scipy.io
import os.path

scipy.io.savemat(os.path.expanduser(OUTFILE),
        {'ResCoeff':ResCoeff, 
        'PHI':PHI, 
        'THETA':THETA, 
        'DEC':DEC,
        'DEG':DEG,
        'NSPK':NSPK,
        'thetaThreshold':thetaThreshold,
        'phiTest':phiTest, 
        'thetaTest':thetaTest, 
        'NPOINTS':NPOINTS, 
        'WfrontVec':WfrontVec, 
        'WplaneVec':WplaneVec,
        'WbinVec':WbinVec, 
        'WremVec':WremVec})
##
###############################


print "Function value for Naive: ", function(Guess0,coeffDir[:]), " for Pinv: ", function(GuessPinv,coeffDir[:]), " for NL-min: ",function(ResCoeff,coeffDir[:])

wait = raw_input()

