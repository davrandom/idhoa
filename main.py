#!/usr/bin/python

'''
* This file is part of IDHOA software.
*
* Copyright (C) 2013 Barcelona Media - www.barcelonamedia.org
* (Written by Davide Scaini <davide.scaini@barcelonamedia.org> for Barcelona Media)
*
* The clone https://github.com/davrandom/idhoa and its modifications are under
* Copyright (C) 2014 Davide Scaini <davide.scaini@upf.edu>
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

import nlopt
# on Mac:
# export PYTHONPATH=$PYTHONPATH:/usr/local/lib/python2.7/site-packages
import time
import os
from argparse import ArgumentParser
import logging
import sys
import numpy as np

import auxiliary as aux
import functions as fct
import plotting as ptt

# parse cli arguments
parser = ArgumentParser(description="Generate Higher Order Ambisonics decoding for Irregular Layouts")
parser.add_argument('configfile', metavar='configfile', type=str,
                    help="Give a configuration file to calculate your coefficients. (Look at the example.ini)")
args = parser.parse_args()

# logging
LOGGING_LEVEL = logging.DEBUG
logging.basicConfig(
        format='%(filename)s: %(lineno)d: %(levelname)s: %(message)s',
        level=LOGGING_LEVEL)

# initialize configuration
configfile = args.configfile
conf = aux.configConst(configfile)
# because of how nlopt call to minimization function is made,
# the easiest way to share the configuration with "function" is using builtins
ambicl = fct.ambisoniClass(conf)
support = fct.Support(conf, ambicl)

logging.info("You choose %s at order: %d " % (conf.DEC, conf.DEG))
logging.info("The number of speakers is %d " % conf.NSPK)
logging.info("Number of sampling points %d " % conf.NPOINTS)

############################
## Saving stuff in a file ##
############################
pyfilename = conf.case + "-" + str(conf.DEG) + "-" + str(conf.DEC) + "-rem" + str(conf.AUTOREM) + "-sym" + str(conf.MATCHSPK)

if (conf.DEC == "basic"): pyfilename += "CP" + str(conf.CP) + "CV" + str(conf.CV) + ".py"
if (conf.DEC == "maxRe"): pyfilename += "CR" + str(conf.CR) + "CT" + str(conf.CT) + "CE" + str(conf.CE) + ".py"
if (conf.DEC == "phase"): pyfilename += "CR" + str(conf.CR) + "CT" + str(conf.CT) + "CPH" + str(conf.CPH) + ".py"

# if os.path.exists(pyfilename):
#     message = "Looks like you already have a configuration/output file with the same name. \nRemove or rename " + pyfilename + " and run again idhoa."
#     sys.exit(message)

# controls(NSPK,DEG) ## checks before running the main program
# SpeakersPlotting(conf.PHI, conf.THETA, 1)
# if conf.WBIN: ptt.SpeakersPlotting(conf.phiTest * conf.WbinVec, conf.thetaTest * conf.WbinVec, 1)
# if conf.WAUTOREM or conf.AUTOREM: ptt.SpeakersPlotting(conf.phiTest * conf.WremVec, conf.thetaTest * conf.WremVec, 1)

# print "\nIf you are happy with the distribution of speakers press a key \nOtherwise abort pressing ctrl+c"
# wait = raw_input()

start = time.time()

## INITIAL PARAMETERS
# calculating ambisonics coefficients in test directions
coeffDir = ambicl.get_ambisonic_coeffs(phi=conf.phiTest, theta=conf.thetaTest, DEC='basic', DEG=conf.DEG, inversion=0)
support.coeffDir = coeffDir
Guess0 = ambicl.get_ambisonic_coeffs(inversion=0)  # initial guess
GuessPinv = ambicl.get_ambisonic_coeffs(inversion=1)  # initial guess

f0  = support.function(Guess0, coeffDir[:])
fPV = support.function(GuessPinv, coeffDir[:])

# FIXME I don't like this, call the initial Guess like InitGuess and don't mix 0 with PV stuff they are two different things
InitGuess = Guess0
if f0 > fPV and np.asarray(
    [abs(GuessPinv) < 1.]).all(): InitGuess = GuessPinv  # This way we choose the starting point closest to the minimum

sij = support.Sij(InitGuess, coeffDir[:], conf.NPOINTS)
pressure0, V0, energyD0, J0, Vradial0, Jradial0, Vtang0, Jtang0 = support.physDir(sij, conf.phiTest, conf.thetaTest)

####################################
## Initial PLOTTING
####################################
phiPl, thetaPl = ptt.PlSpherePtRotat(conf.SEED, conf)
Npt = len(phiPl)
phi = np.asarray(phiPl)
theta = np.asarray(thetaPl)

coeffPl = ambicl.get_ambisonic_coeffs(phiPl, thetaPl, 'basic', conf.DEG, 0)
sij = support.Sij(Guess0, coeffPl, Npt)

pressureZ, VZ, energyDZ, JZ, VradialZ, JradialZ, VtangZ, JtangZ = support.physDir(sij, phiPl, thetaPl)

if conf.DEC != 'basic':
    ptt.Polar(conf, "Naive-Horizontal", phi[theta == 0.], energyDZ[theta == 0.], JradialZ[theta == 0.], JtangZ[theta == 0.],
              ('energy', 'intensity L', 'intensity T'))
if conf.DEC == 'basic':
    ptt.Polar(conf, "Naive-Horizontal", phi[theta == 0.], pressureZ[theta == 0.], VradialZ[theta == 0.], VtangZ[theta == 0.],
              ('pressure', 'velocity L', 'velocity T'))

if conf.nD == '3D':
    if conf.DEC != 'basic':
        ptt.Polar(conf, "Naive-Vertical", theta[phi == 0], energyDZ[phi == 0], JradialZ[phi == 0], JtangZ[phi == 0],
                  ('energy', 'intensity L', 'intensity T'))
    if conf.DEC == 'basic':
        ptt.Polar(conf, "Naive-Vertical", theta[phi == 0], pressureZ[phi == 0], VradialZ[phi == 0], VtangZ[phi == 0],
                  ('pressure', 'velocity L', 'velocity T'))

#####################################
## MINIMIZATION
#####################################
# here we have to shrink the initGuess
InitGuess   = support.downscalematrix(InitGuess, conf.MATCHED)
initvect    = support.mattov(InitGuess)

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
        ResCoeff = support.downscalematrix(ResCoeff, conf.MATCHED)
        initvect = support.mattov(ResCoeff)

    opt = nlopt.opt(nlopt.LN_SBPLX, len(initvect))
    opt.set_min_objective(support.function)
    tol = np.asarray([0.1] * len(initvect))
    # opt.add_equality_mconstraint(eqconstr, tol)
    # opt.add_inequality_mconstraint(inconstr, tol)
    opt.set_initial_step([0.2] * len(initvect))
    if run > 0: opt.set_initial_step([0.9] * len(initvect))
    # opt.set_initial_step(np.random.rand(len(initvect)));

    upbound = [1.] * len(initvect)
    lowbound = [-1.] * len(initvect)

    nodes = ambicl.get_ambisonic_coeffs(conf.PHI, conf.THETA, 'basic', conf.DEG, 0)
    nodes = support.downscalematrix(nodes, conf.MATCHED)
    # calculates the coefficients in the basic or naive scheme.
    # These values are used to figure out which speakers lay in the nodes of some spherical harmonic
    nodes, upbound, lowbound = support.zeroing_bounds(nodes, upbound, lowbound, run)
    if run > 0: nodes, upbound, lowbound = support.zeroing_bounds(res, upbound, lowbound, run)

    # number of non-zero elements
    logging.info( "%d non-zero elements" % len(np.nonzero(initvect)[0]) )
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
    result = opt.last_optimize_result()  # alternative way to get the result

    ResCoeff = support.vtomat(res)
    ResCoeff = support.upscalematrix(ResCoeff, conf.MATCHED)

    print "Function value for Naive: ", support.function(Guess0, coeffDir[:]), " for Pinv: ", \
        support.function(GuessPinv, coeffDir[:]), " for NL-min: ", \
        support.function(ResCoeff, coeffDir[:]), " Elapsed time: ", time.time() - minstart
    prog.append(support.function(ResCoeff, coeffDir[:]))
    tprogress.append(time.time() - minstart)
    minstart = time.time()

    #####################
    ## exit condition
    # if (run>0 and (str(prog[run-1])[0:6]==str(prog[run])[0:6] or prog[run]>min(prog)+1)): break # for PRAXIS use this
    if (run > 0 and ("{:7.3f}".format(prog[run - 1]) == "{:7.3f}".format(prog[run]) or
                             prog[run] > min(prog) + 1) or not conf.MUTESPKRS):
        break  # for SBPLX use this
    run += 1

#####################################
## RESULTS 
#####################################
print "Minimization Results:"
if result == nlopt.SUCCESS:         logging.info("Successful minimization")
if result == nlopt.STOPVAL_REACHED: logging.info("Optimization stopped because stopval was reached")
if result == nlopt.FTOL_REACHED:    logging.info("Optimization stopped because ftol_rel or ftol_abs was reached")
if result == nlopt.XTOL_REACHED:    logging.info("Optimization stopped because xtol_rel or xtol_abs was reached")
if result == nlopt.MAXEVAL_REACHED: logging.info("Optimization stopped because maxeval was reached")
if result == nlopt.MAXTIME_REACHED: logging.info("Optimization stopped because maxtime was reached")
if result == nlopt.FAILURE:         logging.error("Minimization FAILED")
if result == nlopt.INVALID_ARGS:    logging.error("Invalid arguments (e.g. lower bounds are bigger than upper bounds, "
                                                  "an unknown algorithm was specified, etcetera)")
if result == nlopt.OUT_OF_MEMORY:   logging.error("Ran out of memory.")
if result == nlopt.ROUNDOFF_LIMITED:logging.error("Halted because roundoff errors limited progress. "
                                                  "(In this case, the optimization still typically returns a useful result.)")
if result == nlopt.FORCED_STOP:     logging.error("Halted because of a forced termination: "
                                                  "the user called nlopt_force_stop(opt) on the optimization's "
                                                  "nlopt_opt object opt from the user's objective function or constraints.")

ResCoeff = support.vtomat(res)
# np.reshape(res, ((DEG+1)**2,NSPKmatch))
ResCoeff = support.upscalematrix(ResCoeff, conf.MATCHED)

print "\nLF matrix\n", GuessPinv.T
print "\nCoefficients \n", ResCoeff.T
if ResCoeff.T[abs(ResCoeff.T > 1.)].any(): logging.warning("WARNING: You reached a bad minimum.")

##############################
## PLOTTING
##############################
# results plots
coeffPl = ambicl.get_ambisonic_coeffs(phiPl, thetaPl, 'basic', conf.DEG, 0)
sij = support.Sij(ResCoeff, coeffPl, Npt)
pressure, V, energyD, J, Vradial, Jradial, Vtang, Jtang = support.physDir(sij, phiPl, thetaPl)
if conf.DEC != 'basic':
   ptt.Polar(conf, "Horizontal", phi[theta == 0.], energyD[theta == 0.], Jradial[theta == 0.], Jtang[theta == 0.],
             ('energy', 'intensity L', 'intensity T'))
if conf.DEC == 'basic':
    ptt.Polar(conf, "Horizontal", phi[theta == 0.], pressure[theta == 0.], Vradial[theta == 0.], Vtang[theta == 0.],
              ('pressure', 'velocity L', 'velocity T'))

if conf.nD == '3D':
    if conf.DEC != 'basic':
        ptt.Polar(conf, "Vertical", theta[phi == 0], energyD[phi == 0], Jradial[phi == 0], Jtang[phi == 0],
                  ('energy', 'intensity L', 'intensity T'))
    if conf.DEC == 'basic':
        ptt.Polar(conf, "Vertical", theta[phi == 0], pressure[phi == 0], Vradial[phi == 0], Vtang[phi == 0],
                  ('pressure', 'velocity L', 'velocity T'))

###SpeakersPlotting(phi, theta, Jradial)

print "Elapsed time ", time.time() - start

###############################
# output some files
filename = str(conf.DEG) + "-" + str(conf.DEC) + "-rem" + str(conf.AUTOREM)
np.savetxt(filename + ".idhoa", ResCoeff.T, fmt="%f", delimiter="  ")
np.savetxt(filename + "-Gpinv.idhoa", GuessPinv.T, fmt="%f", delimiter="  ")
np.savetxt(filename + "-G0.idhoa", Guess0.T, fmt="%f", delimiter="  ")

##
# Output to a python file for later use (as simple as import the python file)
aux.outp("Layout",  conf.case, pyfilename)
aux.outp("DEC",     conf.DEC, pyfilename)
aux.outp('PHI',     conf.PHI, pyfilename)
aux.outp('THETA',   conf.THETA, pyfilename)
aux.outp('DEC',     conf.DEC, pyfilename)
aux.outp('DEG',     conf.DEG, pyfilename)
aux.outp('nD',      conf.nD, pyfilename)
aux.outp('NSPK',    conf.NSPK, pyfilename)
aux.outp('thetaThreshold', conf.thetaThreshold, pyfilename)
aux.outp('phiTest', conf.phiTest, pyfilename)
aux.outp('thetaTest', conf.thetaTest, pyfilename)
aux.outp('NPOINTS', conf.NPOINTS, pyfilename)
aux.outp('SEED',    conf.SEED, pyfilename)
aux.outp('WfrontVec', conf.WfrontVec, pyfilename)
aux.outp('WplaneVec', conf.WplaneVec, pyfilename)
aux.outp('WbinVec', conf.WbinVec, pyfilename)
aux.outp('WremVec', conf.WremVec, pyfilename)
aux.outp("CP ", conf.CP, pyfilename)
aux.outp("CV ", conf.CV, pyfilename)
aux.outp("CE ", conf.CE, pyfilename)
aux.outp("CR ", conf.CR, pyfilename)
aux.outp("CT ", conf.CT, pyfilename)
aux.outp("CPH", conf.CPH, pyfilename)
aux.outp("MATCHSPK", conf.MATCHSPK, pyfilename)
aux.outp("MATCHTOL", conf.MATCHTOL, pyfilename)
aux.outp("ResCoeff", ResCoeff, pyfilename)
aux.outp("Guess0",   Guess0, pyfilename)
aux.outp("GuessPinv", GuessPinv, pyfilename)
aux.outp("fprogress", prog, pyfilename)
aux.outp("nonzeroprog", nonzeroprog, pyfilename)
aux.outp("tprogress", tprogress, pyfilename)
aux.outp("time", time.time() - start, pyfilename)

##
# save MATLAB output
# patched by AH for compatibility with adt 
# bitbucket.org/ambidecodertoolbox/adt.git
#import scipy.io
#import os.path
#
#scipy.io.savemat(os.path.expanduser(OUTFILE),
#                 {'ResCoeff': ResCoeff,
#                  'Naive': Guess0,
#                  'Pinv': GuessPinv,
#                  'PHI': PHI,
#                  'THETA': THETA,
#                  'DEC': DEC,
#                  'DEG': DEG,
#                  'NSPK': NSPK,
#                  'thetaThreshold': thetaThreshold,
#                  'phiTest': phiTest,
#                  'thetaTest': thetaTest,
#                  'NPOINTS': NPOINTS,
#                  'WfrontVec': WfrontVec,
#                  'WplaneVec': WplaneVec,
#                  'WbinVec': WbinVec,
#                  'WremVec': WremVec})
##
###############################


print "Function value for Naive: ", support.function(Guess0, coeffDir[:]), " for Pinv: ", \
    support.function(GuessPinv, coeffDir[:]), " for NL-min: ", support.function(ResCoeff, coeffDir[:])

wait = raw_input()
