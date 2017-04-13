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
from argparse import ArgumentParser
import logging
import numpy as np
import pickle
import os

import ini_parser as prs
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
cfg = prs.ConfigConstants(configfile)
ambicl = fct.Ambisonic(cfg)
support = fct.Support(cfg, ambicl)

logging.info("You choose %s at order: %d " % (cfg.DEC, cfg.DEG))
logging.info("The number of speakers is %d " % cfg.NSPK)
logging.info("Number of sampling points %d " % cfg.NPOINTS)



# INITIAL PARAMETERS
# calculating ambisonics coefficients in test directions
coffs_in_test_directions = ambicl.get_ambisonic_coeffs(phi=cfg.phiTest, theta=cfg.thetaTest, DEC='basic', DEG=cfg.DEG, inversion=0)
support.coffs_in_test_directions = coffs_in_test_directions
initial_guess_projection = ambicl.get_ambisonic_coeffs(inversion=0)  # initial guess
initial_guess_pseudoinv = ambicl.get_ambisonic_coeffs(inversion=1)  # initial guess

obj_function_at_proj_guess  = support.objective_function(initial_guess_projection, coffs_in_test_directions[:])
obj_function_at_pinv_guess = support.objective_function(initial_guess_pseudoinv, coffs_in_test_directions[:])

initial_best_guess = initial_guess_projection
if obj_function_at_proj_guess > obj_function_at_pinv_guess and np.asarray([abs(initial_guess_pseudoinv) < 1.]).all():
    initial_best_guess = initial_guess_pseudoinv  # This way we choose the starting point closest to the minimum

sij = support.Sij(initial_best_guess, coffs_in_test_directions[:], cfg.NPOINTS)
pressure0, V0, energyD0, J0, Vradial0, Jradial0, Vtang0, Jtang0 = support.directional_components(sij, cfg.phiTest, cfg.thetaTest)

######################
#  Initial PLOTTING  #
######################
fct.controls(cfg.NSPK, cfg.DEG)  # checks before running the main program
# Plot speakers
ptt.threed_polar_plot(cfg.PHI, cfg.THETA, 1, True)
# Plot evaluations points grid
if cfg.exclude_with_theta_binary_mask:
    ptt.threed_polar_plot(cfg.phiTest * cfg.WbinVec, cfg.thetaTest * cfg.WbinVec, 1)
if cfg.autoexclude_regions_with_no_spkrs_smooth or cfg.autoexclude_regions_with_no_spkrs_binary:
    ptt.threed_polar_plot(cfg.phiTest * cfg.WremVec, cfg.thetaTest * cfg.WremVec, 1)

# print "\nIf you are happy with the distribution of speakers press a key \nOtherwise abort pressing ctrl+c"
# wait = raw_input()

phi_plot, theta_plot = ptt.points_over_sphere_plotting(cfg.SEED, cfg)
Npt = len(phi_plot)
phi = np.asarray(phi_plot)
theta = np.asarray(theta_plot)

coeffs_plot = ambicl.get_ambisonic_coeffs(phi_plot, theta_plot, 'basic', cfg.DEG, 0)
sij = support.Sij(initial_guess_projection, coeffs_plot, Npt)

pressureZ, VZ, energyDZ, JZ, VradialZ, JradialZ, VtangZ, JtangZ = support.directional_components(sij, phi_plot, theta_plot)


if cfg.DEC == 'basic':
    ptt.polar_plot(cfg, "Naive-Horizontal", phi[theta == 0.], pressureZ[theta == 0.], VradialZ[theta == 0.], VtangZ[theta == 0.],
                   ('pressure', 'velocity L', 'velocity T'))
else:
    ptt.polar_plot(cfg, "Naive-Horizontal", phi[theta == 0.], energyDZ[theta == 0.], JradialZ[theta == 0.], JtangZ[theta == 0.],
                   ('energy', 'intensity L', 'intensity T'))

if cfg.nD == '3D':
    if cfg.DEC == 'basic':
        ptt.polar_plot(cfg, "Naive-Vertical", theta[phi == 0], pressureZ[phi == 0], VradialZ[phi == 0], VtangZ[phi == 0],
                       ('pressure', 'velocity L', 'velocity T'))
    else:
        ptt.polar_plot(cfg, "Naive-Vertical", theta[phi == 0], energyDZ[phi == 0], JradialZ[phi == 0], JtangZ[phi == 0],
                       ('energy', 'intensity L', 'intensity T'))

##################
#  MINIMIZATION  #
##################
start = time.time()
# here we have to shrink the initGuess
initial_best_guess   = support.downscalematrix(initial_best_guess, cfg.MATCHED)
initvect    = support.mat2vec(initial_best_guess)

# Global Optimization
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
    # Local Optimization
    # LN_COBYLA, LN_BOBYQA, LN_NEWUOA, LN_NEWUOA_BOUND, LN_PRAXIS, LN_NELDERMEAD, LN_SBPLX
    if run > 0:
        ResCoeff = support.downscalematrix(ResCoeff, cfg.MATCHED)
        initvect = support.mat2vec(ResCoeff)

    opt = nlopt.opt(nlopt.LN_SBPLX, len(initvect))
    opt.set_min_objective(support.objective_function)
    tol = np.asarray([0.1] * len(initvect))
    # opt.add_equality_mconstraint(equality_constr, tol)
    # opt.add_inequality_mconstraint(inequality_constr, tol)
    opt.set_initial_step([0.2] * len(initvect))
    if run > 0: opt.set_initial_step([0.9] * len(initvect))
    # opt.set_initial_step(np.random.rand(len(initvect)));

    upbound = [1.] * len(initvect)
    lowbound = [-1.] * len(initvect)

    nodes = ambicl.get_ambisonic_coeffs(cfg.PHI, cfg.THETA, 'basic', cfg.DEG, 0)
    nodes = support.downscalematrix(nodes, cfg.MATCHED)
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

    # getting the NEW minimization matrix
    res = opt.optimize(initvect)
    result = opt.last_optimize_result()  # alternative way to get the result

    ResCoeff = support.vec2mat(res)
    ResCoeff = support.upscalematrix(ResCoeff, cfg.MATCHED)

    print "Function value for Naive: ", support.objective_function(initial_guess_projection, coffs_in_test_directions[:]), " for Pinv: ", \
        support.objective_function(initial_guess_pseudoinv, coffs_in_test_directions[:]), " for NL-min: ", \
        support.objective_function(ResCoeff, coffs_in_test_directions[:]), " Elapsed time: ", time.time() - minstart
    prog.append(support.objective_function(ResCoeff, coffs_in_test_directions[:]))
    tprogress.append(time.time() - minstart)
    minstart = time.time()

    #####################
    # exit condition
    # if (run>0 and (str(prog[run-1])[0:6]==str(prog[run])[0:6] or prog[run]>min(prog)+1)): break # for PRAXIS use this
    if (run > 0 and ("{:7.3f}".format(prog[run - 1]) == "{:7.3f}".format(prog[run]) or
                             prog[run] > min(prog) + 1) or not cfg.mute_small_coeffs):
        break  # for SBPLX use this
    run += 1

#############
#  RESULTS  #
#############
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

ResCoeff = support.vec2mat(res)
# np.reshape(res, ((DEG+1)**2,NSPKmatch))
ResCoeff = support.upscalematrix(ResCoeff, cfg.MATCHED)

print "\nLF matrix\n", initial_guess_pseudoinv.T
print "\nCoefficients \n", ResCoeff.T
if ResCoeff.T[abs(ResCoeff.T > 1.)].any(): logging.warning("WARNING: You reached a bad minimum.")

##############
#  PLOTTING  #
##############
# results plots
coeffs_plot = ambicl.get_ambisonic_coeffs(phi_plot, theta_plot, 'basic', cfg.DEG, 0)
sij = support.Sij(ResCoeff, coeffs_plot, Npt)
pressure, V, energyD, J, Vradial, Jradial, Vtang, Jtang = support.directional_components(sij, phi_plot, theta_plot)
if cfg.DEC == 'basic':
    ptt.polar_plot(cfg, "Horizontal", phi[theta == 0.], pressure[theta == 0.], Vradial[theta == 0.], Vtang[theta == 0.],
                   ('pressure', 'velocity L', 'velocity T'))
else:
    ptt.polar_plot(cfg, "Horizontal", phi[theta == 0.], energyD[theta == 0.], Jradial[theta == 0.], Jtang[theta == 0.],
                  ('energy', 'intensity L', 'intensity T'))


if cfg.nD == '3D':
    if cfg.DEC == 'basic':
        ptt.polar_plot(cfg, "Vertical", theta[phi == 0], pressure[phi == 0], Vradial[phi == 0], Vtang[phi == 0],
                       ('pressure', 'velocity L', 'velocity T'))
    else:
        ptt.polar_plot(cfg, "Vertical", theta[phi == 0], energyD[phi == 0], Jradial[phi == 0], Jtang[phi == 0],
                       ('energy', 'intensity L', 'intensity T'))

###threed_polar_plot(phi, theta, Jradial)

print "Elapsed time ", time.time() - start

###############################
# output some files
filename = "%s-%s-rem%s" % (cfg.DEG, cfg.DEC, cfg.autoexclude_regions_with_no_spkrs_binary)
np.savetxt(filename + ".txt", ResCoeff.T, fmt="%f", delimiter="  ")
np.savetxt(filename + "-Gpinv.txt", initial_guess_pseudoinv.T, fmt="%f", delimiter="  ")
np.savetxt(filename + "-G0.txt", initial_guess_projection.T, fmt="%f", delimiter="  ")

# Save state for later reuse with utilities
cfg.projection_guess_matrix = initial_guess_projection.T
cfg.pseudoinv_guess_matrix = initial_guess_pseudoinv.T
cfg.obj_minimization_matrix = ResCoeff.T
outfile = os.path.split(configfile)[-1]
outfile = outfile.split('.')[0]
outfile = "results_for_%s_%s_%s.idhoa" % (outfile, cfg.DEG, cfg.DEC)
pickle.dump(cfg, open(outfile, "wb"))


# save MATLAB output
# patched by Aaron Heller for compatibility with adt
# bitbucket.org/ambidecodertoolbox/adt.git
if cfg.is_mat_out_file:
    import scipy.io
    import os.path
    
    scipy.io.savemat(os.path.expanduser(cfg.matlab_filename),
                    {'ResCoeff': ResCoeff,
                     'Naive': initial_guess_projection,
                     'Pinv': initial_guess_pseudoinv,
                     'PHI': cfg.PHI,
                     'THETA': cfg.THETA,
                     'DEC': cfg.DEC,
                     'DEG': cfg.DEG,
                     'NSPK': cfg.NSPK,
                     'thetaThreshold': cfg.thetaThreshold,
                     'phiTest': cfg.phiTest,
                     'thetaTest': cfg.thetaTest,
                     'NPOINTS': cfg.NPOINTS,
                     'WfrontVec': cfg.WfrontVec,
                     'WplaneVec': cfg.WplaneVec,
                     'WbinVec': cfg.WbinVec,
                     'WremVec': cfg.WremVec})


print "Function value for Naive: ", support.objective_function(initial_guess_projection, coffs_in_test_directions[:]), " for Pinv: ", \
    support.objective_function(initial_guess_pseudoinv, coffs_in_test_directions[:]), " for NL-min: ", support.objective_function(ResCoeff, coffs_in_test_directions[:])

# wait = raw_input()
