#!/usr/bin/python

"""
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
"""

import numpy as np
from scipy import special, linalg
import math as mh
import logging
import algopy


# checks
def controls(NSPK, DEG):
    if NSPK < (DEG + 1) ** 2:
        raise ValueError("\n##!!!!!!!!!!!!!!!!!!!!!!!!!!##\nYou don't have enough speakers for this decoding order.")


def gauss_legendre(n):
    # http://www.scientificpython.net/1/post/2012/04/gausslegendre1.html
    k = np.arange(1.0, n)
    a_band = np.zeros((2, n))
    a_band[1, 0:n - 1] = k / np.sqrt(4 * k * k - 1)
    x, V = linalg.eig_banded(a_band, lower=True)
    w = 2 * np.real(np.power(V[0, :], 2))
    return x, w


def maxReG1(order):
    x, w = gauss_legendre(order + 1)
    return max(x)


def maxReCoeffs(order):
    coeffs = [special.legendre(i)(maxReG1(order)) for i in range(order + 1)]
    return coeffs


class ambisoniClass:
    def __init__(self, configClass):
        self.cfg = configClass
        self.phi = self.cfg.PHI
        self.theta = self.cfg.THETA
        self.DEC = self.cfg.DEC
        self.DEG = self.cfg.DEG
        self.nD = self.cfg.nD
        self.inversion = 0
        self.coeffs = None
        self.W = None
        self.Y = None; self.Z = None; self.X = None
        self.V = None; self.T = None; self.R = None; self.S = None; self.U = None; self.Q = None
        self.O = None; self.M = None; self.K = None; self.L = None; self.N = None; self.P = None
        self.c12 = None; self.c13 = None; self.c14 = None; self.c15 = None
        self.c16 = None; self.c17 = None; self.c18 = None; self.c19 = None; self.c20 = None
        self.c21 = None; self.c22 = None; self.c23 = None; self.c24 = None; self.c25 = None
        self.c26 = None; self.c27 = None; self.c28 = None; self.c29 = None; self.c30 = None
        self.c31 = None; self.c32 = None; self.c33 = None; self.c34 = None; self.c35 = None

    def get_ambisonic_coeffs(self, phi=None, theta=None, DEC=None, DEG=None, inversion=0):
        if isinstance(phi, np.ndarray) or isinstance(phi, list):
            self.phi = phi
        else:
            self.phi = self.cfg.PHI
        if isinstance(theta, np.ndarray) or isinstance(theta, list):
            self.theta = theta
        else:
            self.theta = self.cfg.THETA
        if DEC:
            self.DEC = DEC
        else:
            self.DEC = self.cfg.DEC
        if DEG:
            self.DEG = DEG
        else:
            self.DEG = self.cfg.DEG
        self.inversion = inversion

        ################################
        ## horizontal Ambisonics only ##
        ################################
        if self.nD == '2D':
            self._2D()
        else:
            self._3D()
        # at the end, please return something ;)
        return self.coeff

    def _2D(self):
            # things to be initializated
            NUM = len(self.phi)
            g0 = [1.] * (self.DEG + 1)
            g1 = [1.] * (self.DEG + 1)
            coeffs = [None] * (len(self.phi))
            self.W = [None] * (len(self.phi))
            if self.DEG >= 1:
                self.X = [None] * (len(self.phi))
                self.Y = [None] * (len(self.phi))
            if self.DEG >= 2:
                self.V = [None] * (len(self.phi))
                self.U = [None] * (len(self.phi))
            if self.DEG >= 3:
                self.Q = [None] * (len(self.phi))
                self.P = [None] * (len(self.phi))
            if self.DEG >= 4:
                self.c8 = [None] * (len(self.phi))
                self.c9 = [None] * (len(self.phi))
            if self.DEG >= 5:
                self.c10 = [None] * (len(self.phi))
                self.c11 = [None] * (len(self.phi))
            if self.DEG >= 6:
                raise ValueError("DEG =", self.DEG, " is not implemented yet\n")

                #####################################################
                #  Calculating the decoding dependent coefficients  #
                #####################################################
            # TOFIX TOBEDONE FIXME
            if self.DEC == 'basic':
                for i in range(0, self.DEG + 1):
                    g0[i] = 1.
                    g1[i] = 1.

            elif self.DEC == 'maxRe':
                g1p = [np.cos(i * np.pi / (2 * self.DEG + 2)) for i in range(0, self.DEG + 1)]
                g0 = [np.sqrt(float(NUM) / (self.DEG + 1))] * 2
                g1 = [g1p[i] * g0[0] for i in range(len(g1p))]

            elif self.DEC == 'phase':
                g1p = [(mh.factorial(self.DEG) ** 2 / float(mh.factorial(self.DEG + i) * mh.factorial(self.DEG - i)))
                       for i in
                       range(0, self.DEG + 1)]  # okf
                Eg1p = g1p[0] ** 2 + 2.0 * sum([g1p[i] ** 2 for i in range(1, self.DEG + 1)])
                g0[0] = np.sqrt(float(NUM) / Eg1p)
                g1 = [x * g0[0] for x in g1p]

            else:
                raise ValueError('Decoding scheme unknow: ', self.DEC, ' Possible ones are: basic, maxRe, phase')

            ##########################################
            ## Calculating the "naive" coefficients ##
            ## in the N2D convention! careful!      ##
            ##########################################
            for i in range(0, len(self.phi)):
                if self.DEG >= 0:
                    self.W[i] = 1.0 / NUM
                    # W has this sqrt(2) factor to take into account that W is not directly the pressure

                if self.DEG >= 1:
                    # 1st order
                    self.X[i] = np.sqrt(2.) * np.cos(self.phi[i]) / NUM
                    self.Y[i] = np.sqrt(2.) * np.sin(self.phi[i]) / NUM

                if self.DEG >= 2:
                    # 2nd order
                    self.U[i] = np.sqrt(2.) * np.cos(2. * self.phi[i]) / NUM
                    self.V[i] = np.sqrt(2.) * np.sin(2. * self.phi[i]) / NUM

                if self.DEG >= 3:
                    # 3rd order
                    self.P[i] = np.sqrt(2.) * np.cos(3. * self.phi[i]) / NUM
                    self.Q[i] = np.sqrt(2.) * np.sin(3. * self.phi[i]) / NUM

                if self.DEG >= 4:
                    # 4th order
                    self.c8[i] = np.sqrt(2.) * np.cos(4. * self.phi[i]) / NUM
                    self.c9[i] = np.sqrt(2.) * np.sin(4. * self.phi[i]) / NUM

                if self.DEG >= 5:
                    # 5th order
                    self.c10[i] = np.sqrt(2.) * np.cos(5. * self.phi[i]) / NUM
                    self.c11[i] = np.sqrt(2.) * np.sin(5. * self.phi[i]) / NUM

                if self.DEG >= 6:
                    # 6th order
                    self.c12[i] = np.sqrt(2.) * np.cos(6. * self.phi[i]) / NUM
                    self.c13[i] = np.sqrt(2.) * np.sin(6. * self.phi[i]) / NUM

                if self.DEG >= 7:
                    # 7th order
                    self.c14[i] = np.sqrt(2.) * np.cos(7. * self.phi[i]) / NUM
                    self.c15[i] = np.sqrt(2.) * np.sin(7. * self.phi[i]) / NUM

                if self.DEG >= 8:
                    # 8th order
                    self.c16[i] = np.sqrt(2.) * np.cos(8. * self.phi[i]) / NUM
                    self.c17[i] = np.sqrt(2.) * np.sin(8. * self.phi[i]) / NUM

                if self.DEG >= 9:
                    # 9th order
                    self.c18[i] = np.sqrt(2.) * np.cos(9. * self.phi[i]) / NUM
                    self.c19[i] = np.sqrt(2.) * np.sin(9. * self.phi[i]) / NUM

                if self.DEG > 9:
                    print "DEG =", self.DEG, " is not implemented yet\n"

            if self.DEG == 1:
                coeffs = np.array([self.W, self.Y, self.X])
            elif self.DEG == 2:
                coeffs = np.array([self.W, self.Y, self.X, self.V, self.U])
            elif self.DEG == 3:
                coeffs = np.array([self.W, self.Y, self.X, self.V, self.U, self.Q, self.P])
            elif self.DEG == 4:
                coeffs = np.array([self.W, self.Y, self.X, self.V, self.U, self.Q, self.P, self.c9, self.c8])
            elif self.DEG == 5:
                coeffs = np.array([self.W, self.Y, self.X, self.V, self.U, self.Q, self.P, self.c9, self.c8,
                                   self.c11, self.c10])
            elif self.DEG == 6:
                coeffs = np.array([self.W, self.Y, self.X, self.V, self.U, self.Q, self.P, self.c9, self.c8,
                                   self.c11, self.c10, self.c13, self.c12])
            elif self.DEG == 7:
                coeffs = np.array([self.W, self.Y, self.X, self.V, self.U, self.Q, self.P, self.c9, self.c8,
                                   self.c11, self.c10, self.c13, self.c12, self.c15, self.c14])
            elif self.DEG == 8:
                coeffs = np.array([self.W, self.Y, self.X, self.V, self.U, self.Q, self.P, self.c9, self.c8,
                                   self.c11, self.c10, self.c13, self.c12, self.c15, self.c14, self.c17, self.c16])
            elif self.DEG == 9:
                coeffs = np.array([self.W, self.Y, self.X, self.V, self.U, self.Q, self.P, self.c9, self.c8,
                                   self.c11, self.c10, self.c13, self.c12, self.c15, self.c14, self.c17, self.c16,
                                   self.c19, self.c18])

            if self.inversion == 1:
                coeffs = np.linalg.pinv(coeffs, rcond=1e-8).T / NUM
                coeffs[abs(coeffs) < 1e-8] = 0.  # because inversion gives some very small values somewhere

            ## MULTIPLYING FOR THE SELECTED DECODING SCHEME
            coeff = np.empty(coeffs.shape, dtype=np.float64)
            g1[0] = g0[0]
            for i in range(self.DEG + 1):
                for jj in range(i * 2 - 1, i * 2 + 1):
                    if jj >= 0:
                        coeff[jj] = g1[i] * coeffs[jj]

            self.coeff = coeff

    def _3D(self):
        # FIXME: this is for nD == '3D' 
        # things to be initializated
        NUM = len(self.phi)
        g1 = [None] * (self.DEG + 1)
        coeffs = [None] * (len(self.phi))
        self.W = [None] * (len(self.phi))
        if self.DEG >= 1:
            self.X = [None] * (len(self.phi))
            self.Y = [None] * (len(self.phi))
            self.Z = [None] * (len(self.phi))
        if self.DEG >= 2:
            self.V = [None] * (len(self.phi))
            self.T = [None] * (len(self.phi))
            self.R = [None] * (len(self.phi))
            self.S = [None] * (len(self.phi))
            self.U = [None] * (len(self.phi))
        if self.DEG >= 3:
            self.Q = [None] * (len(self.phi))
            self.O = [None] * (len(self.phi))
            self.M = [None] * (len(self.phi))
            self.K = [None] * (len(self.phi))
            self.L = [None] * (len(self.phi))
            self.N = [None] * (len(self.phi))
            self.P = [None] * (len(self.phi))
        if self.DEG >= 4:
            self.c16 = [None] * (len(self.phi))
            self.c17 = [None] * (len(self.phi))
            self.c18 = [None] * (len(self.phi))
            self.c19 = [None] * (len(self.phi))
            self.c20 = [None] * (len(self.phi))
            self.c21 = [None] * (len(self.phi))
            self.c22 = [None] * (len(self.phi))
            self.c23 = [None] * (len(self.phi))
            self.c24 = [None] * (len(self.phi))
        if self.DEG >= 5:
            self.c25 = [None] * (len(self.phi))
            self.c26 = [None] * (len(self.phi))
            self.c27 = [None] * (len(self.phi))
            self.c28 = [None] * (len(self.phi))
            self.c29 = [None] * (len(self.phi))
            self.c30 = [None] * (len(self.phi))
            self.c31 = [None] * (len(self.phi))
            self.c32 = [None] * (len(self.phi))
            self.c33 = [None] * (len(self.phi))
            self.c34 = [None] * (len(self.phi))
            self.c35 = [None] * (len(self.phi))
        if self.DEG >= 6:
            raise ValueError("DEG =", self.DEG, " is not implemented yet\n")

        #####################################################
        #  Calculating the decoding dependent coefficients  #
        #####################################################
        if self.DEC == 'basic':
            g1 = self._basic()

        elif self.DEC == 'maxRe':
            g1 = self._maxRe(self.phi)

        elif self.DEC == 'phase':
            g1 = self._phase(self.phi)

        elif self.DEC == 'mixed':
            g1b = self._basic()
            g1m = self._maxRe(self.phi)
            g1p = self._phase(self.phi)
            if self.DEC > 2:
                # for 3rd order 0th and 1st orders are basic, 2nd order should be maxRe, and 3rd should be inphase
                g1 = np.concatenate([g1b[0:4], g1m[4:self.DEG ** 2], g1p[self.DEG ** 2:(self.DEG + 1) ** 2]],
                                    axis=0)
            if self.DEC > 3:
                # for 3rd order 0th and 1st orders are basic, 2nd and 3rd ... unitl maxorder-1 orders should be maxRe,
                # and last order should be inphase
                g1 = np.concatenate([g1b[0:4], g1m[4:4 ** 2], g1p[4 ** 2:(self.DEG + 1) ** 2]], axis=0)

        else:
            raise ValueError('Decoding scheme unknow: ', self.DEC, ' Possible ones are: basic, maxRe, phase')

        ##########################################
        #  Calculating the "naive" coefficients  #
        ##########################################
        for i in range(0, len(self.phi)):
            if self.DEG >= 0:
                self.W[i] = 1.0 / NUM
                # W has this sqrt(2) factor to take into account that W is not directly the pressure

            if self.DEG >= 1:
                # 1st order
                self.X[i] = np.sqrt(3) * np.cos(self.phi[i]) * np.cos(
                    self.theta[i]) / NUM  # factor sqrt(3) comes from sqrt(2m+1) - see tab.3.3 p.156 Jerome Daniel
                self.Y[i] = np.sqrt(3) * np.sin(self.phi[i]) * np.cos(
                    self.theta[
                        i]) / NUM  # The decoding matrix is given for SN3D (aka 'semi-normalized') so to go to N3D
                self.Z[i] = np.sqrt(3) * np.sin(self.theta[i]) / NUM  # you have to use the alphas.

            if self.DEG >= 2:
                # 2nd order
                self.V[i] = np.sqrt(5.) * np.sqrt(3.) / 2. * np.sin(2. * self.phi[i]) * np.cos(self.theta[i]) ** 2 / NUM
                self.T[i] = np.sqrt(5.) * np.sqrt(3.) / 2. * np.sin(self.phi[i]) * np.sin(2. * self.theta[i]) / NUM
                self.R[i] = np.sqrt(5.) * (3. * np.sin(self.theta[i]) ** 2. - 1.) / 2. / NUM
                self.S[i] = np.sqrt(5.) * np.sqrt(3.) / 2. * np.cos(self.phi[i]) * np.sin(2. * self.theta[i]) / NUM
                self.U[i] = np.sqrt(5.) * np.sqrt(3.) / 2. * np.cos(2. * self.phi[i]) * np.cos(self.theta[i]) ** 2 / NUM

            if self.DEG >= 3:
                # 3rd order
                self.Q[i] = np.sqrt(7.) * np.sqrt(5. / 8.) * np.sin(3. * self.phi[i]) * np.cos(self.theta[i]) ** 3 / NUM
                self.O[i] = np.sqrt(7.) * np.sqrt(15.) / 2. * np.sin(2. * self.phi[i]) * np.sin(self.theta[i]) * np.cos(
                    self.theta[i]) ** 2 / NUM
                self.M[i] = np.sqrt(7.) * np.sqrt(3. / 8.) * np.sin(self.phi[i]) * np.cos(self.theta[i]) * (
                    5. * np.sin(self.theta[i]) ** 2 - 1.) / NUM
                self.K[i] = np.sqrt(7.) * np.sin(self.theta[i]) * (5. * np.sin(self.theta[i]) ** 2 - 3.) / 2. / NUM
                self.L[i] = np.sqrt(7.) * np.sqrt(3. / 8.) * np.cos(self.phi[i]) * np.cos(self.theta[i]) * (
                    5. * np.sin(self.theta[i]) ** 2 - 1.) / NUM
                self.N[i] = np.sqrt(7.) * np.sqrt(15.) / 2. * np.cos(2. * self.phi[i]) * np.sin(self.theta[i]) * np.cos(
                    self.theta[i]) ** 2 / NUM
                self.P[i] = np.sqrt(7.) * np.sqrt(5. / 8.) * np.cos(3. * self.phi[i]) * np.cos(self.theta[i]) ** 3 / NUM

            if self.DEG >= 4:
                # 4th order
                self.c16[i] = np.sqrt(9.) * np.sqrt(35. / 2.) * 3. / 8. * np.sin(4. * self.phi[i]) * np.cos(
                    self.theta[i]) ** 4 / NUM  # (4,-4)
                self.c17[i] = np.sqrt(9.) * np.sqrt(35.) * 3. / 4. * np.sin(
                    3. * self.phi[i]) * np.sin(self.theta[i]) * np.cos(
                    self.theta[i]) ** 3 / NUM  # (4,-3)
                self.c18[i] = np.sqrt(9.) * np.sqrt(5. / 2.) / 4. * (
                    -3. * np.sin(2. * self.phi[i]) * np.sin(self.theta[i]) * np.cos(
                        self.theta[i]) ** 3 + 21. * np.sin(
                        2. * self.phi[i]) * np.sin(self.theta[i]) ** 2 * np.cos(self.theta[i]) ** 2) / NUM  # (4,-2)
                self.c19[i] = np.sqrt(9.) * np.sqrt(5.) / 4. * (
                    21. * np.sin(self.phi[i]) * np.cos(self.theta[i]) * np.sin(self.theta[i]) ** 3 - 9. * np.sin(
                        self.phi[i]) * np.cos(
                        self.theta[i]) * np.sin(self.theta[i])) / NUM  # (4,-1)
                self.c20[i] = np.sqrt(9.) * 3. / 64. * (
                    20. * np.cos(2. * self.theta[i]) - 35. * np.cos(4. * self.theta[i]) - 9.) / NUM  # (4,-0)
                self.c21[i] = np.sqrt(9.) * np.sqrt(5.) / 4. * (
                    21. * np.cos(self.phi[i]) * np.cos(self.theta[i]) * np.sin(self.theta[i]) ** 3 - 9. * np.cos(
                        self.phi[i]) * np.cos(
                        self.theta[i]) * np.sin(self.theta[i])) / NUM  # (4, 1)
                self.c22[i] = np.sqrt(9.) * np.sqrt(5. / 2.) / 4. * (
                    -3. * np.cos(2. * self.phi[i]) * np.sin(self.theta[i]) * np.cos(
                        self.theta[i]) ** 3 + 21. * np.cos(
                        2. * self.phi[i]) * np.sin(self.theta[i]) ** 2 * np.cos(self.theta[i]) ** 2) / NUM  # (4, 2)
                self.c23[i] = np.sqrt(9.) * np.sqrt(35.) * 3. / 4. * np.cos(3. * self.phi[i]) * np.sin(
                    self.theta[i]) * np.cos(
                    self.theta[i]) ** 3 / NUM  # (4, 3)
                self.c24[i] = np.sqrt(9.) * np.sqrt(35. / 2.) * 3. / 8. * np.cos(4. * self.phi[i]) * np.cos(
                    self.theta[i]) ** 4 / NUM  # (4, 4)

            if self.DEG >= 5:
                # 5th order
                self.c25[i] = np.sqrt(11.) * 3. / 16. * np.sqrt(77.) * np.sin(5. * self.phi[i]) * \
                         np.cos(self.theta[i]) ** 5 / NUM  # (5,-5)
                self.c26[i] = np.sqrt(11.) * 3. / 8. * np.sqrt(385. / 2.) * np.sin(4. * self.phi[i]) * \
                         np.sin(self.theta[i]) * np.cos(self.theta[i]) ** 4 / NUM  # (5,-4)
                self.c27[i] = np.sqrt(11.) * np.sqrt(385.) / 16. * (
                    9. * np.sin(3. * self.phi[i]) * np.cos(self.theta[i]) ** 3
                    * np.sin(self.theta[i]) ** 2) / NUM  # (5, -3)
                self.c28[i] = np.sqrt(11.) * np.sqrt(1155. / 2.) / 4. * (
                    3. * np.sin(2. * self.phi[i]) * np.cos(self.theta[i]) ** 2 * np.sin(
                        self.theta[i]) ** 3 - np.sin(2. * self.phi[i]) * np.cos(
                        self.theta[i]) ** 2 * np.sin(self.theta[i])) / NUM  # (5,-2)
                self.c29[i] = np.sqrt(11.) * np.sqrt(165. / 2.) / 8. * (
                    np.sin(self.phi[i]) * np.cos(self.theta[i]) - 7. * np.sin(self.phi[i]) * np.cos(
                        self.theta[i]) * np.sin(
                        self.theta[i]) ** 2 + 21. * np.sin(self.phi[i]) * np.cos(self.theta[i]) * np.sin(
                        self.theta[i]) ** 4) / NUM  # (5,-1)
                self.c30[i] = np.sqrt(11.) * np.sqrt(11.) / 8. * (
                    15. * np.sin(self.theta[i]) - 70. * np.sin(self.theta[i]) ** 3 + 63. * np.sin(
                        self.theta[i]) ** 5) / NUM  # (5, 0)
                self.c31[i] = np.sqrt(11.) * np.sqrt(165. / 2.) / 8. * (
                    np.cos(self.phi[i]) * np.cos(self.theta[i]) - 7. * np.cos(self.phi[i]) * np.cos(
                        self.theta[i]) * np.sin(
                        self.theta[i]) ** 2 + 21. * np.cos(self.phi[i]) * np.cos(self.theta[i]) * np.sin(
                        self.theta[i]) ** 4) / NUM  # (5, 1)
                self.c32[i] = np.sqrt(11.) * np.sqrt(1155. / 2.) / 4. * (
                    3. * np.cos(2. * self.phi[i]) * np.cos(self.theta[i]) ** 2 * np.sin(
                        self.theta[i]) ** 3 - np.cos(2. * self.phi[i]) * np.cos(
                        self.theta[i]) ** 2 * np.sin(self.theta[i])) / NUM  # (5, 2)
                self.c33[i] = np.sqrt(11.) * np.sqrt(385.) / 16. * (
                    9. * np.cos(3. * self.phi[i]) * np.cos(self.theta[i]) ** 3 * np.sin(
                        self.theta[i]) ** 2) / NUM  # (5, 3)
                self.c34[i] = np.sqrt(11.) * 3. / 8. * np.sqrt(385. / 2.) * np.cos(4. * self.phi[i]) * np.sin(
                    self.theta[i]) * np.cos(
                    self.theta[i]) ** 4 / NUM  # (5, 4)
                self.c35[i] = np.sqrt(11.) * 3. / 16. * np.sqrt(77.) * np.cos(5. * self.phi[i]) * np.cos(
                    self.theta[i]) ** 5 / NUM  # (5, 5)

            if self.DEG > 6:
                raise ValueError("DEG = %d is not implemented yet\n" % self.DEG)

        if self.DEG == 1:
            coeffs = np.array([self.W, self.Y, self.Z, self.X])
        elif self.DEG == 2:
            coeffs = np.array([self.W, self.Y, self.Z, self.X, self.V, self.T, self.R, self.S, self.U])
        elif self.DEG == 3:
            coeffs = np.array([self.W, self.Y, self.Z, self.X, self.V, self.T, self.R, self.S, self.U,
                               self.Q, self.O, self.M, self.K, self.L, self.N, self.P])
        elif self.DEG == 4:
            coeffs = np.array([self.W, self.Y, self.Z, self.X, self.V, self.T, self.R, self.S, self.U,
                               self.Q, self.O, self.M, self.K, self.L, self.N, self.P,
                               self.c16, self.c17, self.c18, self.c19, self.c20,
                               self.c21, self.c22, self.c23, self.c24])
        elif self.DEG == 5:
            coeffs = np.array([self.W, self.Y, self.Z, self.X, self.V, self.T, self.R, self.S, self.U,
                               self.Q, self.O, self.M, self.K, self.L, self.N, self.P,
                               self.c16, self.c17, self.c18, self.c19, self.c20,
                               self.c21, self.c22, self.c23, self.c24, self.c25,
                               self.c26, self.c27, self.c28, self.c29, self.c30,
                               self.c31, self.c32, self.c33, self.c34, self.c35])

        if self.inversion == 1:
            coeffs = np.linalg.pinv(coeffs, rcond=1e-8).T / NUM
            coeffs[abs(coeffs) < 1e-8] = 0.  # because inversion gives some very small values somewhere

        # MULTIPLYING FOR THE SELECTED DECODING SCHEME
        coeff = np.empty(coeffs.shape, dtype=np.float64)
        for i in range(self.DEG + 1):
            for jj in range(i ** 2, (i + 1) ** 2):
                coeff[jj] = g1[i] * coeffs[jj]

        self.coeff = coeff

    def _basic(self):
        g0 = [None] * (self.DEG + 1)
        g1 = [None] * (self.DEG + 1)
        for i in range(0, self.DEG + 1):
            g0[i] = 1.
            g1[i] = 1.
        g1[0] = g0[0]
        return g1

    def _maxRe(self, phi):
        NUM = len(phi)  # number of speakers
        g0 = [None] * (self.DEG + 1)
        g1 = [None] * (self.DEG + 1)

        g1p = maxReCoeffs(self.DEG)
        Egm = sum([(2. * i + 1.) * g1p[i] ** 2 for i in range(self.DEG + 1)])
        g0 = [np.sqrt(NUM / Egm)] * 4
        g1 = [g1p[i] * g0[0] for i in range(len(g1p))]

        g1[0] = g0[0]
        return g1

    def _phase(self, phi):
        NUM = len(phi)  # number of speakers
        g0 = [None] * (self.DEG + 1)
        g1 = [None] * (self.DEG + 1)

        g0d = np.sqrt(3. * NUM / 4.)  # from Dani 1st order Ambisonics
        g1d = g0d * 1. / 3.  # from Dani 1st order Ambisonics
        g1p = [None] * (self.DEG + 1)

        for i in range(0, self.DEG + 1):
            g0[i] = np.sqrt(NUM * (2 * self.DEG + 1)) / (self.DEG + 1)
            g1[i] = g0[i] * (mh.factorial(self.DEG) * mh.factorial(self.DEG + 1)) / float(
                mh.factorial(self.DEG + i + 1) * mh.factorial(self.DEG - i))  # note that g1(0) is equal to g0
            g1p[i] = (mh.factorial(self.DEG) * mh.factorial(self.DEG + 1)) / float(
                mh.factorial(self.DEG + i + 1) * mh.factorial(self.DEG - i))  # just for debugging purposes

        g1[0] = g0[0]
        return g1


class Support:
    def __init__(self, configClass, ambisoniClass):
        self.cfg = configClass
        self.amb = ambisoniClass
        self.nD  = self.cfg.nD
        self.DEC = self.cfg.DEC
        self.DEG = self.cfg.DEG
        self.PHI = self.cfg.PHI
        self.THETA  = self.cfg.THETA
        self.NSPK   = self.cfg.NSPK
        self.NPOINTS    = self.cfg.NPOINTS
        self.MATCHSPK   = self.cfg.MATCHSPK
        self.NSPKmatch  = self.cfg.NSPKmatch
        self.coeffDir   = None

    ##########################
    ## supporting functions ##
    ##########################
    def Sij(self, coeffSpk, coeffDir, NSphPt):
        NSPK = self.cfg.NSPK
        try:
            sij = np.dot(coeffSpk.T, coeffDir * NSphPt)  # this will have the dimensions of NSPK*NPOINTS
        except ValueError:
            raise ValueError("Wrong dimensions when doing dot product in sij")

        if sij.shape != (NSPK, NSphPt):
            (a, b) = sij.shape
            raise ValueError("Wrong dimensions in Sij. Should be (%i,%i) but it's (%i,%i) \n" % (NSPK, NSphPt, a, b))
        return sij

    def physOmni(self, Sij):
        # LOW FREQUENCIES
        # pressure
        pressure = sum(Sij)
        # velocity
        V = self.velocity(Sij)

        # HIGH FREQUENCIES
        # energy density
        energyD = sum(Sij * Sij)
        # intensity
        J = self.velocity(Sij * Sij)

        return pressure, V, energyD, J

    def velocity(self, Sij):
        PHI     = self.cfg.PHI
        THETA   = self.cfg.THETA
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
            raise ValueError("Something went wrong calculating velocity\n")

        return Vx, Vy, Vz

    @staticmethod
    def Radial(A, B):
        Sum = 0
        if len(A) != len(B): raise ValueError("I'm dying in Radial function. Arrays with different shapes.")
        for i in range(len(A)):
            Sum += A[i] * B[i]

        return Sum

    @staticmethod
    def Tang(A, B):  # aka Cross
        Sum = 0
        if len(A) != len(B): raise ValueError("I'm dying in Tang function. Arrays with different shapes.")
        for i in range(len(A)):
            Sum += (A[i] * B[i - 1] - A[i - 1] * B[i]) ** 2.0

        return np.sqrt(Sum)

    def physDir(self, Sij, phi, theta):
        phi = np.asarray(phi)
        theta = np.asarray(theta)
        Zx = np.cos(phi) * np.cos(theta)
        Zy = np.sin(phi) * np.cos(theta)
        Zz = np.sin(theta)
        # FIXME: rewrite these five lines using a call to ambisoniC (first order)
        # and fixing the coefficient of X,Y,Z

        Z = Zx, Zy, Zz

        press, V, energyD, J = self.physOmni(Sij)

        Vradial = self.Radial(V, Z)
        Jradial = self.Radial(J, Z)

        Vtangential = self.Tang(V, Z)
        Jtangential = self.Tang(J, Z)

        return press, V, energyD, J, Vradial, Jradial, Vtangential, Jtangential

    def oppGain(self, Sij):
        # FOR IN-PHASE DECODING
        oppGain = np.zeros((self.NSPK, self.NPOINTS), dtype=np.float64)
        oppGain[Sij < 0] = Sij[Sij < 0]
        if oppGain.shape != (self.NSPK, self.NPOINTS): raise ValueError("Wrong dimensions in oppGain\n")
        oppGain = sum(oppGain * oppGain)

        return oppGain

    @staticmethod
    def sph2cart(phi, theta):
        """Acoustics convention!"""
        x = np.cos(phi) * np.cos(theta)
        y = np.sin(phi) * np.cos(theta)
        z = np.sin(theta)
        return x, y, z

    @staticmethod
    def Physph2cart(phi, theta):
        """Physics convention!"""
        x = np.cos(phi) * np.sin(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(theta)
        return x, y, z

    @staticmethod
    def PlotOverSphere(phi, theta, rho):
        """Acoustics convention!"""
        x = np.cos(phi) * np.cos(theta) * rho
        y = np.sin(phi) * np.cos(theta) * rho
        z = np.sin(theta) * rho
        return x, y, z

    @staticmethod
    def eval_grad(f, theta):
        theta = algopy.UTPM.init_jacobian(theta)
        return algopy.UTPM.extract_jacobian(f(theta))

    @staticmethod
    def eval_hess(f, theta):
        theta = algopy.UTPM.init_hessian(theta)
        return algopy.UTPM.extract_hessian(len(theta), f(theta))

    def vtomat(self, vector):
        if self.nD == "2D":
            tmp = vector.reshape((self.DEG * 2 + 1), self.NSPKmatch)
        elif self.nD == "3D":
            tmp = vector.reshape(((self.DEG + 1) ** 2, self.NSPKmatch))
        else:
            raise ValueError("Failing during conversion from vector to matrix.")
        return tmp

    def mattov(self, mat):
        if self.nD == "2D":
            tmp = mat.reshape((1, (self.DEG * 2 + 1) * self.NSPKmatch))[0]
        elif self.nD == "3D":
            tmp = mat.reshape((1, (self.DEG + 1) ** 2 * self.NSPKmatch))[0]
        else:
            raise ValueError("Failing during conversion from matrix to vector.")
        return tmp

    def zeroing_bounds(self, nodes, upbound, lowbound, run):
        if self.nD == "3D":
            matrix_size = (self.DEG + 1) ** 2 * self.NSPKmatch
            nodes_bool_vector = np.asarray(abs(nodes) < 3. * 10e-4).reshape(1, matrix_size)
            for i in range(matrix_size):
                if nodes_bool_vector[0, i]:
                    upbound[i] = 0.
                    lowbound[i] = 0.  # putting to zero the speakers that are in a node of a SH

        if self.nD == "2D":
            matrix_size = (self.DEG * 2 + 1) * self.NSPKmatch
            nodes_bool_vector = np.asarray(abs(nodes) < 3. * 10e-4).reshape(1, matrix_size)
            for i in range(matrix_size):
                if nodes_bool_vector[0, i]:
                    upbound[i] = 0.
                    lowbound[i] = 0.  # putting to zero the speakers that are in a node of a SH

        return nodes, upbound, lowbound


    ################################
    # The function to be minimized #
    ################################
    def function(self, VarCoeffSpk, *args):
        if len(args) > 1:
            grad = args[0]
            logging.error("You have to implement the gradient or use another algorithm")

        shapeCoeff = np.shape(VarCoeffSpk)
        # this because the minimization library destroys the shape of VarCoeffSpk
        # in 3D
        if self.cfg.nD == "3D":
            if (shapeCoeff == ((self.cfg.DEG + 1) ** 2 * self.cfg.NSPKmatch,)):
                VarCoeffSpk = self.vtomat(VarCoeffSpk)
            elif (shapeCoeff == ((self.cfg.DEG + 1) ** 2 * self.cfg.NSPK,)):
                raise ValueError("Here the dimension is " + str(shapeCoeff) + " while it should be " + str(
                    ((self.cfg.DEG + 1) ** 2 * self.cfg.NSPKmatch,)))
            elif (shapeCoeff != ((self.cfg.DEG + 1) ** 2, self.cfg.NSPKmatch)) and (shapeCoeff != ((self.cfg.DEG + 1) ** 2, self.cfg.NSPK)):
                raise ValueError("Strange dimensions of VarCoeffSpk in -function-." + str(shapeCoeff))

            if self.cfg.MATCHSPK: VarCoeffSpk = self.upscalematrix(VarCoeffSpk, self.cfg.MATCHED)
            if (np.shape(VarCoeffSpk) != ((self.cfg.DEG + 1) ** 2, self.cfg.NSPK)):
                raise ValueError("Here the matrix should be ((DEG+1)**2,NSPK)")

        # in 2D
        if self.cfg.nD == "2D":
            if (shapeCoeff == ((self.cfg.DEG * 2 + 1) * self.cfg.NSPKmatch,)):
                VarCoeffSpk = self.vtomat(VarCoeffSpk)
            elif (shapeCoeff == ((self.cfg.DEG * 2 + 1) * self.cfg.NSPK,)):
                raise ValueError(
                    "Here the dimension is " + str(shapeCoeff) + " while it should be " + str(
                        ((self.cfg.DEG * 2 + 1) * self.cfg.NSPKmatch,)))
            elif (shapeCoeff != ((self.cfg.DEG * 2 + 1), self.cfg.NSPKmatch)) and (shapeCoeff != ((self.cfg.DEG * 2 + 1), self.cfg.NSPK)):
                raise ValueError("Strange dimensions of VarCoeffSpk in -function-." + str(shapeCoeff))

            if self.cfg.MATCHSPK: VarCoeffSpk = self.upscalematrix(VarCoeffSpk, self.cfg.MATCHED)
            if (np.shape(VarCoeffSpk) != ((self.cfg.DEG * 2 + 1), self.cfg.NSPK)):
                raise ValueError("Here the matrix should be ((DEG*2+1),NSPK): (%d,%d)" % ((self.cfg.DEG * 2 + 1), self.cfg.NSPK))


        """The function to be minimized"""
        Wj = np.ones(self.coeffDir.shape[1])  # biasing factor dependent on direction j
        Mj = np.ones(self.coeffDir.shape[1])  # biasing factor dependent on direction j

        sij = self.Sij(VarCoeffSpk, self.coeffDir, self.cfg.NPOINTS)
        pressure, V, energyD, J, Vradial, Jradial, Vtang, Jtang = self.physDir(sij, self.cfg.phiTest, self.cfg.thetaTest)

        # weighting functions
        if self.cfg.WFRONT: Wj = Wj * self.cfg.WfrontVec  # maybe it's better if you calculate Wfront Wplane and Wbinary outside the function, so that you calculate them only once...
        if self.cfg.WPLANE: Wj = Wj * self.cfg.WplaneVec
        if self.cfg.WBIN:   Wj = Wj * self.cfg.WbinVec
        if self.cfg.WAUTOREM: Mj = Mj * self.cfg.WremVec
        if self.cfg.WFRONT or self.cfg.WPLANE or self.cfg.WBIN: Wj = Wj * float(len(Wj)) / sum(Wj)

        if self.cfg.DEC == 'basic':
            Tpressure = ((1. - pressure) ** 2 * Wj * Mj) / self.cfg.NPOINTS
            TVlon = ((1. - Vradial) ** 2 * Wj * Mj) / self.cfg.NPOINTS
            TVtang = ((Vtang) ** 2 * Wj * Mj) / self.cfg.NPOINTS
            Tvar = 0
            Tpressure2 = 0
            TVlon2 = 0
            TVtang2 = 0

            TenergyD = ((1. - energyD) ** 2 * Wj * Mj) / self.cfg.NPOINTS
            # TenergyD = ((1.-energyD)**2*Wj*Mj*2)/self.cfg.NPOINTS
            TJrad = ((1. - Jradial) ** 2 * Wj * Mj) / self.cfg.NPOINTS
            TJtang = ((Jtang) ** 2 * Wj * Mj) / self.cfg.NPOINTS

            if self.cfg.PREFHOM:
                Tvar = np.var(VarCoeffSpk[0]) / (np.mean(VarCoeffSpk[0])) ** 2

            if (self.cfg.WAUTOREM or self.cfg.WBIN):
                Tpressure2 = (~self.cfg.WremVec) * ((0.5 - pressure) ** 2 * Wj) / self.cfg.NPOINTS  # Pressure -3 dB
                TVlon2 = (~self.cfg.WremVec) * ((1. - Vradial) ** 2 * Wj) / self.cfg.NPOINTS
                TVtang2 = (~self.cfg.WremVec) * ((Vtang) ** 2 * Wj) / self.cfg.NPOINTS

            target = np.sum(self.cfg.CP * (Tpressure + Tpressure2) +
                            + self.cfg.CR * (TVlon + 0.1 * TVlon2)
                            + self.cfg.CT * (TVtang + 0.1 * TVtang2)
                            + self.cfg.CPH * TJrad + self.cfg.CPH * TJtang + self.cfg.CPH * TenergyD) + self.cfg.CV * Tvar

        elif self.cfg.DEC == 'maxRe' or self.cfg.DEC == 'mixed':
            TenergyD = ((1. - energyD) ** 2 * Wj * Mj) / self.cfg.NPOINTS
            TJrad = ((1. - Jradial) ** 2 * Wj * Mj) / self.cfg.NPOINTS
            TJtang = ((Jtang) ** 2 * Wj * Mj) / self.cfg.NPOINTS
            Tvar = 0
            TJrad2 = 0
            TenergyD2 = 0
            TJtang2 = 0
            if self.cfg.PREFHOM:
                Tvar = np.var(VarCoeffSpk[0]) / (np.mean(VarCoeffSpk[0])) ** 2
            if (self.cfg.WAUTOREM or self.cfg.WBIN):
                TenergyD2 = (~self.cfg.WremVec) * ((0.5 - energyD) ** 2 * Wj) / self.cfg.NPOINTS  # Energy -3 dB
                TJrad2 = (~self.cfg.WremVec) * ((1. - Jradial) ** 2 * Wj) / self.cfg.NPOINTS
                TJtang2 = (~self.cfg.WremVec) * ((Jtang) ** 2 * Wj) / self.cfg.NPOINTS
                # ~array : negation of an array of bools

            target = np.sum(
                self.cfg.CE * (TenergyD + TenergyD2) + self.cfg.CR * (TJrad + 0.1 * TJrad2) + self.cfg.CT * (TJtang + 0.1 * TJtang2)) + self.cfg.CV * Tvar

        elif self.cfg.DEC == 'phase':
            TenergyD = ((1. - energyD) ** 2 * Wj) / self.cfg.NPOINTS
            TJrad = ((1. - Jradial) ** 2 * Wj * Mj) / self.cfg.NPOINTS
            TJtang = ((Jtang) ** 2 * Wj * Mj) / self.cfg.NPOINTS
            Tvar = 0
            if self.cfg.PREFHOM:
                Tvar = np.var(VarCoeffSpk[0]) / (np.mean(VarCoeffSpk[0])) ** 2
            ToppGain = np.linalg.norm(self.oppGain(sij)) ** 2 / self.cfg.NPOINTS

            target = np.sum(self.cfg.CE * TenergyD + self.cfg.CR * TJrad + self.cfg.CT * TJtang) \
                            + self.cfg.CPH * ToppGain * self.cfg.CV * Tvar
            # FIXME: missing some extra factors (see dani's paper page 5)
        else:
            target = 0
            raise ValueError("The decoding you chose is not implemented.")
        return target

    # with nlopt you can also write linear and non linear constraints... have a look at the reference
    # http://ab-initio.mit.edu/wiki/index.php/NLopt_Python_Reference
    def eqconstr(self, result, x, grad):
        if grad.size > 0:
            print "gradient to be implemented"

        nodes = self.amb.get_ambisonic_coeff(self.PHI, self.THETA, 'basic', self.DEG, 0)  # calculates the coefficients in the basic or naive scheme.
        # These values are used to figure out which speakers lay in the nodes of some spherical harmonic
        for i in range((self.DEG + 1) ** 2 * self.NSPK):
            if np.asarray(abs(nodes) < 10e-8).reshape(1, ((self.DEG + 1) ** 2 * self.NSPK))[0, i]:
                result[i] = x[i]  # putting to zero the speakers that are in a node of a SH
        return

    def inconstr(self, result, x, grad):
        if grad.size > 0:
            print "gradient to be implemented"

        nodes = self.amb.get_ambisonic_coeff(self.PHI, self.THETA, 'basic', self.DEG, 0)  # calculates the coefficients in the basic or naive scheme.
        # These values are used to figure out which speakers lay in the nodes of some spherical harmonic
        for i in range((self.DEG + 1) ** 2 * self.NSPK):
            if np.asarray(nodes > 0.).reshape(1, ((self.DEG + 1) ** 2 * self.NSPK))[0, i]: result[i] > x[i]  # keeping the sign
            if np.asarray(nodes < 0.).reshape(1, ((self.DEG + 1) ** 2 * self.NSPK))[0, i]: result[i] < x[i]  # keeping the sign
        return

    def downscalematrix(self, CoeffMat, matched):
        '''
        * here we need two things:
        * - the matrix of coefficients that we want to downscale
        * - the vector of paired speakers
        *
        * The downscaled (shrinked) matrix will be our output here
        '''
        ResShrinked = CoeffMat.T
        for idx, spki in enumerate(matched):
            if (spki >= 0):
                # the new matrix has to have the idx rows only
                if (idx == 0):
                    ResShrinked = CoeffMat.T[idx]
                else:
                    ResShrinked = np.vstack((ResShrinked, CoeffMat.T[idx]))
        if self.MATCHSPK:
            return ResShrinked.T
        else:
            return CoeffMat

    def upscalematrix(self, ResShrinked, matched):
        """
        * here we need three things:
        * - the matrix of coefficients that we want to downscale
        * - the vector of paired speakers
        * - the vector of SH signs that change with left/right symmetry
        *
        * The upscaled matrix will be our output here
        """

        ResUp = ResShrinked.T
        if (np.shape(ResShrinked) == ((self.DEG + 1) ** 2, self.NSPK)) or (np.shape(ResShrinked) == ((self.DEG * 2 + 1), self.NSPK)):
            return ResShrinked

        else:
            counter = 0  # you don't strictly need a counter, you can use len(shrinked) ... I guess

            for idx, spki in enumerate(matched):
                if (idx == 0):
                    ResUp = ResShrinked.T[counter]
                    counter += 1
                    continue

                if (spki >= 0):
                    ResUp = np.vstack((ResUp, ResShrinked.T[counter]))
                    counter += 1

                elif (spki < 0):
                    # search for the matched speaker
                    # (which is the one labelled with positive spki)
                    positive_match = matched[matched >= 0]
                    (wanted_idx, wanted_val) = [(jdx, val) for jdx, val in enumerate(positive_match) if val == abs(spki)][0]
                    ResUp = np.vstack((ResUp, ResShrinked.T[wanted_idx] * self.cfg.SIGNSvec))  # multiply with vector of signs
                else:
                    raise ValueError("Ehm ... this is impossible! [upscalematrix function]")

            if self.MATCHSPK:
                return ResUp.T
            else:
                return ResShrinked
