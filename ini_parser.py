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

from ConfigParser import SafeConfigParser, ParsingError
import json
import numpy as np

import auxiliary as aux


class ConfigConstants:
    def __init__(self, configfile):
        self.configfile = configfile
        self._parse()
        self.projection_guess_matrix = None
        self.pseudoinv_guess_matrix = None
        self.obj_minimization_matrix = None

    def _parse(self):
        # Require values
        try:
            parser = SafeConfigParser()
            parser.read(self.configfile)

            self.case = parser.get('Layout', 'name')
            self.PHI = json.loads(parser.get('Layout', 'PHI'))
            self.THETA = json.loads(parser.get('Layout', 'THETA'))

            self.matlab_filename = parser.get('Outfiles', 'matlab_filename')
            self.is_mat_out_file = parser.get('Outfiles', 'save_matlab_for_ambisonics_toolbox')

            self.DEC = parser.get('Ambisonics', 'DEC')
            self.DEG = parser.getint('Ambisonics', 'DEG')
            self.nD = parser.get('Ambisonics', 'nD')

            self.mute_small_coeffs = parser.getboolean('Flags', 'mute_small_coeffs')
            self.autoexclude_regions_with_no_spkrs_binary = parser.getboolean('Flags', 'autoexclude_regions_with_no_spkrs_binary')
            self.prefer_homogeneous_coeffs = parser.getboolean('Flags', 'prefer_homogeneous_coeffs')
            self.prefer_front = parser.getboolean('Flags', 'prefer_front')
            self.prefer_horiz_plane = parser.getboolean('Flags', 'prefer_horiz_plane')
            self.exclude_with_theta_binary_mask = parser.getboolean('Flags', 'exclude_with_theta_binary_mask')
            self.autoexclude_regions_with_no_spkrs_smooth = parser.getboolean('Flags', 'autoexclude_regions_with_no_spkrs_smooth')
            self.thetaThreshold = parser.getint('Flags', 'thetaThreshold')
            self.match_symmetric_spkrs = parser.getboolean('Flags', 'match_symmetric_spkrs')
            self.match_tolerance = parser.getint('Flags', 'match_tolerance')

            self.CP = parser.getint('Minimiz', 'CP')
            self.CV = parser.getint('Minimiz', 'CV')
            self.CE = parser.getint('Minimiz', 'CE')
            self.CR = parser.getint('Minimiz', 'CR')
            self.CT = parser.getint('Minimiz', 'CT')
            self.CPH = parser.getint('Minimiz', 'CPH')

            self.SEED = parser.getint('Other', 'SEED')

        except ParsingError, err:
            print 'Could not parse:', err

        # number of speakers
        self.NSPK = len(self.PHI)

        self.PHI, self.THETA = aux.fix_phi(self.PHI, self.THETA)

        # set: phiTest, thetaTest, NPOINTS, WfrontVec, WplaneVec, WbinVec, WremVec
        self._autoinit()

        self.SIGNSvec = aux.SHsigns(self.DEG)
        if self.match_symmetric_spkrs:
            self.MATCHED, self.PHI, self.THETA = aux.analyzeLayout(self.PHI, self.THETA, self.match_tolerance * np.pi / 180)
            self.NSPKmatch = len(self.MATCHED[self.MATCHED >= 0]);
            self.MAP = aux.MapMatched(self.MATCHED)
        else:
            self.NSPKmatch = self.NSPK
            self.MATCHED = np.zeros(self.NSPK)

    def _autoinit(self):
        ###############################
        #  automatic initializations  #
        ###############################
        # Threshold for binary masking the lower region
        self.thetaThreshold *= np.pi / 180.0

        self.phiTest, self.thetaTest = aux.get_points_over_sphere(self.SEED, self.nD)  # generating sampling space points

        if self.autoexclude_regions_with_no_spkrs_binary:
            self.phiTest, self.thetaTest = self._autoremoval()  # autoremoving test points without speakers around

        self.WfrontVec = self._Wfront(self.phiTest, self.thetaTest)
        self.WplaneVec = self._Wplane(self.thetaTest)
        self.WbinVec = self._Wbinary(self.thetaTest, self.thetaThreshold)
        self.WremVec = self._Wautoremoval()  # only necessary if autoexclude_regions_with_no_spkrs_smooth is true...
        self.NPOINTS = len(self.phiTest)

    #############
    #  weights  #
    #############

    @staticmethod
    def _Wfront(phi, the):
        # wvect = 1. + np.cos(phi) * np.cos(the)/2.
        wvect = np.cos(phi) ** 8 * np.cos(the) ** 8 / 2.
        return wvect

    @staticmethod
    def _Wplane(the):
        wvect = 1. + np.cos(the) ** 2.
        return wvect

    @staticmethod
    def _Wbinary(the, theta_threshold):
        return the > theta_threshold

    def _spkrs_distance(self):
        dvec = np.array([])
        mins = np.array([])
        for i, valpi in enumerate(self.PHI):
            for j, valpj in enumerate(self.PHI):  # you could do probably something like: for j in range(i+1,len(PHI))
                dvec = np.append(dvec, [aux.angular_distance(valpi, self.THETA[i], valpj, self.THETA[j])])

            mins = np.append(mins, [min(dvec[dvec != 0])])
            dvec = np.array([])  # reset dvec
        mean = np.mean(mins)  # calculates the mean only of the smalles values - the closest speakers
        return mean

    def _autoremoval(self):
        mean_spk_dist = self._spkrs_distance()
        phit = []
        thetat = []
        for i, valpi in enumerate(self.phiTest):
            for j, valpj in enumerate(self.PHI):
                if aux.angular_distance(self.phiTest[i], valpi, valpj, self.THETA[j]) < mean_spk_dist * 1.0:
                    phit.append(self.phiTest[i])
                    thetat.append(self.thetaTest[i])
                    break
        return phit, thetat

    def _Wautoremoval(self):
        mean_spk_dist = self._spkrs_distance()
        wvect = []
        for i, valpi in enumerate(self.phiTest):
            for j, valpj in enumerate(self.PHI):
                if aux.angular_distance(self.phiTest[i], valpi, valpj, self.THETA[j]) < mean_spk_dist * 1.0:
                    temp = True
                    # temp = 1
                    break
                else:
                    temp = False
                    # temp = 0.3
            wvect.append(temp)
        return np.asarray(wvect)
