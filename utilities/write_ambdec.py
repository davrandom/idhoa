'''
* This file is part of IDHOA software.
*
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

import argparse
import pickle
import numpy as np
import os
import sys
lib_path = os.path.abspath(os.path.join(os.path.dirname(__file__),".."))
print lib_path
sys.path.append(lib_path)

from auxiliary import Conventions
from plotting import threed_polar_plot


# Parsing arguments
parser = argparse.ArgumentParser(description='Generate Ambdec config file from "idhoa" files.')
parser.add_argument('-LF', '--lf-matrix', type=str, dest='lf_results_file', default='utilities/results_for_ccrma_2_basic.idhoa',
                    help='Path to "idhoa" results file generate LF decoding matrix. ')
parser.add_argument('-pinv', '--use_pinv_for_LF', dest='pinv', action='store_true',
                    help='Will use pinv guess for low frequencies instead of the NL minimum. ')
parser.add_argument('-HF', '--hf-matrix', type=str, dest='hf_results_file', default='utilities/results_for_ccrma_2_phase.idhoa',
                    help='Path to "idhoa" results file generate HF decoding matrix.')
parser.add_argument('-jr', '--jack_routing', type=str, dest='speaker_order', default='',
                    help='Change the order of speakers from the one you entered in ini file. '
                         'A possible string could be "4 1 3 2 5"')
parser.add_argument('-os', '--output_scale', type=str, dest='output_scale', default='N3D',
                    help='Default normalization is N3D. You can ask to change it by specifying '
                         'the destination convention. '
                         '\nAvailable: "N3D" "SN3D" "MaxN" "FuMa" "SN2D" "N2D" ')
parser.add_argument('-is', '--input_scale', type=str, dest='input_scale', default='fuma',
                    help='Default normalization is N3D. You can ask to change it by specifying '
                         'the destination convention. '
                         '\nAvailable: "N3D" "SN3D" "MaxN" "FuMa" "SN2D" "N2D" ')
parser.add_argument('-of', '--output_filename', type=str, dest='output_filename', default='idhoa.ambdec',
                    help='Name of the Ambdec configuration file.')

args = parser.parse_args()


# threed_polar_plot(az,el,ra)
# threed_polar_plot(az,el,1)


# Reading ini configurations files
LF = pickle.load(open(args.lf_results_file, 'rb'))
HF = pickle.load(open(args.hf_results_file, 'rb'))

if args.speaker_order != '' :
    spkord = args.speaker_order.split(' ')
    spkord = [int(j) for i, j in enumerate(spkord)]
else :
    spkord = range(1, len(LF.PHI)+1)


# Writing ambdec configuration file
outf = open(args.output_filename, 'w')

if LF.DEG != HF.DEG:
    raise ValueError("The two idhoa files produce decodings at different orders (LF:%d, HF:%d)."
                     "It is not possible to combine them in a single Ambdec configuration file." % (LF.DEG, HF.DEG))
if LF.nD != HF.nD:
    raise ValueError("The two idhoa files produce decodings with different dimensionality (LF:%s, HF:%s)."
                     "It is not possible to combine them in a single Ambdec configuration file." % (LF.nD, HF.nD))
if LF.NSPK != HF.NSPK:
    raise ValueError("The two idhoa files produce decodings for two different speakers' layout (LF:%d, HF:%d)."
                     "It is not possible to combine them in a single Ambdec configuration file." % (LF.NSPK, HF.NSPK))

DEG = HF.DEG
nD  = HF.nD
NSPK= HF.NSPK

lf_mat = LF.obj_minimization_matrix
hf_mat = HF.obj_minimization_matrix
if args.pinv:
    lf_mat = LF.pseudoinv_guess_matrix

# from rad to grad
LF.grad_az = [180.0 * ii / np.pi for ii in LF.PHI]
LF.grad_el = [180.0 * ii / np.pi for ii in LF.THETA]
HF.grad_az = [180.0 * ii / np.pi for ii in HF.PHI]
HF.grad_el = [180.0 * ii / np.pi for ii in HF.THETA]
try:
    hflabel = HF.label
    lflabel = LF.label
except:
    LF.label = [str(i + 1) for i, j in enumerate(LF.PHI)]

try:
    lfradius = LF.radius
    hfradius = HF.radius
except:
    LF.radius = np.ones(len(LF.PHI))

# convention change (if needed)
conv = Conventions(DEG)
# TODO
# e.g.# hf_mat = ResCoeff.T*conv.shrink(conv.fuma.n2d)

mask = ""
if nD == '2D':
    if DEG == 1:
        mask = "b"
    if DEG == 2:
        mask = "11b"
    if DEG == 3:
        mask = "831b"
    if DEG > 3:
        mask = "dontknow!"
elif nD == '3D':
    if DEG == 1:
        mask = "f"
    if DEG == 2:
        mask = "1ff"
    if DEG == 3:
        mask = "ffff"
else:
    raise ValueError("The dimensionality of the problem is unknown: %s. Should be 2D or 3D." % nD)

h1 = '''# AmbDec configuration
# Written by write_ambdec.py

/description     5.1 idhoa decoder 

/version          3

/dec/chan_mask    %s
/dec/freq_bands   2
/dec/speakers     %d
/dec/coeff_scale  n3d

/opt/input_scale  %s
/opt/nfeff_comp   input
/opt/delay_comp   on
/opt/level_comp   off
/opt/xover_freq    400
/opt/xover_ratio   0.0

/speakers/{ ''' % (mask, NSPK, args.input_scale)

print h1
outf.write(h1+"\n")

string = ""
for jj,spkname in enumerate(spkord):
    print "add_spkr    %s    %.3f    %.5f    %.5f    system:playback_%d" % (LF.label[jj], LF.radius[jj], LF.grad_az[jj], LF.grad_el[jj], spkname)
    string = "add_spkr    %s    %.3f    %.5f    %.5f    system:playback_%d" % (LF.label[jj], LF.radius[jj], LF.grad_az[jj], LF.grad_el[jj], spkname)
    outf.write(string+"\n")
    string = ""

h2 = '''/}

/lfmatrix/{
order_gain     1.00000  1.00000  1.00000  1.00000'''
print h2
outf.write(h2+"\n")

string = ""
for row in range(len(lf_mat)):
    for column,val in  enumerate(lf_mat[row]):
        string += "   %.6f" % val
    print "add_row " +string
    outf.write("add_row " +string+"\n")
    string = ""

h3 = '''/}

/hfmatrix/{
order_gain     1.00000  1.00000  1.00000  1.00000'''
print h3
outf.write(h3+"\n")

for row in range(len(hf_mat)):
    for column, val in  enumerate(hf_mat[row]):
        string += "   %.6f" % val
    print "add_row " +string
    outf.write("add_row " +string+"\n")
    string = ""


h4 = '''/}

/end'''
print h4
outf.write(h4+"\n")
