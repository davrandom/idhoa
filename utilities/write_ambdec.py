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


import os,sys
lib_path = os.path.abspath(os.path.join(os.path.dirname(__file__),".."))
print lib_path
sys.path.append(lib_path)

from auxiliary import cart2sph, fix_phi, Conventions
from plotting import threed_polar_plot
import numpy as np
from ConfigParser import SafeConfigParser, ParsingError
import json
import argparse
from numpy import array


class Ini_Parser:
    def __init__(self,ini_configfile):
        try:
            parser = SafeConfigParser()
            parser.read(ini_configfile)
            
            self.case = parser.get('Layout','name')
            self.az = json.loads(parser.get('Layout', 'PHI'))
            self.el = json.loads(parser.get('Layout', 'THETA'))
            try:
                self.ra = json.loads(parser.get('Layout', 'radius'))
            except:
                self.ra = np.ones(len(self.az))
            try:
                self.label = json.loads(parser.get('Layout', 'chlabel'))
            except:
                self.label = [str(i+1) for i,j in enumerate(self.az)]
            
            self.DEC = parser.get('Ambisonics','DEC')
            self.DEG = parser.getint('Ambisonics','DEG')
          
            self.case = parser.get('Layout','name')
            self.autoexclude_regions_with_no_spkrs_binary =  parser.getboolean('Flags','autoexclude_regions_with_no_spkrs_binary')
            self.match_symmetric_spkrs = parser.getboolean('Flags','match_symmetric_spkrs')
            self.CP = parser.getint('Minimiz','CP') 
            self.CV = parser.getint('Minimiz','CV')
            self.CE = parser.getint('Minimiz','CE')
            self.CR = parser.getint('Minimiz','CR')
            self.CT = parser.getint('Minimiz','CT')
            self.CPH= parser.getint('Minimiz','CPH')
        
        except ParsingError, err:
            print 'Could not parse:', err
        
        print "Azimut and elevation of ", ini_configfile, " please check everything is fine."
        self.az, self.el = fix_phi(self.az,self.el)
        print "\nradius (cm), elevation, azimut (radiants)"
        print self.ra
        print self.el
        print self.az
        
        # from rad to grad
        self.gaz = [180.0*ii/np.pi for ii in self.az]
        self.gel = [180.0*ii/np.pi for ii in self.el]
        
        print "\nradius (cm), elevation, azimut (grads)"
        print self.ra
        print self.gel
        print self.gaz
        print self.label
        
        #threed_polar_plot(az,el,ra)
        #threed_polar_plot(az,el,1)




# MAIN
# Parsing arguments
parser = argparse.ArgumentParser(description='Generate Ambdec config file from idhoa ini and py files.')
parser.add_argument('-L', '--lf-matrix', type=str, dest='lf_ini_configfile', default='init_files/example.ini', 
                    help='Path to .ini init file used to generate LF decoding matrix. It assumes that the .py file with results is located in main idhoa folder.')
parser.add_argument('-H', '--hf-matrix', type=str, dest='hf_ini_configfile', default='init_files/example.ini',  
                    help='Path to .ini init file used to generate HF decoding matrix. It assumes that the .py file with results is located in main idhoa folder.')
parser.add_argument('-r', '--jack_routing', type=str, dest='speaker_order', default='',
                    help='Change the order of speakers from the one you entered in ini file. A possible string could be "4 1 3 2 5"')
parser.add_argument('-c', '--output_convention', type=str, dest='convention', default='N3D',
                    help='Default normalization is N3D. You can ask to change it by specifying the destination convention. \nAvailable: "N3D" "SN3D" "MaxN" "FuMa" "SN2D" "N2D" ')
parser.add_argument('-o', '--output_filename', type=str, dest='output_filename', default='idhoa.ambdec',
                    help='Name of the Ambdec configuration file.')

args = parser.parse_args()

input_scale = "fuma" # possible ones: fuma, n3d

# Reading ini configurations files
LF = Ini_Parser(args.lf_ini_configfile)
HF = Ini_Parser(args.hf_ini_configfile)

if args.speaker_order != '' :
    spkord = args.speaker_order.split(' ')
    spkord = [int(j) for i,j in enumerate(spkord)]
else :
    spkord = range(1,len(LF.az)+1)


# Writing ambdec configuration file
outf = open(args.output_filename, 'w')

if LF.DEG != HF.DEG :
    raise ValueError("The two ini files produce decodings at different orders. It is not possible to combine them in a single Ambdec configuration file, sorry.")
    
# starting with LF
lf_pyfilename = LF.case+"-"+str(LF.DEG)+"-"+str(LF.DEC)+"-rem"+str(LF.autoexclude_regions_with_no_spkrs_binary)+"-sym"+str(LF.match_symmetric_spkrs)

if (LF.DEC=="basic"): lf_pyfilename += "CP"+str(LF.CP)+"CV"+str(LF.CV)+".py" 
if (LF.DEC=="maxRe"): lf_pyfilename += "CR"+str(LF.CR)+"CT"+str(LF.CT)+"CE"+str(LF.CE)+".py"
if (LF.DEC=="phase"): lf_pyfilename += "CR"+str(LF.CR)+"CT"+str(LF.CT)+"CPH"+str(LF.CPH)+".py"
execfile(lf_pyfilename)
lf_mat = ResCoeff.T

# then HF
hf_pyfilename = HF.case+"-"+str(HF.DEG)+"-"+str(HF.DEC)+"-rem"+str(HF.autoexclude_regions_with_no_spkrs_binary)+"-sym"+str(HF.match_symmetric_spkrs)

if (HF.DEC=="basic"): hf_pyfilename += "CP"+str(HF.CP)+"CV"+str(HF.CV)+".py" 
if (HF.DEC=="maxRe"): hf_pyfilename += "CR"+str(HF.CR)+"CT"+str(HF.CT)+"CE"+str(HF.CE)+".py"
if (HF.DEC=="phase"): hf_pyfilename += "CR"+str(HF.CR)+"CT"+str(HF.CT)+"CPH"+str(HF.CPH)+".py"
execfile(hf_pyfilename)
hf_mat = ResCoeff.T


# convention change (if needed)
conv = Conventions(DEG)
# TODO
#e.g.# hf_mat = ResCoeff.T*conv.shrink(conv.fuma.n2d)


if nD == '2D':
    if DEG == 1:
        mask = "b"
    if DEG == 2:
        mask = "11b"
    if DEG == 3:
        mask = "831b"
    if DEG > 3:
        mask = "dontknow!"

if nD.replace(" ","") == '3D':
    if DEG == 1:
        mask = "f"
    if DEG == 2:
        mask = "1ff"
    if DEG == 3:
        mask = "ffff"
    

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

/speakers/{ ''' % (mask, NSPK, input_scale)

print h1
outf.write(h1+"\n")

string = ""
for jj,spkname in enumerate(spkord):
    print "add_spkr    %s    %.3f    %.5f    %.5f    system:playback_%d" % (LF.label[jj], LF.ra[jj], LF.gaz[jj], LF.gel[jj], spkname)
    string = "add_spkr    %s    %.3f    %.5f    %.5f    system:playback_%d" % (LF.label[jj], LF.ra[jj], LF.gaz[jj], LF.gel[jj], spkname)
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
    for column,val in  enumerate(hf_mat[row]):
        string += "   %.6f" % val
    print "add_row " +string
    outf.write("add_row " +string+"\n")
    string = ""


h4 = '''/}

/end'''
print h4
outf.write(h4+"\n")
