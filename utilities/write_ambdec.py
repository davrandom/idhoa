from auxiliary import cart2sph, fixPHI, Conventions
from plotting import SpeakersPlotting
import numpy as np
from ConfigParser import SafeConfigParser,ParsingError
import json


# 5.0 for paper layout
configfile = "./init_files/fivedotzero.ini"
try:
    parser = SafeConfigParser()
    parser.read(configfile)
    
    case = parser.get('Layout','name')
    az = json.loads(parser.get('Layout', 'PHI'))
    el = json.loads(parser.get('Layout', 'THETA'))
    ra = json.loads(parser.get('Layout', 'radius'))
    label = json.loads(parser.get('Layout', 'chlabel'))
    
    DEC = parser.get('Ambisonics','DEC')
    DEG = parser.getint('Ambisonics','DEG')
  
    case = parser.get('Layout','name')
    AUTOREM =  parser.getboolean('Flags','AUTOREM') 
    MATCHSPK = parser.getboolean('Flags','MATCHSPK') 
    CP = parser.getint('Minimiz','CP') 
    CV = parser.getint('Minimiz','CV')
    CE = parser.getint('Minimiz','CE')
    CR = parser.getint('Minimiz','CR')
    CT = parser.getint('Minimiz','CT')
    CPH= parser.getint('Minimiz','CPH')

except ParsingError, err:
    print 'Could not parse:', err
 
spkord = [4, 1, 3, 2, 5] ## this is important to get the connections right in jack


az, el = fixPHI(az,el)
print "\nradius (cm), elevation, azimut (radiants)"
print ra
print el
print az

# from rad to grad
gaz = [180.0*ii/np.pi for ii in az]
gel = [180.0*ii/np.pi for ii in el]

print "\nradius (cm), elevation, azimut (grads)"
print ra
print gel
print gaz

#SpeakersPlotting(az,el,ra)
#SpeakersPlotting(az,el,1)

######################################################################################
######################################################################################
######################################################################################

from numpy import array

pyfilename = case+"-"+str(DEG)+"-"+str(DEC)+"-rem"+str(AUTOREM)+"-sym"+str(MATCHSPK)

if (DEC=="basic"): pyfilename += "CP"+str(CP)+"CV"+str(CV)+".py" 
if (DEC=="maxRe"): pyfilename += "CR"+str(CR)+"CT"+str(CT)+"CE"+str(CE)+".py"
if (DEC=="phase"): pyfilename += "CR"+str(CR)+"CT"+str(CT)+"CPH"+str(CPH)+".py"

#execfile("torrente-fz-2-maxRe-remFalse-symFalseCR100CT1400CE400.py")
execfile(pyfilename)

horiz_deg = int((ResCoeff.T.shape[1]-1)/2.)

conv = Conventions(horiz_deg)
lf_mat = GuessPinv.T#*conv.shrink(conv.fuma.n2d)
hf_mat = ResCoeff.T#*conv.shrink(conv.fuma.n2d)


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
/dec/coeff_scale  fuma

/opt/input_scale  fuma
/opt/nfeff_comp   input
/opt/delay_comp   off
/opt/level_comp   off
/opt/xover_freq    400
/opt/xover_ratio   0.0

/speakers/{ ''' % (mask, NSPK)

print h1


for jj,spkname in enumerate(spkord):
    print "add_spkr    %s    %.3f    %.5f    %.5f    system:playback_%d" % (label[jj], ra[jj]/100, gaz[jj], gel[jj], spkname)

h2 = '''/}

/lfmatrix/{
order_gain     1.00000  1.00000  1.00000  1.00000'''
print h2

string = ""
for row in range(len(lf_mat)):
    for column,val in  enumerate(lf_mat[row]):
        string += "   %.6f" % val
    print "add_row " +string
    string = ""

h3 = '''/}

/hfmatrix/{
order_gain     1.00000  1.00000  1.00000  1.00000'''
print h3
for row in range(len(hf_mat)):
    for column,val in  enumerate(hf_mat[row]):
        string += "   %.6f" % val
    print "add_row " +string
    string = ""


h4 = '''/}

/end'''
print h4
