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


import __builtin__ as bt
from ConfigParser import SafeConfigParser,ParsingError
import json
import numpy as np
import os
import sys
import warnings
import string
import math as mh ## this is used only for factorial...

from functions import ambisoniC


## defining points over the sphere
def SpherePt(NP):
    if nD == "2D":
        phiok = [2.0*np.pi*(float(i)/NP) for i in range(NP)]
        thetaok = [0.0 for i in range(NP)]

    if nD == "3D":
        thetaPrev = [(np.pi/2.0-(float(i)/NP)*np.pi) for i in range(NP)]
        theta = []
        phi = []
        
        for i, theP in enumerate(thetaPrev):
            n = max(int(2*NP*np.cos(theP)),1)
            phi.append( [(float(jj)/n*2)*np.pi for jj in range(n)] )
            temp = [theP] *(n)
            theta.append(temp)

        phiok = [item for sublist in phi for item in sublist]
        thetaok = [item for sublist in theta for item in sublist]

    if len(phiok)!=len(thetaok): sys.exit("Died generating points on the sphere")
    return np.asarray(phiok), np.asarray(thetaok)




def configConst(configfile):
    # Requre values
    try:
        parser = SafeConfigParser()
        parser.read(configfile)
    
        bt.case = parser.get('Layout','name')
        bt.PHI = json.loads(parser.get('Layout', 'PHI'))
        bt.THETA = json.loads(parser.get('Layout', 'THETA'))
    
        bt.OUTFILE = parser.get('Outfiles','matlab')
    
        bt.DEC = parser.get('Ambisonics','DEC')
        bt.DEG = parser.getint('Ambisonics','DEG')
        bt.nD = parser.get('Ambisonics','nD')
    
        bt.MUTESPKRS =parser.getboolean('Flags','MUTESPKRS')
        bt.AUTOREM =  parser.getboolean('Flags','AUTOREM') 
        bt.PREFHOM =  parser.getboolean('Flags','PREFHOM') 
        bt.WFRONT =   parser.getboolean('Flags','WFRONT')  
        bt.WPLANE =   parser.getboolean('Flags','WPLANE')  
        bt.WBIN =     parser.getboolean('Flags','WBIN') 
        bt.WAUTOREM = parser.getboolean('Flags','WAUTOREM') 
        bt.thetaThreshold = parser.getint('Flags','thetaThreshold')
        bt.MATCHSPK = parser.getboolean('Flags','MATCHSPK') 
        bt.MATCHTOL = parser.getint('Flags','MATCHTOL') 
    
        bt.CP = parser.getint('Minimiz','CP') 
        bt.CV = parser.getint('Minimiz','CV')
        bt.CE = parser.getint('Minimiz','CE')
        bt.CR = parser.getint('Minimiz','CR')
        bt.CT = parser.getint('Minimiz','CT')
        bt.CPH= parser.getint('Minimiz','CPH')
    
        bt.SEED= parser.getint('Other','SEED')
         
    
    except ParsingError, err:
        print 'Could not parse:', err
    
        
    ## number of speakers
    bt.NSPK = len(bt.PHI)
    
    bt.PHI, bt.THETA = fixPHI(bt.PHI,bt.THETA)
    
    bt.phiTest, bt.thetaTest, bt.NPOINTS, bt.WfrontVec, bt.WplaneVec, bt.WbinVec, bt.WremVec = autoinit(bt.PHI,bt.THETA,bt.SEED,bt.AUTOREM,bt.thetaThreshold)
    
    bt.coeffDir = ambisoniC(bt.phiTest,bt.thetaTest,'basic',bt.DEG,0) # calculating ambisonics coefficients in test directions
    
    
    bt.SIGNSvec = SHsigns(bt.DEG)
    bt.MATCHED, bt.PHI, bt.THETA = analyzeLayout(bt.PHI,bt.THETA,bt.MATCHTOL*np.pi/180)
    bt.MAP = MapMatched(bt.MATCHED)
    
    if bt.MATCHSPK: bt.NSPKmatch = len(bt.MATCHED[bt.MATCHED>=0])
    else: bt.NSPKmatch = bt.NSPK



#############
## weights ##
#############

def Wfront(phiT,thetaT):
    #wvect = 1. + np.cos(phiT) * np.cos(thetaT)/2.
    wvect = np.cos(phiT)**8 * np.cos(thetaT)**8 /2.
    return wvect

def Wplane(thetaT):
    wvect = 1. + np.cos(thetaT)**2.
    return wvect

def Wbinary(thetaT,tThresh):
    return  thetaT > tThresh


def angDist(phi1,e1,phi2,e2):
    import numpy as np
    # Implement haversine formula (see Wikipedia)
    dist = 2.0*np.arcsin(np.sqrt( (np.sin((e1 - e2)/2.0))**2 + np.cos(e1)*np.cos(e2)*(np.sin((phi1 - phi2)/2.0))**2 ) )
    return dist

def spkDist():
    import numpy as np
    PHI = bt.PHI
    THETA = bt.THETA
    
    dvec = np.array([])
    mins = np.array([])
    for i in range(len(PHI)):
        for j in range(len(PHI)): # you could do probably something like: for j in range(i+1,len(PHI))
            dvec= np.append(dvec,[angDist(PHI[i],THETA[i],PHI[j],THETA[j])])

        mins = np.append(mins,[min(dvec[dvec!=0])])
        dvec = np.array([]) # reset dvec
    mean = np.mean(mins) # calculates the mean only of the smalles values - the closest speakers
    return mean



def autoremoval(PHI,THETA,phiT,thetaT):

    meanSpkDist = spkDist()
    phit = []
    thetat = []
    for i in range(len(phiT)):
        for j in range(len(PHI)):
            if angDist(phiT[i],thetaT[i],PHI[j],THETA[j]) < meanSpkDist*1.0:
                phit.append(phiT[i])
                thetat.append(thetaT[i])
                break
    return phit, thetat



def Wautoremoval(PHI,THETA,phiT,thetaT):

    meanSpkDist = spkDist()
    
    wvect = []
    for i in range(len(phiT)):
        for j in range(len(PHI)):
            if angDist(phiT[i],thetaT[i],PHI[j],THETA[j]) < meanSpkDist*1.0:
                temp = True
                #temp = 1
                break
            else: temp = False
            #else: temp = 0.3
        wvect.append(temp)

    return np.asarray(wvect)



def autoinit(PHI,THETA,SEED,AUTOREM,thetaThreshold):
    ###############################
    ## automatic initializations ##
    ###############################
    # Threshold for binary masking the lower region
    thetaThreshold *= np.pi/180.0

    phiTest,thetaTest = SpherePt(SEED) ## generating sampling space points

    if AUTOREM: phiTest, thetaTest = autoremoval(PHI,THETA,phiTest,thetaTest) # autoremoving test points without speakers around

    WfrontVec = Wfront(phiTest,thetaTest)
    WplaneVec = Wplane(thetaTest)
    WbinVec = Wbinary(thetaTest,thetaThreshold)
    WremVec = Wautoremoval(PHI,THETA,phiTest,thetaTest) # only necessary if WAUTOREM is true...

    NPOINTS = len(phiTest)

    return phiTest, thetaTest, NPOINTS, WfrontVec, WplaneVec, WbinVec, WremVec



def cart2sph(x,y,z):
    """From cartesian to spherical coordinates"""
    """Maths convention!"""
    XsqPlusYsq = x**2 + y**2
    r = np.sqrt(XsqPlusYsq + z**2)               # r
    elev = np.arctan2(z,np.sqrt(XsqPlusYsq))     # theta
    az = np.arctan2(y,x)                           # phi
    # this is in Maths convention...
    # probably to convert to Acoustics convention

    #  math.atan2(y, x)
    # Returns atan(y / x), in radians. The result is between -pi and pi
    return r, elev, az



def outp(name, value, filename):
    string = name + " = " + repr(value) + "\n"
    with open(filename, "a") as myfile:
        myfile.write(string)
        myfile.close()


def analyzeLayout(PHI,THETA,tol):
    mp = len(PHI)
    paired = np.zeros(mp)
    tdx = 0

    for idx, phi in enumerate(PHI):
        if (paired[idx]>=1): continue

        for jdx in range(idx+1,mp):
            if (paired[jdx]>=1): continue
            phj = PHI[jdx]
            if (angDist(phi,THETA[idx],-phj,THETA[jdx]) < tol):
                # pair found !
                tdx += 1

                if (phi>0 and phj<0): paired[idx] = tdx; paired[jdx] = -tdx;
                elif (phi<0 and phj>0): paired[idx] = -tdx; paired[jdx] = tdx;
                else: sys.exit("Well... there is a problem in analyzeLayout function: \nthe pair is not left/right symmetric. phi = %.2f phj = %.2f" % (phi,phj) )
                
                # maybe do this after the matching
                phiabs = abs(PHI[idx]-PHI[jdx])/2
                PHI[idx] = phiabs 
                PHI[jdx] = -phiabs
                THETA[idx]=THETA[jdx]=(THETA[idx]+THETA[jdx])/2
    
    print "Paired %d loudspeakers with tolerance %.2f" % ((len(paired[abs(paired)>=1])), tol)

    if (len(paired[abs(paired)>=1])%2 != 0): sys.exit("Something went wrong when matching pairs in analyzeLayout. The number of matched pairs is: %d" % (len(paired[paired>=1])))
    return paired, PHI, THETA

def fixPHI(PHI, THETA):
    # convert in radiants if necessary
    if (np.asarray(PHI)>2.0*np.pi).any(): 
        ## conversion in radiants
        PHI = [PHI[i]*np.pi/180.0 for i in range(len(PHI))]
        THETA = [THETA[i]*np.pi/180.0 for i in range(len(THETA))]
    
    # define positive PHIs
    if (np.asarray(PHI)<0).any():
        PHI = [2*np.pi+PHI[i] if PHI[i]<0 else PHI[i] for i in range(len(PHI))]
    # then fold them again between -pi,+pi 
    if (np.asarray(PHI)>np.pi).any():
        PHI = [-np.pi+(PHI[i]%np.pi) if PHI[i]>np.pi else PHI[i] for i in range(len(PHI))]

    return PHI, THETA

def SHsigns(DEG):
    signs = []
    for deg in range(DEG+1):
        signs.append([-1 if i<deg else 1 for i in range(2*deg+1)])

    signs = [item for sublist in signs for item in sublist]

    return signs


def MapMatched(matched):
    # the idea here is to speed up the upscalematrix function
    # because this 
    # wantedrow = np.argwhere(matched[matched>=0] == abs(matched[idx]))[0][0]
    # is not very efficient
    matched = matched[matched>=0]
    mapped = np.zeros(np.shape(matched))

    for idx, val in enumerate(matched):
        if (val != 0): 
            mapped[val-1] = idx

    return mapped

def fac2(val) :
    ## stupid function to calculate factorial2 not to import another library (scipy)
    if val <= 0: return 1.0
    else: 
        tmp = val
        for ii in range(val,1,-2):
            if ii < val: tmp = tmp*ii

        return float(tmp)


################################################################
# documentation and tools to change from different conventions #
################################################################
class Conventions :
    # documentation and tools to change from different conventions #

    def ch(self) :
        # Ambisonic Channel Number by proposed AmbiX standard
        #  0 (0,0)  W
        #  1 (1,-1) Y
        #  2 (1,0)  Z
        #  3 (1,1)  X
        #  4 (2,-2) V
        #  5 (2,-1) T
        #  6 (2,0)  R
        #  7 (2,1)  S
        #  8 (2,2)  U
        #  9 (3,-3) Q
        # 10 (3,-2) O
        # 11 (3,-1) M
        # 12 (3,0)  K
        # 13 (3,1)  L
        # 14 (3,2)  N
        # 15 (3,3)  P
        # (up to third order tehre are some assigned letters, and then follows with similar pattern)
        self.ACN3D = "W Y Z X V T R S U Q O M K L N P"
        self.ACN2D = "W Y X V U Q P"



    # conversions
    # example of use: Y(SN2D) = sn2D.sn3D * Y(SN3D)

    def conv(self) :
        # look p156 Daniel's thesis
        # maxn/gerz and furse malham are essentially the same except for the 1/sqrt(2)m in W
        class N3D :
            pass
        n3d = N3D()

        class SN3D :
            pass
        sn3d = SN3D()

        class MaxN :
            pass
        maxn = MaxN()

        class FuMa :
            pass
        fuma = FuMa()

        class SN2D :
            pass
        sn2d = SN2D()

        class N2D :
            pass
        n2d = N2D()


        ### calculations
        #
        SN2DSN3D = np.array([])
        for m in range(self.deg+1):
            mcoeff = np.sqrt( (2.0**(2*m-1) * np.math.factorial(m)**2 ) / (1.0*np.math.factorial(2*m)) )
            tmp = [mcoeff for n in range(m**2,(m+1)**2)]
            SN2DSN3D = np.append(SN2DSN3D,tmp,axis=0)
        
        if SN2DSN3D.any() == 0.0: raise ValueError("Conversion between SN3D and N3D failed.")

        #
        N3DSN3D = np.array([])
        for m in range(self.deg+1):
            mcoeff = np.sqrt(2*m+1)
            tmp = [mcoeff for n in range(m**2,(m+1)**2)]
            N3DSN3D = np.append(N3DSN3D,tmp,axis=0)

        if N3DSN3D.any() == 0.0: raise ValueError("Conversion between SN3D and N3D failed.")
        
        # this is for l=m or l=-m ONLY!!!
        N2DN3D = np.array([])
        for m in range(self.deg+1):
            mcoeff = (mh.factorial(m) * 2**(m-1)) / fac2(2*m+1) # c**(-2) it's the 1/N3D coeff
            mcoeff = np.sqrt(2.) * np.sqrt(mcoeff)
            tmp = [mcoeff for n in range(m**2,(m+1)**2)]
            N2DN3D = np.append(N2DN3D,tmp,axis=0)
        
        if N2DN3D.any() == 0.0: raise ValueError("Conversion between SN3D and N3D failed.")
 

        #
        #MAXNtoN3D = np.array([np.sqrt(2), np.sqrt(3), np.sqrt(3), np.sqrt(3), np.sqrt(15)/2, np.sqrt(15)/2, np.sqrt(15), np.sqrt(15)/2, np.sqrt(15)/2, np.sqrt(35./8.), np.sqrt(35)/3, np.sqrt(224./45.), np.sqrt(7), np.sqrt(224./45.), np.sqrt(35)/3, np.sqrt(35./8.)])
        MAXNtoN3D = np.array([1.0, np.sqrt(3), np.sqrt(3), np.sqrt(3), np.sqrt(15)/2, np.sqrt(15)/2, np.sqrt(5), np.sqrt(15)/2, np.sqrt(15)/2, np.sqrt(35./8.), np.sqrt(35)/3, np.sqrt(224./45.), np.sqrt(7), np.sqrt(224./45.), np.sqrt(35)/3, np.sqrt(35./8.)])
        ## http://en.wikipedia.org/wiki/Ambisonic_data_exchange_formats
        MAXNtoN3D = np.resize( MAXNtoN3D,(1,(self.deg+1)**2) ).flatten()
        if len(MAXNtoN3D) != (self.deg+1)**2: raise ValueError("Conversion between MaxN and N3D is badly implemented: %s vs %s elements." % (len(MAXNtoN3D),(self.deg+1)**2) )

        FUMAtoN3D = MAXNtoN3D*1.0
        FUMAtoN3D[0] = np.sqrt(2)

        sn2d.sn3d = SN2DSN3D        # not completely correct since 2D exists only for n= m ... TDB
        sn3d.sn2d = 1.0/SN2DSN3D

        n3d.sn3d = N3DSN3D
        sn3d.n3d = 1.0/N3DSN3D

        maxn.n3d = MAXNtoN3D
        n3d.maxn = 1.0/MAXNtoN3D

        fuma.n3d = FUMAtoN3D
        n3d.fuma = 1.0/FUMAtoN3D    # don't touch it's correct with Daniel's thesis convention p156


        ### derived from these
        n2d.n3d = N2DN3D
        n3d.n2d = 1.0/N2DN3D
        #n2d.sn3d sn3d.n3d
        n2d.fuma = n3d.fuma / n2d.n3d
        fuma.n2d = 1.0/n2d.fuma
        #n2d.n3d n3d.fuma


        ### filling
        self.n3d  = n3d 
        self.sn3d = sn3d 
        self.maxn = maxn 
        self.fuma = fuma 
        self.sn2d = sn2d 
        self.n2d  = n2d 


    def __init__(self,deg) :
        self.deg = deg
        self.ch()
        self.conv()

        if self.deg>3 : warnings.warn("FuMa and MaxN are implementend only up to third order.")

        # I just want to warn, not to raise an exception...
        # like https://docs.python.org/2/library/warnings.html
        # or http://stackoverflow.com/questions/3891804/how-to-raise-a-warning-in-python-without-stopping-interrupting-the-program


    def shrink(self, vector) :
        ## reduction of the 2d ones from (deg+1)**2 to 2*deg+1
        if len(vector) != (self.deg+1)**2 : raise ValueError("Can't shrink this! Wrong dimension.")

        shr = np.asarray([])
        for ii in range(self.deg+1):
            aa = ii**2+1   -1
            bb = (ii+1)**2 -1
            shr = np.append(shr,vector[aa]) 
            if aa != bb: shr = np.append(shr,vector[bb])

        return shr





class Parser :

    # function for parsing the data
    def _data_parser(self, text, dic, area):
        ## ignores all the lines that do not contain as firts word the one carried by "area"
        ## all the others are processed following the rules in dictionary (dic) 
        ## the text is then splitted into words and returned as an array
        tmp = text.split()
        if len(tmp)>0 and tmp[0] == area :
            for i, j in dic.iteritems():
                text = text.replace(i,j)
            text = text.split()
            return text
        else:
            return None
    

    def ambdec(self, my_text):
        ## method that does the actual processing 
        dic = {'add_spkr':' ','add_row':' ','^I':' '}
        
        ## processing the first part of ambdec config file
        ## where the information on the layout is stored
        tmp = np.array(["spk_id","dist","az","el","connect"])
        for line in my_text:
            layout = self._data_parser(text=line, dic=dic, area="add_spkr")
            if layout != None: tmp = np.vstack((tmp, layout))
        
        self.layout = tmp
        self.layout_data = tmp[1:,1:4].astype(np.float)
        
        ## processing the second part of the ambdec config file
        ## where the decoding coefficients are stored
        tmp = np.array([])
        count = 0
        for line in my_text:
            parsed = self._data_parser(text=line, dic=dic, area="add_row")
            if parsed != None:
                if count == 0: tmp = parsed; count += 1
                else:          tmp = np.vstack((tmp,parsed))
        
        ## splitting into the two matrices: lf and hf
        self.lf_mat = tmp[:self.layout_data.shape[0]].astype(np.float)
        self.hf_mat = tmp[self.layout_data.shape[0]:].astype(np.float)
        if self.lf_mat.shape != self.hf_mat.shape : raise ValueError("Something went wrong while parsing ambdec file.")


    def __init__(self, filename):
        self.inputfile = open(filename)
        my_text = self.inputfile.readlines()
        self.ambdec(my_text)

#### Auxiliary funcitons for evaluation

def analyze_results(vect, varname, dB = None) :
    mean = np.mean(vect)
    std = np.std(vect)
    print "%s mean: %.2f; std: %.2f"%(varname, mean, std)
    
    if dB is None:
        return (mean,std)
    
    elif dB == 10 or dB == 20:
        mean_dB = dB * np.log10(mean)
        std_dB = dB * np.log10(std)
        print "%s mean (dB): %.1f; std (dB): %.1f"%(varname, mean_dB, std_dB)        
        return (mean,std,mean_dB,std_dB)

    else:
        raise ValueError("dB should either be 10, 20 or None")

