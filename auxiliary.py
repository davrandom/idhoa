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
    
        
        #######################
        ## Speakers position ##
        #######################
        case = 'studio'
    
        if (case=="studio"):
            # studio
            PHI = [27.0415460692, 2.04283140615, 332.958453931, 1.03372198009, 90.5036008348, 269.496399165, 166.80995235 , 193.19004765 , 54.0379157019, 302.679984191, 124.309021931, 235.690978069, 92.4265791571, 267.573420843, 35.1695965277, 324.830403472, 15.3392440494, 344.660755951, 30.5463301933, 329.453669807, 143.681356818, 216.318643182, 301.738222232]          # position of the speakers # AZIMUTH
            THETA = [-1.35939557882, -1.39079054354, -1.35939557882, 25.2255330043 , -2.05829450955, 2.05829450955 , 6.89524790513 , 6.89524790513 , -4.49039914903, -1.34181305608, 2.73330540595 , 2.73330540595 , 27.4212399195 , 27.4212399195 , 21.065228491  , 21.065228491  , -22.9266013259, -22.9266013259, 30.6450023912 , 30.6450023912 , 27.2300782225 , 27.2300782225 , 88.0711642016]        # position of the speakers # ELEVATION
        
        elif(case=="CCRMA"):
            # CCRMA Dome
            PHI   = [22.500000, 67.500000, 112.500000, 157.500000, -157.500000, -112.500000, -67.500000, -22.500000, 30.000000, 90.000000, 150.000000, -150.000000, -90.000000, -30.000000, 0.000000];
            THETA = [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 40.000000, 40.000000, 40.000000, 40.000000, 40.000000, 40.000000, 90.000000];
        
        elif(case=="icos"):
            # icosahedron
            PHI = [45.0, -45.0,-135.0, 135.0, 45.0, -45.0,-135.0, 135.0, 90.0,  90.0, -90.0, -90.0, 0.0,   0.0, 180.0, 180.0, 69.1, -69.1, 110.9,-110.9]
            THETA = [35.3, 35.3, 35.3, 35.3,-35.3,-35.3,-35.3,-35.3, 69.1,-69.1, 69.1,-69.1, 20.9,-20.9, 20.9,-20.9, 0.0, 0.0,  0.0, 0.0]
        
        elif(case=="CARPA"):
            # CARPA
            X = [0., 0., 0., -163, -163,-163, 163, 163, 163, -485, -485, -485, 485, 485, 485]
            Y = [240, 86, -218, 240, 86, -218, 240, 86, -218, 240, -116, -218, 240., -116, -227]
            Z = [-3, 348, 281, -3, 348, 281, -3, 348, 281, -24, 213, 7, -24, 213, 7]
            if len(X) != len(Y): print "Something wrong in CARPA"
            if len(X) != len(Z): print "Something wrong in CARPA"
        
            pts = [(X[a],Y[a],Z[a]) for a in range(len(X))]
            relp = np.array([aux.cart2sph(x,y,z) for x,y,z in pts])
        
            R =  [relp[i][0] for i in range(len(relp))]
            el =  [relp[i][1] for i in range(len(relp))]
            az =  [relp[i][2] for i in range(len(relp))]
        
            PHI = az
            THETA = el
        
        elif(case=="Karlsruhe"):
            #Karlsruhe setup 
            PHI = [-0.2381, 0.4675, -0.5145, -0.1827, 0.2954, -0.3345, 0.163, 0.6278, 0.8405, 1.0381, 1.2482, 1.3916, 1.6292, 1.6912, 1.9513, 2.0034, 2.1096, 2.3073, 2.5128, 2.6254, 2.8467, 3.088, -3.0899, -2.9089, -2.6912, -2.6518, -2.5348, -2.2267, -2.225, -2.007, -1.8741, -1.7502, 0, 0.1626, 0.4083, 0.7239, 1.2337, 1.1019, 1.4266, 1.7326, 1.9609, 2.172, 2.626, 2.7546, -2.9764, 3.1416, 2.9442, -2.4255, -2.3719, -1.9198, -1.5708, -1.4515, -1.036, -0.6246, -0.5373, -0.0837, -1.5126, -1.4486, -1.2073, -1.1685, -0.8848, -0.8347, -0.8199]
            THETA = [0.1174, 0.1079, -0.2007, -0.2259, -0.22, -0.6564, -0.6413, -0.6439, -0.1307, -0.6366, -0.3503, 0.0475, -0.6406, -0.1268, -0.2497, 0.1671, -0.6385, 0.0049, -0.593, -0.1309, -0.5175, -0.1308, -0.6392, 0.0767, -0.261, -0.6411, 0.029, -0.6084, -0.1019, 0.175, -0.6413, -0.1355, 1.4449, 0.2858, 0.9461, 0.4669, 0.9792, 0.1745, 0.5625, 0.2838, 0.9628, 0.5267, 0.3842, 0.7852, 0.2128, 1.1425, 0.529, 0.3259, 0.7764, 0.5088, 1.0263, 0.3971, 0.5015, 0.9397, 0.4339, 0.6627, -0.6391, -0.039, -0.4489, 0.1024, -0.6441, -0.2665, 0.1452]
        
        elif(case=="trirectangle"):
            #tri-rectangle layout
            PHI = [30.0, 150.0, 210.0, 330.0, 0.0, 0.0, 90.0, 90.0, 180.0, 180.0, 270.0, 270.0]
            THETA = [0.0, 0.0, 0.0, 0.0, 30.0, -30.,30.0, -30.,30.0, -30.,30.0, -30.]
        
        ################
        ## Parameters ##
        ################
        OUTFILE = '~/Documents/MATLAB/AmbiDecoderToolbox/Karlsruhe_DOME.mat'
        DEC = 'maxRe' # available: basic, maxRe, phase
        DEG = 3
        
        
        # parameters
        SEED = 17   # number of maximum horizontal points in the sampling function (divided by two)
                    # from 14 to 20 it's a reasonable number
        
        # FLAGS
        MUTESPKRS = 1 # tries to mute "unnecessary" speakers
        AUTOREM = 0 # removes some points from the sphere sampling (the farthest from the speakers)
        PREFHOM = 0
        # weights
        WFRONT = 0
        WPLANE = 0
        WBIN = 1
        WAUTOREM = 0 # the same as AUTOREM but with weights
        thetaThreshold = -10 # degrees
        MATCHSPK = 1
        MATCHTOL = 5
     
        ###################################
        ## Internals of the minimization ##
        ###################################
        CP = 400. # arbitrary coefficients
        CV = 25.
        CE = 400.
        CR = 100.
        CT = 100.
        CPH= 1000.
      
    
    
    
    
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
    wvect = 1. + np.cos(phiT) * np.cos(thetaT)/2.
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


