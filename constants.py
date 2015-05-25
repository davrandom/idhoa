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


import __builtin__ as bt
from ConfigParser import SafeConfigParser,ParsingError
import json
import numpy as np
import os
import sys
import auxiliary as aux

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
    
    bt.PHI, bt.THETA = aux.fixPHI(bt.PHI,bt.THETA)
    
    bt.phiTest, bt.thetaTest, bt.NPOINTS, bt.WfrontVec, bt.WplaneVec, bt.WbinVec, bt.WremVec = aux.autoinit(bt.PHI,bt.THETA,bt.SEED,bt.AUTOREM,bt.thetaThreshold)
    
    
    
    bt.SIGNSvec = aux.SHsigns(bt.DEG)
    bt.MATCHED, bt.PHI, bt.THETA = aux.analyzeLayout(bt.PHI,bt.THETA,bt.MATCHTOL*np.pi/180)
    bt.MAP = aux.MapMatched(bt.MATCHED)
    
    if MATCHSPK: bt.NSPKmatch = len(bt.MATCHED[bt.MATCHED>=0])
    else: bt.NSPKmatch = bt.NSPK


