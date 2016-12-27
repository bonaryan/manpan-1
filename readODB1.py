import numpy as np
import os
import string
import sys
import odbAccess

from odbAccess import *
from abaqusConstants import *

NameOfFile='maxforcedispldata.dat'
out = open(NameOfFile,'w')


for b in (3, 4, 5):
	for l in (1, 2, 3):
		for k in (3, 5):
			for j in (2, 3, 4):

				# Variables holding information of the current profile
				name='RIKS-N-1-'+str(j)+'-'+str(k)+'-'+str(l)+'-'+str(b)+'.odb'
                
				nameOfStep='RIKS'
				myOdb = odbAccess.openOdb(path=name)
				RIKS= myOdb.steps[nameOfStep]
                				
				rp1key = RIKS.historyRegions.keys()[1]
				ho1key = RIKS.historyRegions[rp1key].historyOutputs.keys()[0]
				rp2key = RIKS.historyRegions.keys()[2]
				ho2key = RIKS.historyRegions[rp2key].historyOutputs.keys()[0]
				load_hist = RIKS.historyRegions[rp1key].historyOutputs[ho1key].data
				disp_hist = RIKS.historyRegions[rp2key].historyOutputs[ho2key].data
				maxpos = load_hist.index(max(load_hist, key=lambda x: x[1]))
				load = load_hist[maxpos][1]
				disp = -disp_hist[maxpos][1]
				
				out.write(str(load)+'\t'+str(disp)+'\n')
				myOdb.close()
out.close()