import numpy as np
import os
import string
import sys
import odbAccess

from odbAccess import *
from abaqusConstants import *




for i in range(1):
	for j in range(1,3,1):
		for k in range(2,5,2):
			for l in range(0,3,1):
				for b in range(2,5,1):

					# Variables holding information of the current profile
					names=['RIKS-N-'+str(i+1)+'-'+str(j+1)+'-'+str(k+1)+'-'+str(l+1)+'-'+str(b+1)]


					for y in range(len(names)):
						nameOfStep='RIKS'
						NameOfFile='maxforcedispldata-'+'.txt'
						FileResultsX=open(NameOfFile,'w')
						Name=names[y]+'.odb'
						myOdb = odbAccess.openOdb(path=Name)
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
						
						
					
						for x in range (len(Name)):
							FileResultsX.write('%10.8E\t %10.8E\t\n' % (load, disp))
							
					FileResultsX.write('\n')

	myOdb.close()	
	FileResultsX.close()
	
	