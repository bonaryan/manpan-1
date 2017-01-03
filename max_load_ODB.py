import numpy as np
import os
import string
import sys
import odbAccess

from odbAccess import *
from abaqusConstants import *

# Import pickle to load the .pkl database
import pickle

# Open and read the database
profiles_file = open("profiles.pkl",'rb')
profiles = pickle.load(profiles_file)
profiles_file.close()

profiles_file = open("meta.pkl",'rb')
profiles_meta = pickle.load(profiles_file)
profiles_file.close()


NameOfFile='maxforcedispldata-N.txt'
out = open(NameOfFile,'w')



for b in (3, 4, 5):
	for l in (1, 2, 3):
		for k in (3, 5):
			for j in (2, 3, 4):
				for i in range (1,2,1):

					# Variables holding information of the current profile
					name='RIKS-N-1-'+str(j)+'-'+str(k)+'-'+str(l)+'-'+str(b)+'.odb'
					current_d = float(profiles_meta[i-1][j-1][k-1][l-1][0][0])
					current_t = float(profiles_meta[i-1][j-1][k-1][l-1][1][0])
					current_l = float(profiles_meta[i-1][j-1][k-1][l-1][7][0])
					current_fy = float(profiles_meta[i-1][j-1][k-1][l-1][3][0])
					current_area = float(profiles_meta[i-1][j-1][k-1][l-1][4][0])
					current_Iy = float(profiles_meta[i-1][j-1][k-1][l-1][5][0])
					E=210000
					epsilon=(sqrt(235/current_fy))*(sqrt(235/current_fy))
					current_lambda = current_l/(pi*sqrt(E*current_Iy/(current_area*current_fy)))
					
					nameOfStep='RIKS'
					myOdb = odbAccess.openOdb(path=name)
					RIKS= myOdb.steps[nameOfStep]
									
					rp1key = RIKS.historyRegions.keys()[1]
					ho1key = RIKS.historyRegions[rp1key].historyOutputs.keys()[0]
					rp2key = RIKS.historyRegions.keys()[2]
					ho2key = RIKS.historyRegions[rp2key].historyOutputs.keys()[0]
					load_hist = RIKS.historyRegions[rp1key].historyOutputs[ho1key].data
					disp_hist = RIKS.historyRegions[rp2key].historyOutputs[ho2key].data
					maxpos = load_hist.index(max(load_hist,key=lambda x:x[1]))
					load = load_hist[maxpos][1]
					disp = -disp_hist[maxpos][1]
					
					out.write(str(current_d)+'\t'+str(round(current_d/current_t,2))+'\t'+str(current_lambda)+'\t'+str(b)+'\t'+str(load)+'\t'+
					str(disp)+'\t'+str(current_d/(current_t*epsilon))+'\n')
					myOdb.close()
out.close()