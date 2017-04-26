import numpy as np
import os
import string
import sys
import odbAccess

from odbAccess import *
from abaqusConstants import *

NameOfFile='maxload.txt'
out = open(NameOfFile,'w')

for j_sides in range(5, 18):
    sides_string = "%02d" % (j_sides,)
    for i_classification in range(30, 49, 3):
        myOdb = odbAccess.openOdb(path=sides_string+'/'+str(j_sides)+'-'+str(i_classification)+'/'+sides_string+'-'+str(i_classification)+'.odb')
        RIKSstep = myOdb.steps['RIKS']
        rp1key = RIKSstep.historyRegions.keys()[1]
        ho1key = RIKSstep.historyRegions[rp1key].historyOutputs.keys()[0]
        rp2key = RIKSstep.historyRegions.keys()[2]
        ho2key = RIKSstep.historyRegions[rp2key].historyOutputs.keys()[0]
        asskey = RIKSstep.historyRegions.keys()[0]
        hoasse = RIKSstep.historyRegions[asskey].historyOutputs.keys()[-1]
        load_hist = RIKSstep.historyRegions[rp1key].historyOutputs[ho1key].data
        disp_hist = RIKSstep.historyRegions[rp2key].historyOutputs[ho2key].data
        lpf_hist = RIKSstep.historyRegions[asskey].historyOutputs[hoasse].data
        maxpos = load_hist.index(max(load_hist,key=lambda x:x[1]))
        load = load_hist[maxpos][1]
        disp = -disp_hist[maxpos][1]
        lpf = lpf_hist[maxpos][1]
        out.write(str(j_sides)+'-'+str(i_classification)+'    '+str(lpf)+'    '+str(load)+'    '+str(disp)+'\n')


out.close()
