import numpy as np
import os
import string
import sys
import odbAccess

from odbAccess import *
from abaqusConstants import *

from itertools import product

out_filename = 'maxload-summary.txt'
out = open(out_filename,'w')


prj_name = 'GNIA'
combinations = list(
    product(
        range(6, 31, 2), 
        range(27, 55, 3)
	    )
    )

# Loop through the combinations of the given input values
for model_params in combinations:
    IDstring = (''.join("%03d-"%e for e in model_params) + prj_name)
    IDstring = IDstring.translate(None, '.')
    myOdb = odbAccess.openOdb(path=IDstring+'/RIKS'+IDstring+'.odb')
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
    out.write(str(model_params)+'    '+str(lpf)+'    '+str(load)+'    '+str(disp)+'    '+'\n')
    odbAccess.closeOdb(myOdb)

out.close()
