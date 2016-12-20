import numpy as np
import string
import sys
import os
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
#from visualization import *
from connectorBehavior import *
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
# Import profiles database and profile metadata --------------------------------------------------------------

# Import pickle to load the .pkl database
import pickle

# Define a method to get the block number of a specific string in the keywords

def GetBlockPosition(model,blockPrefix):
	pos = 0
	for block in model.keywordBlock.sieBlocks:
		if string.lower(block[0:len(blockPrefix)])==string.lower(blockPrefix):
			return pos
		pos=pos+1
	return -1

# Open and read the database
profiles_file = open("./profiles.pkl",'rb')
profiles = pickle.load(profiles_file)
profiles_file.close()

profiles_file = open("./meta.pkl",'rb')
profiles_meta = pickle.load(profiles_file)
profiles_file.close()

# number of corners
#i = int(sys.argv[-5])
i = 0

# diameter of the profile
#j = int(sys.argv[-4])
j = 0

# Profile slenderness
#k = int(sys.argv[-3])
k = 0

# Member slenderness
#l = int(sys.argv[-2])
l = 0

# bolt spacing to diameter ratio (s/d)
#b = int(sys.argv[-1])
b = 6

# Variables holding information of the current profile
buckle_model = 'BCKL-'+str(i+1)+'-'+str(j+1)+'-'+str(k+1)+'-'+str(l+1)
current_d = float(profiles_meta[i][j][k][l][0][0])
current_t = float(profiles_meta[i][j][k][l][1][0])
current_tg = float(profiles_meta[i][j][k][l][2][0])
current_fy = float(profiles_meta[i][j][k][l][3][0])
current_l = float(profiles_meta[i][j][k][l][7][0])
current_llip = sqrt((profiles[i][j][k][0][0]-profiles[i][j][k][0][2])**2+(profiles[i][j][k][1][0]-profiles[i][j][k][1][2])**2)
area = profiles_meta[i][j][k][l][4][0]
current_Iy = float(profiles_meta[i][j][k][l][5][0])

# Buckling model ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

mdb.Model(modelType=STANDARD_EXPLICIT, name=buckle_model)
c_model = mdb.models[buckle_model]


# Delete initial model
del mdb.models['Model-1']

# Create Parts ----------------------------------------------------------------------------------------------------------------------

# Sector

# -Profile sketch for sector
sector_sketch = c_model.ConstrainedSketch(name='sector', sheetSize=1200.0)

# -Sketch sector lines
for n in range(profiles[i][j][k].shape[1]-1):
	sector_sketch.Line(
		point1=(profiles[i][j][k][0][n], profiles[i][j][k][1][n]), 
		point2=(profiles[i][j][k][0][n+1], profiles[i][j][k][1][n+1])
		)

# -Extrude sector part
l_tot = 2*current_l + 3*current_d
sector_part = c_model.Part(
	dimensionality=THREE_D,
	name='sector',
	type=DEFORMABLE_BODY
	)
sector_part.BaseShellExtrude(
	depth=l_tot,
	sketch=sector_sketch
	)

# Calculate bolt positions

# -Distance on the width
bolts_w = current_llip/2

# -Distances on the length
current_b =  b
s = current_b*current_d
(n0, s0) = divmod(current_l, s)
s1 = (s0 + s)/2

bolts_z1 = np.concatenate([[bolts_w], bolts_w + ((current_d - current_llip)/5) * np.linspace(1, 4, 4), [current_d - bolts_w]])
bolts_z2 = np.concatenate([[current_d + s1], s1 + current_d + (s * np.linspace(1, n0-1, n0-1))])
bolts_z3 = bolts_z1 + (current_l + current_d)
bolts_z4 = bolts_z2 + (current_l + current_d)
bolts_z5 = bolts_z3 + (current_l + current_d)

bolts_z = np.concatenate([bolts_z1, bolts_z2, bolts_z3, bolts_z4, bolts_z5])

# Washer diameter
d_washer = 30

# Initiate list to store datum planes
datum_p=[]
					
# Make holes
					
for o in range(int(bolts_z.shape[0])):
							
	sector_part.HoleBlindFromEdges(
		depth=1.0,
		diameter=d_washer,
		distance1=bolts_z[o],
		distance2=bolts_w,
		edge1=sector_part.edges.getClosest(coordinates=((profiles[i][j][k][0][1], profiles[i][j][k][1][1], 0),))[0][0],
		edge2=sector_part.edges.getClosest(coordinates=((profiles[i][j][k][0][0], profiles[i][j][k][1][0], 1),))[0][0],
		plane=sector_part.faces.getClosest(coordinates=((profiles[i][j][k][0][0], profiles[i][j][k][1][0], 0),))[0][0],
		planeSide=SIDE1
		)
	
	sector_part.HoleBlindFromEdges(
		depth=1.0,
		diameter=d_washer,
		distance1=bolts_z[o],
		distance2=bolts_w,
		edge1=sector_part.edges.getClosest(coordinates=((profiles[i][j][k][0][-2], profiles[i][j][k][1][-2], 0),))[0][0],
		edge2=sector_part.edges.getClosest(coordinates=((profiles[i][j][k][0][-1], profiles[i][j][k][1][-1], 1),))[0][0],
		plane=sector_part.faces.getClosest(coordinates=((profiles[i][j][k][0][-1], profiles[i][j][k][1][-1], 0),))[0][0],
		planeSide=SIDE1
		)
	
	# Create datum planes to be used for partitioning the sector
	
	datum1=sector_part.DatumPlaneByPrincipalPlane(
		offset=bolts_z[o]-bolts_w, 
		principalPlane=XYPLANE
		)
	datum2=sector_part.DatumPlaneByPrincipalPlane(
		offset=bolts_z[o]+bolts_w, 
		principalPlane=XYPLANE
		)
	datum_p.append(datum1)
	datum_p.append(datum2)

# Partition the sector

# -Number of datum planes
n_dat = int(len(sector_part.datums))

# cut all the faces using the datum planes
					
for o in range((n_dat-2)):
#for o in range(2):
	sector_part.PartitionFaceByDatumPlane(
		datumPlane=sector_part.datums.items()[o+1][1],
		faces=sector_part.faces[:]
		)

# Gusset

# -Profile sketch for gusset
gusset_sketch=c_model.ConstrainedSketch(name='__profile__', sheetSize=1200.0)

# -Sketch gusset lines
# First point of the first sector
x0 = profiles[i][j][k][0][0]
y0 = profiles[i][j][k][1][0]

# Angle of the first gusset fin
phi = pi*5/6

# Calculate the end point of the gusset's first fin as an orthogonal projection of the sector's first point on the line of the gusset plate
gp1 = np.array([(x0*cos(phi)+y0*sin(phi))*cos(phi), (x0*cos(phi)+y0*sin(phi))*sin(phi)])

# Rotation matrix to multiply the previous point in order to get the points of the other 2 gusset fins
Rmat = np.array([[cos(-2*pi/3), -sin(-2*pi/3)], [sin(-2*pi/3), cos(-2*pi/3)]])

# Calculate the end points of the other 2 gusset fins by multiplying with the 120 degrees rotation matrix
gp2 = gp1.dot(Rmat)
gp3 = gp2.dot(Rmat)

# Draw lines for the sketch of the gusset plate between 0, 0 and the calculated points gp1, gp2, gp3
gusset_sketch.Line(
	point1=(0.0, 0.0),
	point2=(gp1[0], gp1[1])
	)
gusset_sketch.Line(
	point1=(0.0, 0.0), 
	point2=(gp2[0], gp2[1])
	)
gusset_sketch.Line(
	point1=(0.0, 0.0), 
	point2=(gp3[0], gp3[1])
	)

# -Extrude gusset part
gusset_part=c_model.Part(
	dimensionality=THREE_D,
	name='gusset',
	type=DEFORMABLE_BODY
	)
gusset_part.BaseShellExtrude(
	depth=current_d,
	sketch=gusset_sketch
	)

# -Holes
for o in range(int(bolts_z1.shape[0])):
	gusset_part.HoleBlindFromEdges(
		depth=1.0,
		diameter=d_washer,
		distance1=bolts_z1[o],
		distance2=bolts_w,
		edge1=gusset_part.edges.getClosest(coordinates=((gp1[0]/2, gp1[1]/2, 0),))[0][0],
		edge2=gusset_part.edges.getClosest(coordinates=((gp1[0], gp1[1], 1),))[0][0], 
		plane=gusset_part.faces.getClosest(coordinates=((gp1[0], gp1[1], 0),))[0][0], 
		planeSide=SIDE1
		)
	
	gusset_part.HoleBlindFromEdges(
		depth=1.0,
		diameter=d_washer,
		distance1=bolts_z1[o],
		distance2=bolts_w,
		edge1=gusset_part.edges.getClosest(coordinates=((gp2[0]/2, gp2[1]/2, 0),))[0][0],
		edge2=gusset_part.edges.getClosest(coordinates=((gp2[0], gp2[1], 1),))[0][0],
		plane=gusset_part.faces.getClosest(coordinates=((gp2[0], gp2[1], 0),))[0][0],
		planeSide=SIDE1
		)
	
	gusset_part.HoleBlindFromEdges(
		depth=1.0,
		diameter=d_washer,
		distance1=bolts_z1[o],
		distance2=bolts_w,
		edge1=gusset_part.edges.getClosest(coordinates=((gp3[0]/2, gp3[1]/2, 0),))[0][0],
		edge2=gusset_part.edges.getClosest(coordinates=((gp3[0], gp3[1], 1),))[0][0],
		plane=gusset_part.faces.getClosest(coordinates=((gp3[0], gp3[1], 0),))[0][0],
		planeSide=SIDE1
		)

# Partition gusset

gusset_part.DatumPointByCoordinate((gp1[0]-current_llip*cos(5*pi/6), gp1[1]-current_llip*sin(5*pi/6), 0),)
gusset_part.DatumPointByCoordinate((gp1[0]-current_llip*cos(5*pi/6), gp1[1]-current_llip*sin(5*pi/6), current_d),)

gusset_part.PartitionFaceByShortestPath(
	faces=gusset_part.faces.getClosest(coordinates=((gp1[0], gp1[1], 0),))[0][0],
	point1=gusset_part.datum.items()[0][1],
	point2=gusset_part.datum.items()[1][1],
	)

gusset_part.DatumPointByCoordinate((gp2[0]-current_llip*cos(-pi/2), gp2[1]-current_llip*sin(-pi/2), 0),)
gusset_part.DatumPointByCoordinate((gp2[0]-current_llip*cos(-pi/2), gp2[1]-current_llip*sin(-pi/2), current_d),)

gusset_part.PartitionFaceByShortestPath(
	faces=gusset_part.faces.getClosest(coordinates=((gp2[0], gp2[1], 0),))[0][0],
	point1=gusset_part.datum.items()[2][1],
	point2=gusset_part.datum.items()[3][1],
	)

gusset_part.DatumPointByCoordinate((gp3[0]-current_llip*cos(pi/6), gp3[1]-current_llip*sin(pi/6), 0),)
gusset_part.DatumPointByCoordinate((gp3[0]-current_llip*cos(pi/6), gp3[1]-current_llip*sin(pi/6), current_d),)

gusset_part.PartitionFaceByShortestPath(
	faces=gusset_part.faces.getClosest(coordinates=((gp3[0], gp3[1], 0),))[0][0],
	point1=gusset_part.datum.items()[4][1],
	point2=gusset_part.datum.items()[5][1],
	)

# Material ----------------------------------------------------------------------------------------------------------------------

c_model.Material(name='optim355')
c_model.materials['optim355'].Elastic(table=((210000.0, 0.3), ))

# Create sections ---------------------------------------------------------------------------------------------------------------

# -for sector
c_model.HomogeneousShellSection(
	idealization=NO_IDEALIZATION, 
	integrationRule=SIMPSON,
	material='optim355',
	name='sector',
	numIntPts=5,
	poissonDefinition=DEFAULT,
	preIntegrate=OFF,
	temperature=GRADIENT, 
	thickness=current_t,
	thicknessField='',
	thicknessModulus=None,
	thicknessType=UNIFORM,
	useDensity=OFF
	)
	
# -for gusset
c_model.HomogeneousShellSection(
	idealization=NO_IDEALIZATION, 
	integrationRule=SIMPSON,
	material='optim355',
	name='gusset',
	numIntPts=5,
	poissonDefinition=DEFAULT,
	preIntegrate=OFF,
	temperature=GRADIENT, 
	thickness=current_tg,
	thicknessField='',
	thicknessModulus=None,
	thicknessType=UNIFORM,
	useDensity=OFF
	)

# Assign sections ---------------------------------------------------------------------------------------------------------------

# -for sector
sector_part.Set(
	faces=sector_part.faces[:],
	name='AllSectorFaces'
	)
sector_part.SectionAssignment(
	offset=0.0,
	offsetField='',
	offsetType=MIDDLE_SURFACE,
	region=sector_part.sets['AllSectorFaces'],
	sectionName='sector', 
	thicknessAssignment=FROM_SECTION
	)

# -for gusset
gusset_part.Set(
	faces=gusset_part.faces[:],
	name='AllGussetFaces')
gusset_part.SectionAssignment(
	offset=0.0,
	offsetField='',
	offsetType=MIDDLE_SURFACE,
	region=gusset_part.sets['AllGussetFaces'],
	sectionName='gusset', 
	thicknessAssignment=FROM_SECTION
	)

# Meshing -----------------------------------------------------------------------------------------------------------------------

# Global seeding in mm
seedsize = 30

# -Sector
sector_part.setMeshControls(
	algorithm=MEDIAL_AXIS,
	elemShape=QUAD, 
	regions=sector_part.faces[:]
	)
sector_part.seedPart(
	deviationFactor=0.1, 
	minSizeFactor=0.1,
	size=seedsize
	)
sector_part.generateMesh()

# -Gusset
gusset_part.setMeshControls(
	algorithm=MEDIAL_AXIS,
	elemShape=QUAD,
	regions=gusset_part.faces[:]
	)
gusset_part.seedPart(
	deviationFactor=0.1, 
	minSizeFactor=0.1,
	size=seedsize
	)
gusset_part.generateMesh()

# Create assembly ---------------------------------------------------------------------------------------------------------------

c_assembly=c_model.rootAssembly
c_assembly.DatumCsysByDefault(CARTESIAN)

# -Sectors
s1_instance=c_assembly.Instance(
	dependent=ON,
	name='sector-1',
	part=sector_part
	)
c_assembly.DatumAxisByPrincipalAxis(
	principalAxis=ZAXIS
	)
s3_instance=c_assembly.RadialInstancePattern(
	axis=(0.0, 0.0, 1.0), 
	instanceList=('sector-1', ),
	number=3, point=(0.0, 0.0, 0.0),
	totalAngle=240.0
	)

s2_instance=s3_instance[0]
s3_instance=s3_instance[1]

s_instance = (s1_instance,s2_instance ,s3_instance)

# -Gusset plate

# --Create the instances
g1_instance=c_assembly.Instance(
	dependent=ON,
	name='gusset-1',
	part=gusset_part
	)
g2_instance=c_assembly.Instance(
	dependent=ON,
	name='gusset-2',
	part=gusset_part
	)
g3_instance=c_assembly.Instance(
	dependent=ON,
	name='gusset-3',
	part=gusset_part
	)

# --Translate them to the right position
g2_instance.translate(
	vector=(0.0, 0.0, (current_l + current_d))
	)
g3_instance.translate(
	vector=(0.0, 0.0, 2*(current_l + current_d))
	)

# Interactions ------------------------------------------------------------------------------------------------------------------

# Create sets node regions to be used for the tie and coupling constraints
# initiate variables to store points for findAt

holes11=()
holes12=()
holes21=()
holes22=()
holes31=()
holes32=()
gholes1=()
gholes2=()
gholes3=()

# Position of the holes on the cross-section
sh11 = np.array([profiles[i][j][k][0][1], profiles[i][j][k][1][1]])
sh12 = np.array([profiles[i][j][k][0][-2], profiles[i][j][k][1][-2]])

gh1 = (gp1[0]-bolts_w*cos(5*pi/6), gp1[1]-bolts_w*sin(5*pi/6))
gh2 = (gp2[0]-bolts_w*cos(-pi/2), gp2[1]-bolts_w*sin(-pi/2))
gh3 = (gp3[0]-bolts_w*cos(pi/6), gp3[1]-bolts_w*sin(pi/6))

gh=(gh1, gh2, gh3)

# Rotation matrix to multiply the previous point in order to get the points of the other 2 gusset fins
Rmat = np.array([[cos(-2*pi/3), -sin(-2*pi/3)], [sin(-2*pi/3), cos(-2*pi/3)]])

# Calculate the end points of the other 2 gusset fins by multiplying with the 120 degrees rotation matrix
sh21 = sh11.dot(Rmat)
sh22 = sh12.dot(Rmat)

sh31 = sh21.dot(Rmat)
sh32 = sh22.dot(Rmat)

sh = ((sh11, sh12), (sh21, sh22), (sh31, sh32))

# Create reference points for the bolt ridig body couplings

# Create the necessary sets and the tie constraints for all the bolts

# End 1 connection
for oo in (range(3)):
	ii=1
	for o in tuple(bolts_z1):
		
		c_assembly.ReferencePoint((gh[oo-3][0], gh[oo-3][1], float(o)))
		
		c_assembly.Set(
			edges=s_instance[oo-3].edges.findAt(((sh[oo-3][0][0], sh[oo-3][0][1], float(o)-d_washer/2), ), )+\
			s_instance[oo-2].edges.findAt(((sh[oo-2][1][0], sh[oo-2][1][1], float(o)-d_washer/2), ), )+\
			g1_instance.edges.findAt(((gh[oo-3][0], gh[oo-3][1], float(o)-d_washer/2), ), ),
			name='b'+str(ii)+str(oo)+'set1'
			)
		
		c_model.RigidBody(
			name='b1'+str(ii)+str(oo)+'joint1',
			refPointRegion=Region(referencePoints=(c_assembly.referencePoints.findAt((gh[oo-3][0], gh[oo-3][1], float(o))), )),
		    tieRegion=c_assembly.sets['b'+str(ii)+str(oo)+'set1']
			)
		
		ii+=1
			
# Span 1

for oo in (range(3)):
	ii=1
	for o in tuple(bolts_z2):
		
		c_assembly.ReferencePoint((gh[oo-3][0], gh[oo-3][1], float(o)))
		
		c_assembly.Set(
			edges=s_instance[oo-3].edges.findAt(((sh[oo-3][0][0], sh[oo-3][0][1], float(o)-d_washer/2), ), )+\
			s_instance[oo-2].edges.findAt(((sh[oo-2][1][0], sh[oo-2][1][1], float(o)-d_washer/2), ), ),
			name='b'+str(ii)+str(oo)+'-set2'
			)
		
		c_model.RigidBody(
			name='b1'+str(ii)+str(oo)+'span1',
			refPointRegion=Region(referencePoints=(c_assembly.referencePoints.findAt((gh[oo-3][0], gh[oo-3][1], float(o))), )),
		    tieRegion=c_assembly.sets['b'+str(ii)+str(oo)+'-set2']
			)
		
		ii+=1

# middle connection


for oo in (range(3)):
	ii=1
	for o in tuple(bolts_z3):
		
		c_assembly.ReferencePoint((gh[oo-3][0], gh[oo-3][1], float(o)))
		
		c_assembly.Set(
			edges=s_instance[oo-3].edges.findAt(((sh[oo-3][0][0], sh[oo-3][0][1], float(o)-d_washer/2), ), )+\
			s_instance[oo-2].edges.findAt(((sh[oo-2][1][0], sh[oo-2][1][1], float(o)-d_washer/2), ), )+\
			g2_instance.edges.findAt(((gh[oo-3][0], gh[oo-3][1], float(o)-d_washer/2), ), ),
			name='b'+str(ii)+str(oo)+'set3'
			)
		
		c_model.RigidBody(
			name='b1'+str(ii)+str(oo)+'joint2',
			refPointRegion=Region(referencePoints=(c_assembly.referencePoints.findAt((gh[oo-3][0], gh[oo-3][1], float(o))), )),
		    tieRegion=c_assembly.sets['b'+str(ii)+str(oo)+'set3']
			)
		
		ii+=1

# Span 2

for oo in (range(3)):
	ii=1
	for o in tuple(bolts_z4):
		
		c_assembly.ReferencePoint((gh[oo-3][0], gh[oo-3][1], float(o)))
		
		c_assembly.Set(
			edges=s_instance[oo-3].edges.findAt(((sh[oo-3][0][0], sh[oo-3][0][1], float(o)-d_washer/2), ), )+\
			s_instance[oo-2].edges.findAt(((sh[oo-2][1][0], sh[oo-2][1][1], float(o)-d_washer/2), ), ),
			name='b'+str(ii)+str(oo)+'-set4'
			)
		
		c_model.RigidBody(
			name='b1'+str(ii)+str(oo)+'span2',
			refPointRegion=Region(referencePoints=(c_assembly.referencePoints.findAt((gh[oo-3][0], gh[oo-3][1], float(o))), )),
		    tieRegion=c_assembly.sets['b'+str(ii)+str(oo)+'-set4']
			)
		
		ii+=1

# End 2 connection
for oo in (range(3)):
	ii=1
	for o in tuple(bolts_z5):
		
		c_assembly.ReferencePoint((gh[oo-3][0], gh[oo-3][1], float(o)))
		
		c_assembly.Set(
			edges=s_instance[oo-3].edges.findAt(((sh[oo-3][0][0], sh[oo-3][0][1], float(o)-d_washer/2), ), )+\
			s_instance[oo-2].edges.findAt(((sh[oo-2][1][0], sh[oo-2][1][1], float(o)-d_washer/2), ), )+\
			g3_instance.edges.findAt(((gh[oo-3][0], gh[oo-3][1], float(o)-d_washer/2), ), ),
			name='b'+str(ii)+str(oo)+'set5'
			)
		
		c_model.RigidBody(
			name='b1'+str(ii)+str(oo)+'joint3',
			refPointRegion=Region(referencePoints=(c_assembly.referencePoints.findAt((gh[oo-3][0], gh[oo-3][1], float(o))), )),
		    tieRegion=c_assembly.sets['b'+str(ii)+str(oo)+'set5']
			)
		
		ii+=1

# Create reference points for BCs/loads.

# -RPs for the faces at the two ends of the columns
c_assembly.ReferencePoint((0.0, 0.0, 0.0))
c_assembly.ReferencePoint((0.0, 0.0, (2*current_l + 3*current_d)))

# - RP at the middle
c_assembly.ReferencePoint((0.0, 0.0, (current_l + 1.5*current_d)))


# - End face couplings to reference points

# End 1
c_assembly.Set(
	name='RP-1-set', 
	referencePoints=(c_assembly.referencePoints.findAt((0, 0, 0)), )
	)

c_assembly.Set(
	edges=g1_instance.edges.getByBoundingBox(-current_d,-current_d,0,current_d,current_d,0)+\
	s_instance[0].edges.getByBoundingBox(-current_d,-current_d,0,current_d,current_d,0)+\
	s_instance[1].edges.getByBoundingBox(-current_d,-current_d,0,current_d,current_d,0)+\
	s_instance[2].edges.getByBoundingBox(-current_d,-current_d,0,current_d,current_d,0),
	name='end1-face',
	)

c_model.Coupling(
	controlPoint=c_assembly.sets['RP-1-set'], 
	couplingType=KINEMATIC,
	influenceRadius=WHOLE_SURFACE,
	localCsys=None,
	name='end1-coupling', 
	surface=c_assembly.sets['end1-face'],
	u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
	)

# End 2

c_assembly.Set(
	name='RP-2-set',
	referencePoints=(c_assembly.referencePoints.findAt((0, 0, 2*(current_l+1.5*current_d))), )
	)

c_assembly.Set(
	edges=g3_instance.edges.getByBoundingBox(-current_d,-current_d,2*(current_l+1.5*current_d),current_d,current_d,2*(current_l+1.5*current_d))+\
	s_instance[0].edges.getByBoundingBox(-current_d,-current_d,2*(current_l+1.5*current_d),current_d,current_d,2*(current_l+1.5*current_d))+\
	s_instance[1].edges.getByBoundingBox(-current_d,-current_d,2*(current_l+1.5*current_d),current_d,current_d,2*(current_l+1.5*current_d))+\
	s_instance[2].edges.getByBoundingBox(-current_d,-current_d,2*(current_l+1.5*current_d),current_d,current_d,2*(current_l+1.5*current_d)),
	name='end2-face'
	)
	
c_model.Coupling(
	controlPoint=c_assembly.sets['RP-2-set'], 
	couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE,
	localCsys=None,
	name='end2-coupling', 
	surface=c_assembly.sets['end2-face'],
	u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
	)

# Middle

c_assembly.Set(
	name='RP-Mid-set', 
	referencePoints=(c_assembly.referencePoints.findAt((0.0, 0.0, (current_l + 1.5*current_d))), )
	)
	
c_assembly.Set(
	edges=g2_instance.edges.findAt(((0, 0, (current_l + 1.5*current_d)), ), ),
	name='gusset-fin-interface',
	)

c_model.Coupling(
	controlPoint=c_assembly.sets['RP-Mid-set'], 
	couplingType=KINEMATIC,
	influenceRadius=WHOLE_SURFACE,
	localCsys=None,
	name='Mid-coupling', 
	surface=c_assembly.sets['gusset-fin-interface'],
	u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
	)
	
# Step -----------------------------------------------------------------------------------

c_model.BuckleStep(
	maxIterations=300,
	name='Buckling',
	numEigen=4,
	previous='Initial',
	vectors=10
	)

# Boundary Conditions --------------------------------------------------------------

# BCs					
end1_BC=c_model.DisplacementBC(
	amplitude=UNSET, 
	createStepName='Initial', 
	distributionType=UNIFORM, 
	fieldName='', 
	localCsys=None, 
	name='fix-end1', 
	region=Region(referencePoints=(c_assembly.referencePoints.findAt((0, 0, 0)), )), 
	u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=SET
	)
	
end2_BC=c_model.DisplacementBC(
	amplitude=UNSET,
	createStepName='Initial', 
	distributionType=UNIFORM, 
	fieldName='', 
	localCsys=None, 
	name='fix-end2', 
	region=Region(referencePoints=(c_assembly.referencePoints.findAt((0, 0, 2*(current_l+1.5*current_d))), )), 
	u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=SET
	)

middle_BC=c_model.DisplacementBC(
	amplitude=UNSET,
	createStepName='Initial', 
	distributionType=UNIFORM, 
	fieldName='', 
	localCsys=None, 
	name='fix-middle', 
	region=Region(referencePoints=(c_assembly.referencePoints.findAt((0, 0, current_l+1.5*current_d)), )), 
	u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET
	)


# Load -------------------------------------------------------------------------------

c_model.ConcentratedForce(
	cf3=-1.0,
	createStepName='Buckling',
	distributionType=UNIFORM,
	field='',
	localCsys=None,
	name='compression',
	region=c_assembly.sets['RP-2-set']
	)

# Create the job -------------------------------------------------------------------

c_job=mdb.Job(
	atTime=None,
	contactPrint=OFF,
	description='',
	echoPrint=OFF, 
    explicitPrecision=SINGLE,
	getMemoryFromAnalysis=True,
	historyPrint=OFF, 
    memory=90,
	memoryUnits=PERCENTAGE,
	model=buckle_model,
	modelPrint=OFF, 
    multiprocessingMode=DEFAULT,
	name=buckle_model,
	nodalOutputPrecision=SINGLE, 
    numCpus=1,
	numGPUs=0,
	queue=None,
	resultsFormat=ODB,
	scratch='',
	type=ANALYSIS,
	userSubroutine='',
	waitHours=0,
	waitMinutes=0
	)

# Edit the keywords to output translations on the output file
c_model.keywordBlock.synchVersions(storeNodesAndElements=False)
c_model.keywordBlock.insert(GetBlockPosition(c_model,'*End Step')-1, '*NODE FILE\nU')

# Write the input file
mdb.jobs[buckle_model].writeInput()

# RIKS model, Only axial ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

riks_model_N = 'RIKS-N-'+str(i+1)+'-'+str(j+1)+'-'+str(k+1)+'-'+str(l+1)

# Copy model from buckling analysis
r_model_N=mdb.Model(
	name=riks_model_N,
	objectToCopy=c_model
	)

# Delete buckling step
del r_model_N.steps['Buckling']

# Create RIKS step
r_model_N.StaticRiksStep(
	name='RIKS',
	previous='Initial',
	nlgeom=ON,
	maxNumInc=30,
	initialArcInc=0.2
	)

# Change to plastic material, optim355
r_model_N.materials['optim355'].Plastic(table=((381.1, 0.0), (
    391.2, 0.0053), (404.8, 0.0197), (418.0, 0.0228), (444.2, 0.0310), (499.8, 
    0.0503), (539.1, 0.0764), (562.1, 0.1009), (584.6, 0.1221), (594.4, 
    0.1394)))

# Apply concentrated force
N_pl_rd = 510*area

r_model_N.ConcentratedForce(
	cf3=-N_pl_rd,
	createStepName='RIKS',
	distributionType=UNIFORM,
	field='',
	localCsys=None,
	name='compression',
	region=r_model_N.rootAssembly.sets['RP-2-set']
	)

# Field and History output requests

r_model_N.historyOutputRequests.changeKey(
	fromName='H-Output-1',
	toName='load'
	)
	
r_model_N.historyOutputRequests['load'].setValues(
	rebar=EXCLUDE,
	region=r_model_N.rootAssembly.sets['RP-1-set'], 
    sectionPoints=DEFAULT, variables=('RF3', )
	)
	
r_model_N.HistoryOutputRequest(
	createStepName='RIKS',
	name='disp',
	rebar=EXCLUDE,
	region=r_model_N.rootAssembly.sets['RP-2-set'],
	sectionPoints=DEFAULT,
	variables=('U3', )
	)

r_model_N.HistoryOutputRequest(
	createStepName='RIKS',
	name='moment',
	rebar=EXCLUDE,
	region=r_model_N.rootAssembly.sets['RP-Mid-set'],
	sectionPoints=DEFAULT,
	variables=('UR1', )
	)

r_model_N.fieldOutputRequests.changeKey(
	fromName='F-Output-1', 
    toName='fields'
	)
r_model_N.fieldOutputRequests['fields'].setValues(
	variables=('S', 'MISES', 'E', 'PEEQ', 'U')
	)

# Delete keyword nodefile
r_model_N.keywordBlock.synchVersions(storeNodesAndElements=False)
r_model_N.keywordBlock.replace(GetBlockPosition(r_model_N,'*End Step')-1, '\n')

# Change keywords to include initial imperfections file (filename was given wrong initially and corrected later)
amp_impf = s/2000			
#r_model_N.keywordBlock.synchVersions(storeNodesAndElements=False)
r_model_N.keywordBlock.replace(GetBlockPosition(r_model_N, '*step')-1, 
'\n** ----------------------------------------------------------------\n** \n**********GEOMETRICAL IMPERFECTIONS\n*IMPERFECTION,FILE='
+ str(buckle_model) +',STEP=1\n1,'+ str(float(amp_impf)) +'\n2,'+ str(float(amp_impf)) +'\n3,'+ str(float(amp_impf)) +'\n4,'+ str(float(amp_impf)) +'\n**')

# Create Job
mdb.Job(
	atTime=None,
	contactPrint=OFF,
	description='',
	echoPrint=OFF, 
    explicitPrecision=SINGLE,
	getMemoryFromAnalysis=True,
	historyPrint=OFF, 
    memory=90,
	memoryUnits=PERCENTAGE,
	model='RIKS-N-1-1-1-1',
	modelPrint=OFF, 
    multiprocessingMode=DEFAULT,
	name='RIKS-N-1-1-1-1',
	nodalOutputPrecision=SINGLE,
	numCpus=1,
	numGPUs=0,
	queue=None,
	resultsFormat=ODB,
	scratch='', 
    type=ANALYSIS,
	userSubroutine='',
	waitHours=0,
	waitMinutes=0
	)
	
# Write the input file
mdb.jobs[riks_model_N].writeInput()

# RIKS model, Axial snd bending ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

riks_model_NM = 'RIKS-NM-'+str(i+1)+'-'+str(j+1)+'-'+str(k+1)+'-'+str(l+1)

# Copy model from buckling analysis
r_model_NM=mdb.Model(
	name=riks_model_NM,
	objectToCopy=r_model_N
	)

# Apply bending moment at the mid-connection
# Calculate the magnitude of moment as 10% of moment resistance
W = current_Iy/(current_d/2)
M_resist = W*current_fy
M = 0.1*M_resist

r_model_NM.Moment(
	cm1=-M,
	createStepName='RIKS',
	distributionType=UNIFORM,
	field='',
	localCsys=None,
	name='moment',
	region=c_assembly.sets['RP-Mid-set']
	)

# Create Job
mdb.Job(
	atTime=None,
	contactPrint=OFF,
	description='',
	echoPrint=OFF, 
    explicitPrecision=SINGLE,
	getMemoryFromAnalysis=True,
	historyPrint=OFF, 
    memory=90,
	memoryUnits=PERCENTAGE,
	model='RIKS-NM-1-1-1-1',
	modelPrint=OFF, 
    multiprocessingMode=DEFAULT,
	name='RIKS-NM-1-1-1-1',
	nodalOutputPrecision=SINGLE,
	numCpus=1,
	numGPUs=0,
	queue=None,
	resultsFormat=ODB,
	scratch='', 
    type=ANALYSIS,
	userSubroutine='',
	waitHours=0,
	waitMinutes=0
	)
	
# Write the input file
mdb.jobs[riks_model_NM].writeInput()


# Save the model -------------------------------------------------------------------------------------------------------
mdb.saveAs(pathName=os.getcwd()+'\\'+str(i+1)+'-'+str(j+1)+'-'+str(k+1)+'-'+str(l+1)+'.cae')

