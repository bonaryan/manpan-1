# # Import necessary libraries --------------------------------------------------------------------

import numpy as np
import sys
import os
import meshEdit
import abq_toolset as xtr
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
from connectorBehavior import *
from shutil import copyfile
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)


# # Input variables -------------------------------------------------------------------------------

# Number of corners
n = 6

# Diameter of prescribed circle in mm
d = 1000.

# Cross-section slenderness defined by the prescribed circle lambda1=(d/(t^2*epsilon))
cs_slenderness = 120.

# Member slenderness for overall column buckling lambda2= sqrt(A*fy/Ncr)
mb_slenderness = 1.

# Radius if the bended corners of the polygon given as a ratio to the thickness r=rcoef*t
# It regards to the bends of the polygon. The arc radious of the lips' bends is half this value
rcoef = 6.

# Number of elements along the arc length of the bended corners
nbend = 3

# length of the lips given as a ratio to the prescribed circle diameter l=d*l_ratio
l_ratio = 0.14

# Thickness of the gusset plates given as a ratio to the profile thickness tgusset=t_ratio*t
t_ratio = 1.2

# Yield strength in MPa. Used for epsilon, not for the modelling material properties
fy = 355.

# Young's modulus
E_young = 210000.

# Bolt spacing given as a ratio to the prescribed circle diameter, b=s/d.
b = 1

# Imperfection factor for overall bowing of the column u1=l/impamp
flx_imp = 250

# Direction angle for the overall bowing, 0 is on global y axis in rads
theta = pi/2

# Imperfection factor for distortional imperfections u2=s/dist_imp
dist_imp = 250

# Bolt diameter in mm
M_bolt = 16

# Clearence from the washer to the edge of the lip and the start of the bending arc in mm
clearence = 3

# Calculated characteristics
R = d/2
epsilon = sqrt(fy/235)
t = (epsilon**2 * d / cs_slenderness)
tg = (t_ratio*t)
d_washer = ceil(1.5893*M_bolt+5.1071)

# Model specific ID string. Used for save filename and for the jobnames
IDstring = str(int(n))+'-'+str(int(d))+'-'+str(int(b))+'-'+str(int(cs_slenderness))+'-'+str(int(100*mb_slenderness))+'-'+str(int(fy))

# Make a new subdirectory for the current session
os.mkdir(IDstring)

# Copy necessary files to the new directory

copyfile('abq_toolset.py', '.\\'+IDstring+'\\abq_toolset.py')

# Change the working directory
os.chdir('.\\'+IDstring)

# Create a new model ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

static_model_name = 'IMP'

stc_mdl = mdb.Model(modelType=STANDARD_EXPLICIT, name=static_model_name)


# Calculate cross section coordinates for one sector --------------------------------------------------------------------------------

# Angle corresponding to one face of the polygon
theta = 2*pi/n;

# Angles of radii (measured from x-axis)
phi=np.linspace(5*pi/6, pi/6, n/3+1)

# xy coords of the polygon's corners
x = R*np.cos(phi);
y = R*np.sin(phi);

## Bends
# Bending radius
rbend = rcoef*t

# Distance between bending centre and corner
lc = rbend/np.cos(theta/2)

# Centers of bending arcs
xc  = x[1:-1] - lc*np.cos(phi[1:-1])
yc  = y[1:-1] - lc*np.sin(phi[1:-1])

# Bending arc angle
theta_b = pi - theta

# Angles of the edges' midlines (measured from x-axis)
phi_mids = phi[0:-1] - theta/2 

# xy coords of the arc's points
xarc = [[0 for j in range(nbend+1)] for i in range(int(n/3 -1))]
yarc = [[0 for j in range(nbend+1)] for i in range(int(n/3 -1))] 
for i in range(int(n/3 -1)):
    for j in range(nbend+1):
        xarc[i][j] = xc[i] + rbend*np.cos(phi_mids[i]-(j)*(theta/nbend))
        yarc[i][j] = yc[i] + rbend*np.sin(phi_mids[i]-(j)*(theta/nbend))

## Start-end extensions
# Bending radius
rs = rbend/2
xcs = [0, 0]
ycs = [0, 0]

# First bend
v1 = phi_mids[0]-pi/2
v2 = (phi[0]+phi_mids[0]-pi/2)/2
l1 = (t+tg)/(2*np.cos(phi[0]-phi_mids[0]))
l2 = rs/np.sin(v2-phi_mids[0]+pi/2)
x1 = x[0]+l1*np.cos(v1)
y1 = y[0]+l1*np.sin(v1)

# First bend centre coords
xcs[0] = x1+l2*np.cos(v2)
ycs[0] = y1+l2*np.sin(v2)

# Last bend
v1 = phi_mids[-1]+pi/2
v2 = (v1+phi[-1])/2
l1 = (t+tg)/(2*np.cos(v1-phi[-1]-pi/2))
l2 = rs/np.sin(v2-phi[-1])
x1 = x[-1]+l1*np.cos(v1)
y1 = y[-1]+l1*np.sin(v1)

# Last bend centre coords
xcs[1] = x1+l2*np.cos(v2)
ycs[1] = y1+l2*np.sin(v2)

# First and last bend arc points coords
xsarc = [[0 for j in range(nbend+1)] for j in [0,1]]
ysarc = [[0 for j in range(nbend+1)] for j in [0,1]] 
for j in range(nbend+1):
    xsarc[0][j] = xcs[0] + rs*np.cos(4*pi/3+(j)*((phi_mids[0]-pi/3)/nbend))
    ysarc[0][j] = ycs[0] + rs*np.sin(4*pi/3+(j)*((phi_mids[0]-pi/3)/nbend))
    xsarc[1][j] = xcs[1] + rs*np.cos(phi_mids[-1]+pi+(j)*((phi[-1]+pi/2-phi_mids[-1])/nbend))
    ysarc[1][j] = ycs[1] + rs*np.sin(phi_mids[-1]+pi+(j)*((phi[-1]+pi/2-phi_mids[-1])/nbend))


## Points of the lips

# Lip length according to bolt washer diameter
l_lip = d_washer+(2*clearence)

# First lip
xstart = [xsarc[0][0] + l_lip*np.cos(phi[0]), xsarc[0][0] + l_lip*np.cos(phi[0])/2]
ystart = [ysarc[0][0] + l_lip*np.sin(phi[0]), ysarc[0][0] + l_lip*np.sin(phi[0])/2]


# Last point
xend = [xsarc[1][-1] + l_lip*np.cos(phi[-1])/2, xsarc[1][-1] + l_lip*np.cos(phi[-1])]
yend = [ysarc[1][-1] + l_lip*np.sin(phi[-1])/2, ysarc[1][-1] + l_lip*np.sin(phi[-1])]

## Collect the x, y values in a sorted 2xn array
xarcs, yarcs=[],[]
for i in range(len(phi)-2):
	xarcs=xarcs+xarc[i][:]
	yarcs=yarcs+yarc[i][:]

x_sector = xstart+xsarc[0][:]+xarcs[:]+xsarc[1][:]+xend
y_sector = ystart+ysarc[0][:]+yarcs[:]+ysarc[1][:]+yend

# Calculate cross sectional characteristics ----------------------------------------------------------------------------------------

l=0
for i in range(len(x_sector)-1):
    l = l+sqrt((x_sector[i+1]-x_sector[i])**2+(y_sector[i+1]-y_sector[i])**2)
    Area = 3*t*l

I_prexcribed = pi * (d**4 - (d - 2 * t)**4) / 64
I_inscribed = pi * ((d*cos(theta/2))**4 - (d*cos(theta/2) - 2 * t)**4) / 64

Iy = (pi * (d**4 - (d - 2 * t)**4) / 64) * l/(pi*d)

length = mb_slenderness*pi*sqrt(E_young*Iy/(Area*fy))

# Create Parts ----------------------------------------------------------------------------------------------------------------------

# Sector

# -Profile sketch for sector
sector_sketch = stc_mdl.ConstrainedSketch(name='sector', sheetSize=1200.0)

# -Sketch sector lines
for n in range(len(x_sector)-1):
	sector_sketch.Line(
		point1=(x_sector[n], y_sector[n]), 
		point2=(x_sector[n+1], y_sector[n+1])
		)

# -Extrude sector part
l_tot = 2*length + 3*d
sector_part = stc_mdl.Part(
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
bolts_w = l_lip/2

# -Distances on the length
s = b*d
(n0, s0) = divmod(length, s)
s1 = (s0 + s)/2

bolts_z1 = np.concatenate([[bolts_w], bolts_w + ((d - l_lip)/5) * np.linspace(1, 4, 4), [d - bolts_w]])
bolts_z2 = np.concatenate([[d + s1], s1 + d + (s * np.linspace(1, n0-1, n0-1))])
bolts_z3 = bolts_z1 + (length + d)
bolts_z4 = bolts_z2 + (length + d)
bolts_z5 = bolts_z3 + (length + d)

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
		edge1=sector_part.edges.getClosest(coordinates=((x_sector[1],y_sector[1], 0),))[0][0],
		edge2=sector_part.edges.getClosest(coordinates=((x_sector[0], y_sector[0], 1),))[0][0],
		plane=sector_part.faces.getClosest(coordinates=((x_sector[0], y_sector[0], 0),))[0][0],
		planeSide=SIDE1
		)
	
	sector_part.HoleBlindFromEdges(
		depth=1.0,
		diameter=d_washer,
		distance1=bolts_z[o],
		distance2=bolts_w,
		edge1=sector_part.edges.getClosest(coordinates=((x_sector[-2], y_sector[-2], 0),))[0][0],
		edge2=sector_part.edges.getClosest(coordinates=((x_sector[-1], y_sector[-1], 1),))[0][0],
		plane=sector_part.faces.getClosest(coordinates=((x_sector[-1], y_sector[-1], 0),))[0][0],
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
	sector_part.PartitionFaceByDatumPlane(
		datumPlane=sector_part.datums.items()[o+1][1],
		faces=sector_part.faces[:]
		)

# Gusset

# -Profile sketch for gusset
gusset_sketch=stc_mdl.ConstrainedSketch(name='__profile__', sheetSize=1200.0)

# -Sketch gusset lines
# First point of the first sector
x0 = x_sector[0]
y0 = y_sector[0]

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
gusset_part=stc_mdl.Part(
	dimensionality=THREE_D,
	name='gusset',
	type=DEFORMABLE_BODY
	)
gusset_part.BaseShellExtrude(
	depth=d,
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

gusset_part.DatumPointByCoordinate((gp1[0]-l_lip*cos(5*pi/6), gp1[1]-l_lip*sin(5*pi/6), 0),)
gusset_part.DatumPointByCoordinate((gp1[0]-l_lip*cos(5*pi/6), gp1[1]-l_lip*sin(5*pi/6), d),)

gusset_part.PartitionFaceByShortestPath(
	faces=gusset_part.faces.getClosest(coordinates=((gp1[0], gp1[1], 0),))[0][0],
	point1=gusset_part.datum.items()[0][1],
	point2=gusset_part.datum.items()[1][1],
	)

gusset_part.DatumPointByCoordinate((gp2[0]-l_lip*cos(-pi/2), gp2[1]-l_lip*sin(-pi/2), 0),)
gusset_part.DatumPointByCoordinate((gp2[0]-l_lip*cos(-pi/2), gp2[1]-l_lip*sin(-pi/2), d),)

gusset_part.PartitionFaceByShortestPath(
	faces=gusset_part.faces.getClosest(coordinates=((gp2[0], gp2[1], 0),))[0][0],
	point1=gusset_part.datum.items()[2][1],
	point2=gusset_part.datum.items()[3][1],
	)

gusset_part.DatumPointByCoordinate((gp3[0]-l_lip*cos(pi/6), gp3[1]-l_lip*sin(pi/6), 0),)
gusset_part.DatumPointByCoordinate((gp3[0]-l_lip*cos(pi/6), gp3[1]-l_lip*sin(pi/6), d),)

gusset_part.PartitionFaceByShortestPath(
	faces=gusset_part.faces.getClosest(coordinates=((gp3[0], gp3[1], 0),))[0][0],
	point1=gusset_part.datum.items()[4][1],
	point2=gusset_part.datum.items()[5][1],
	)

# Material ----------------------------------------------------------------------------------------------------------------------

stc_mdl.Material(name='elastic')
stc_mdl.materials['elastic'].Elastic(table=((E_young, 0.3), ))

# Create sections ---------------------------------------------------------------------------------------------------------------

# -for sector
stc_mdl.HomogeneousShellSection(
	idealization=NO_IDEALIZATION, 
	integrationRule=SIMPSON,
	material='elastic',
	name='sector',
	numIntPts=5,
	poissonDefinition=DEFAULT,
	preIntegrate=OFF,
	temperature=GRADIENT, 
	thickness=t,
	thicknessField='',
	thicknessModulus=None,
	thicknessType=UNIFORM,
	useDensity=OFF
	)
	
# -for gusset
stc_mdl.HomogeneousShellSection(
	idealization=NO_IDEALIZATION, 
	integrationRule=SIMPSON,
	material='elastic',
	name='gusset',
	numIntPts=5,
	poissonDefinition=DEFAULT,
	preIntegrate=OFF,
	temperature=GRADIENT, 
	thickness=tg,
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

assmbl=stc_mdl.rootAssembly
assmbl.DatumCsysByDefault(CARTESIAN)

# -Sectors
s1_instance=assmbl.Instance(
	dependent=ON,
	name='sector-1',
	part=sector_part
	)
assmbl.DatumAxisByPrincipalAxis(
	principalAxis=ZAXIS
	)
s3_instance=assmbl.RadialInstancePattern(
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
g1_instance=assmbl.Instance(
	dependent=ON,
	name='gusset-1',
	part=gusset_part
	)
g2_instance=assmbl.Instance(
	dependent=ON,
	name='gusset-2',
	part=gusset_part
	)
g3_instance=assmbl.Instance(
	dependent=ON,
	name='gusset-3',
	part=gusset_part
	)

# --Translate them to the right position
g2_instance.translate(
	vector=(0.0, 0.0, (length + d))
	)
g3_instance.translate(
	vector=(0.0, 0.0, 2*(length + d))
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
sh11 = np.array([x_sector[1], y_sector[1]])
sh12 = np.array([x_sector[-2], y_sector[-2]])

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
		
		assmbl.ReferencePoint((gh[oo-3][0], gh[oo-3][1], float(o)))
		
		assmbl.Set(
			edges=s_instance[oo-3].edges.findAt(((sh[oo-3][0][0], sh[oo-3][0][1], float(o)-d_washer/2), ), )+\
			s_instance[oo-2].edges.findAt(((sh[oo-2][1][0], sh[oo-2][1][1], float(o)-d_washer/2), ), )+\
			g1_instance.edges.findAt(((gh[oo-3][0], gh[oo-3][1], float(o)-d_washer/2), ), ),
			name='b'+str(ii)+str(oo)+'set1'
			)
		
		stc_mdl.RigidBody(
			name='b1'+str(ii)+str(oo)+'joint1',
			refPointRegion=Region(referencePoints=(assmbl.referencePoints.findAt((gh[oo-3][0], gh[oo-3][1], float(o))), )),
		    tieRegion=assmbl.sets['b'+str(ii)+str(oo)+'set1']
			)
		
		ii+=1
			
# Span 1

for oo in (range(3)):
	ii=1
	for o in tuple(bolts_z2):
		
		assmbl.ReferencePoint((gh[oo-3][0], gh[oo-3][1], float(o)))
		
		assmbl.Set(
			edges=s_instance[oo-3].edges.findAt(((sh[oo-3][0][0], sh[oo-3][0][1], float(o)-d_washer/2), ), )+\
			s_instance[oo-2].edges.findAt(((sh[oo-2][1][0], sh[oo-2][1][1], float(o)-d_washer/2), ), ),
			name='b'+str(ii)+str(oo)+'-set2'
			)
		
		stc_mdl.RigidBody(
			name='b1'+str(ii)+str(oo)+'span1',
			refPointRegion=Region(referencePoints=(assmbl.referencePoints.findAt((gh[oo-3][0], gh[oo-3][1], float(o))), )),
		    tieRegion=assmbl.sets['b'+str(ii)+str(oo)+'-set2']
			)
		
		ii+=1

# middle connection


for oo in (range(3)):
	ii=1
	for o in tuple(bolts_z3):
		
		assmbl.ReferencePoint((gh[oo-3][0], gh[oo-3][1], float(o)))
		
		assmbl.Set(
			edges=s_instance[oo-3].edges.findAt(((sh[oo-3][0][0], sh[oo-3][0][1], float(o)-d_washer/2), ), )+\
			s_instance[oo-2].edges.findAt(((sh[oo-2][1][0], sh[oo-2][1][1], float(o)-d_washer/2), ), )+\
			g2_instance.edges.findAt(((gh[oo-3][0], gh[oo-3][1], float(o)-d_washer/2), ), ),
			name='b'+str(ii)+str(oo)+'set3'
			)
		
		stc_mdl.RigidBody(
			name='b1'+str(ii)+str(oo)+'joint2',
			refPointRegion=Region(referencePoints=(assmbl.referencePoints.findAt((gh[oo-3][0], gh[oo-3][1], float(o))), )),
		    tieRegion=assmbl.sets['b'+str(ii)+str(oo)+'set3']
			)
		
		ii+=1

# Span 2

for oo in (range(3)):
	ii=1
	for o in tuple(bolts_z4):
		
		assmbl.ReferencePoint((gh[oo-3][0], gh[oo-3][1], float(o)))
		
		assmbl.Set(
			edges=s_instance[oo-3].edges.findAt(((sh[oo-3][0][0], sh[oo-3][0][1], float(o)-d_washer/2), ), )+\
			s_instance[oo-2].edges.findAt(((sh[oo-2][1][0], sh[oo-2][1][1], float(o)-d_washer/2), ), ),
			name='b'+str(ii)+str(oo)+'-set4'
			)
		
		stc_mdl.RigidBody(
			name='b1'+str(ii)+str(oo)+'span2',
			refPointRegion=Region(referencePoints=(assmbl.referencePoints.findAt((gh[oo-3][0], gh[oo-3][1], float(o))), )),
		    tieRegion=assmbl.sets['b'+str(ii)+str(oo)+'-set4']
			)
		
		ii+=1

# End 2 connection
for oo in (range(3)):
	ii=1
	for o in tuple(bolts_z5):
		
		assmbl.ReferencePoint((gh[oo-3][0], gh[oo-3][1], float(o)))
		
		assmbl.Set(
			edges=s_instance[oo-3].edges.findAt(((sh[oo-3][0][0], sh[oo-3][0][1], float(o)-d_washer/2), ), )+\
			s_instance[oo-2].edges.findAt(((sh[oo-2][1][0], sh[oo-2][1][1], float(o)-d_washer/2), ), )+\
			g3_instance.edges.findAt(((gh[oo-3][0], gh[oo-3][1], float(o)-d_washer/2), ), ),
			name='b'+str(ii)+str(oo)+'set5'
			)
		
		stc_mdl.RigidBody(
			name='b1'+str(ii)+str(oo)+'joint3',
			refPointRegion=Region(referencePoints=(assmbl.referencePoints.findAt((gh[oo-3][0], gh[oo-3][1], float(o))), )),
		    tieRegion=assmbl.sets['b'+str(ii)+str(oo)+'set5']
			)
		
		ii+=1

# Create reference points for BCs/loads.

# -RPs for the faces at the two ends of the columns
assmbl.ReferencePoint((0.0, 0.0, 0.0))
assmbl.ReferencePoint((0.0, 0.0, (2*length + 3*d)))

# - RP at the middle
assmbl.ReferencePoint((0.0, 0.0, (length + 1.5*d)))


# - End face couplings to reference points

# End 1
assmbl.Set(
	name='RP-1-set', 
	referencePoints=(assmbl.referencePoints.findAt((0, 0, 0)), )
	)

assmbl.Set(
	edges=g1_instance.edges.getByBoundingBox(-d,-d,0,d,d,0)+\
	s_instance[0].edges.getByBoundingBox(-d,-d,0,d,d,0)+\
	s_instance[1].edges.getByBoundingBox(-d,-d,0,d,d,0)+\
	s_instance[2].edges.getByBoundingBox(-d,-d,0,d,d,0),
	name='end1-face',
	)

stc_mdl.Coupling(
	controlPoint=assmbl.sets['RP-1-set'], 
	couplingType=KINEMATIC,
	influenceRadius=WHOLE_SURFACE,
	localCsys=None,
	name='end1-coupling', 
	surface=assmbl.sets['end1-face'],
	u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
	)

# End 2

assmbl.Set(
	name='RP-2-set',
	referencePoints=(assmbl.referencePoints.findAt((0, 0, 2*(length+1.5*d))), )
	)

assmbl.Set(
	edges=g3_instance.edges.getByBoundingBox(-d,-d,2*(length+1.5*d),d,d,2*(length+1.5*d))+\
	s_instance[0].edges.getByBoundingBox(-d,-d,2*(length+1.5*d),d,d,2*(length+1.5*d))+\
	s_instance[1].edges.getByBoundingBox(-d,-d,2*(length+1.5*d),d,d,2*(length+1.5*d))+\
	s_instance[2].edges.getByBoundingBox(-d,-d,2*(length+1.5*d),d,d,2*(length+1.5*d)),
	name='end2-face'
	)
	
stc_mdl.Coupling(
	controlPoint=assmbl.sets['RP-2-set'], 
	couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE,
	localCsys=None,
	name='end2-coupling', 
	surface=assmbl.sets['end2-face'],
	u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
	)

# Middle

assmbl.Set(
	name='RP-Mid-set', 
	referencePoints=(assmbl.referencePoints.findAt((0.0, 0.0, (length + 1.5*d))), )
	)
	
assmbl.Set(
	edges=g2_instance.edges.findAt(((0, 0, (length + 1.5*d)), ), ),
	name='gusset-fin-interface',
	)

stc_mdl.Coupling(
	controlPoint=assmbl.sets['RP-Mid-set'], 
	couplingType=KINEMATIC,
	influenceRadius=WHOLE_SURFACE,
	localCsys=None,
	name='Mid-coupling', 
	surface=assmbl.sets['gusset-fin-interface'],
	u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
	)

# Step -----------------------------------------------------------------------------------
stc_mdl.StaticStep(
	name='IMP',
	previous='Initial')


#stc_mdl.StaticRiksStep(
#	name='RIKS',
#	previous='Initial',
#	nlgeom=ON,
#	maxNumInc=10,
#	initialArcInc=0.2
#	)
	
	
# Boundary Conditions --------------------------------------------------------------

# BCs					
end1_BC=stc_mdl.DisplacementBC(
	amplitude=UNSET, 
	createStepName='Initial', 
	distributionType=UNIFORM, 
	fieldName='', 
	localCsys=None, 
	name='fix-end1', 
	region=Region(referencePoints=(assmbl.referencePoints.findAt((0, 0, 0)), )), 
	u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=SET
	)
	
end2_BC=stc_mdl.DisplacementBC(
	amplitude=UNSET,
	createStepName='Initial', 
	distributionType=UNIFORM, 
	fieldName='', 
	localCsys=None, 
	name='fix-end2', 
	region=Region(referencePoints=(assmbl.referencePoints.findAt((0, 0, 2*(length+1.5*d))), )), 
	u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=SET
	)

middle_BC=stc_mdl.DisplacementBC(
	amplitude=UNSET,
	createStepName='Initial', 
	distributionType=UNIFORM, 
	fieldName='', 
	localCsys=None, 
	name='fix-middle', 
	region=Region(referencePoints=(assmbl.referencePoints.findAt((0, 0, length+1.5*d)), )), 
	u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET
	)


# Load -------------------------------------------------------------------------------
# local imperfection amplitude
l_flx_imp = 200
u_dist = s/l_flx_imp

# Start-end points of the 3 sectors
ss11=np.array([x_sector[0], y_sector[0]])
ss12=np.array([x_sector[-1], y_sector[-1]])

#rotate to get start-end ponts of the other two sectors
ss21=ss11.dot(Rmat)
ss22=ss12.dot(Rmat)
ss31=ss21.dot(Rmat)
ss32=ss22.dot(Rmat)

for i in np.concatenate([np.array(range(1, len(bolts_z2))),(len(bolts_z2)+np.array(range(1, len(bolts_z4))))]):

	assmbl.Surface(
		name='lip1'+str(i),
		side1Edges=assmbl.instances['sector-1'].edges.findAt(((ss11[0], ss11[1], (np.concatenate([bolts_z2,bolts_z4]))[i]-2*d_washer), ), )+\
			assmbl.instances['sector-1-rad-2'].edges.findAt(((ss22[0], ss22[1], (np.concatenate([bolts_z2,bolts_z4]))[i]-2*d_washer), ), )+\
			assmbl.instances['sector-1-rad-2'].edges.findAt(((ss21[0], ss21[1], (np.concatenate([bolts_z2,bolts_z4]))[i]-2*d_washer), ), )+\
			assmbl.instances['sector-1-rad-3'].edges.findAt(((ss32[0], ss32[1], (np.concatenate([bolts_z2,bolts_z4]))[i]-2*d_washer), ), )
		)
	
	assmbl.Surface(
		name='lip2'+str(i),
		side1Edges=assmbl.instances['sector-1-rad-3'].edges.findAt(((ss31[0], ss31[1], (np.concatenate([bolts_z2,bolts_z4]))[i]-2*d_washer), ), )+\
			assmbl.instances['sector-1'].edges.findAt(((ss12[0], ss12[1], (np.concatenate([bolts_z2,bolts_z4]))[i]-2*d_washer), ), )
		)
	
	mdb.models['IMP'].ShellEdgeLoad(
		createStepName='IMP',
		distributionType=UNIFORM, 
		field='',
		localCsys=None,
		magnitude=1.0*(-1)**i,
		name='Load-1'+str(i),
		region=mdb.models['IMP'].rootAssembly.surfaces['lip1'+str(i)]
		)
	
	mdb.models['IMP'].ShellEdgeLoad(
		createStepName='IMP',
		distributionType=UNIFORM, 
		field='',
		localCsys=None,
		magnitude=1.0*(-1)**(i+1),
		name='Load-2'+str(i),
		region=mdb.models['IMP'].rootAssembly.surfaces['lip2'+str(i)]
		)

# Field output request. Request displacements 'U'

stc_mdl.fieldOutputRequests.changeKey(
	fromName='F-Output-1', 
    toName='fields'
	)
stc_mdl.fieldOutputRequests['fields'].setValues(
	variables=('U',)
	)


# RIKS model, Only axial ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

riks_model_name = 'RIKS'

# Copy the previous static model
riks_mdl=mdb.Model(
	name=riks_model_name,
	objectToCopy=stc_mdl
	)

# Delete buckling step
del riks_mdl.steps['IMP']

# Create RIKS step
riks_mdl.StaticRiksStep(
	name='RIKS',
	previous='Initial',
	nlgeom=ON,
	maxNumInc=10,
	extrapolation=PARABOLIC
	)

# Rename the material
riks_mdl.materials.changeKey(fromName='elastic', toName='optim355')

# Change to plastic material, optim355
riks_mdl.materials['optim355'].Plastic(
	table=((381.1, 0.0), (
    391.2, 0.0053), (404.8, 0.0197), (418.0, 0.0228), (444.2, 0.0310), (499.8, 
    0.0503), (539.1, 0.0764), (562.1, 0.1009), (584.6, 0.1221), (594.4, 
    0.1394))
	)

# Change the section material name accordingly
riks_mdl.sections['gusset'].setValues(
	idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON,
	material='optim355',
	numIntPts=5,
	preIntegrate=OFF,
	thickness=18.1277,
	thicknessField='',
	thicknessType=UNIFORM
	)

riks_mdl.sections['sector'].setValues(
	idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON,
	material='optim355',
	numIntPts=5,
	preIntegrate=OFF,
	thickness=15.1064,
	thicknessField='',
	thicknessType=UNIFORM
	)

# Overall buckling Imperfections

# Make assmblly instances independent o parts
assmbl_riks=riks_mdl.rootAssembly

assmbl_riks.makeIndependent(instances=(
    assmbl_riks.instances['sector-1'], 
    assmbl_riks.instances['sector-1-rad-2'], 
    assmbl_riks.instances['sector-1-rad-3'], 
    assmbl_riks.instances['gusset-1'], 
    assmbl_riks.instances['gusset-2'], 
    assmbl_riks.instances['gusset-3']))

# Translate nodes for global imperfections

for i in assmbl_riks.allInstances.items():
	for j in range(len(i[1].nodes)):
		xi = i[1].nodes[j].coordinates[0]
		yi = i[1].nodes[j].coordinates[1]
		zi = i[1].nodes[j].coordinates[2]
		assmbl_riks.editNode(
			nodes=i[1].nodes[j],
			offset1=(l_tot/flx_imp)*sin(2*pi*zi/l_tot)*cos(theta),
			offset2=(l_tot/flx_imp)*sin(2*pi*zi/l_tot)*sin(theta)
			)


## The following commented code displaces the three sectors at the bolts' positions
#for i in [2, 4]:
#	for ii in range(1, len(bolts_z2)+1):
#		
#		riks_mdl.DisplacementBC(
#			amplitude=UNSET, 
#			createStepName='IMP',
#			distributionType=UNIFORM,
#			fieldName='',
#			fixed=OFF, 
#			localCsys=None,
#			name='b1'+str(i)+str(ii),
#			region=assmbl_riks.sets['b'+str(ii)+'0-set'+str(i)], 
#			u1=cos(5*pi/6)*u_dist, u2=sin(5*pi/6)*u_dist, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET
#			)
#		
#		riks_mdl.DisplacementBC(
#			amplitude=UNSET, 
#			createStepName='IMP',
#			distributionType=UNIFORM,
#			fieldName='',
#			fixed=OFF, 
#			localCsys=None,
#			name='b2'+str(i)+str(ii),
#			region=assmbl_riks.sets['b'+str(ii)+'1-set'+str(i)], 
#			u1=UNSET, u2=u_dist, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET
#			)
#			
#		riks_mdl.DisplacementBC(
#			amplitude=UNSET, 
#			createStepName='IMP',
#			distributionType=UNIFORM,
#			fieldName='',
#			fixed=OFF, 
#			localCsys=None,
#			name='b3'+str(i)+str(ii),
#			region=assmbl_riks.sets['b'+str(ii)+'2-set'+str(i)], 
#			u1=cos(pi/6)*u_dist, u2=sin(pi/6)*u_dist, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET
#			)


# Apply concentrated force
N_pl_rd = 510*Area

riks_mdl.ConcentratedForce(
	cf3=-N_pl_rd,
	createStepName='RIKS',
	distributionType=UNIFORM,
	field='',
	localCsys=None,
	name='compression',
	region=riks_mdl.rootAssembly.sets['RP-2-set']
	)

# Field and History output requests

riks_mdl.historyOutputRequests.changeKey(
	fromName='H-Output-1',
	toName='load'
	)
	
riks_mdl.historyOutputRequests['load'].setValues(
	rebar=EXCLUDE,
	region=riks_mdl.rootAssembly.sets['RP-1-set'], 
    sectionPoints=DEFAULT, variables=('RF3', )
	)
	
riks_mdl.HistoryOutputRequest(
	createStepName='RIKS',
	name='disp',
	rebar=EXCLUDE,
	region=riks_mdl.rootAssembly.sets['RP-2-set'],
	sectionPoints=DEFAULT,
	variables=('U3', )
	)

riks_mdl.HistoryOutputRequest(
	createStepName='RIKS',
	name='moment',
	rebar=EXCLUDE,
	region=riks_mdl.rootAssembly.sets['RP-Mid-set'],
	sectionPoints=DEFAULT,
	variables=('UR1', )
	)

riks_mdl.fieldOutputRequests.changeKey(
	fromName='F-Output-1', 
    toName='fields'
	)
riks_mdl.fieldOutputRequests['fields'].setValues(
	variables=('S', 'MISES', 'E', 'PEEQ', 'U')
	)

# Edit keywords, create and submit jobs
# The static analysis is first run. Then the maximum displacement is found in the static model
# results and used as a normalization factor for the imperfection amplitude

# definition of a method to search for a keyword position
def GetBlockPosition(model,blockPrefix):
	import string
	pos = 0
	for block in model.keywordBlock.sieBlocks:
		if string.lower(block[0:len(blockPrefix)])==string.lower(blockPrefix):
			return pos
		pos=pos+1
	return -1

# Edit the keywords for the static model to write 'U' magnitude displacements
stc_mdl.keywordBlock.synchVersions(storeNodesAndElements=False)
stc_mdl.keywordBlock.insert(GetBlockPosition(stc_mdl,'*End Step')-1, '*NODE FILE\nU')

# Delete initial model
del mdb.models['Model-1']

# Create and submit the jobs -------------------------------------------------------------------

# Static model
stc_job=mdb.Job(
	atTime=None,
	contactPrint=OFF,
	description='',
	echoPrint=OFF, 
    explicitPrecision=SINGLE,
	getMemoryFromAnalysis=True,
	historyPrint=OFF, 
    memory=90,
	memoryUnits=PERCENTAGE,
	model=static_model_name,
	modelPrint=OFF, 
    multiprocessingMode=DEFAULT,
	name=IDstring+'-imp',
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

# Submit the static job and wait for the results
stc_job.submit() 
stc_job.waitForCompletion()

# find the maximum displacement from the static analysis
stc_odb = open_odb(IDstring+'-imp.odb')
Umax = xtr.max_result(stc_odb, ['U', 'Magnitude'])

# Calculate the imperfection amplitude based on max displacement
a_factor = s/(dist_imp*Umax)

# Edit the keywords for the compression riks model to include imperfections from the static analysis
riks_mdl.keywordBlock.synchVersions(storeNodesAndElements=False)
riks_mdl.keywordBlock.replace(GetBlockPosition(riks_mdl, '*step')-1, 
'\n** ----------------------------------------------------------------\n** \n**********GEOMETRICAL IMPERFECTIONS\n*IMPERFECTION,FILE='
+ IDstring+'-imp' +',STEP=1\n1,'+ str(a_factor)+'\n**')


# Riks model
riks_job=mdb.Job(
	atTime=None,
	contactPrint=OFF,
	description='',
	echoPrint=OFF, 
    explicitPrecision=SINGLE,
	getMemoryFromAnalysis=True,
	historyPrint=OFF, 
    memory=90,
	memoryUnits=PERCENTAGE,
	model=riks_model_name,
	modelPrint=OFF, 
    multiprocessingMode=DEFAULT,
	name=IDstring+'-riks',
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

# Submit the static job and wait for the results
riks_job.submit()

## Save the model -------------------------------------------------------------------------------------------------------
mdb.saveAs(pathName=os.getcwd()+'\\'+IDstring+'.cae')

# Return to parent directory
os.chdir('..')