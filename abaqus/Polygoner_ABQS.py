# -*- coding: mbcs -*-
import pickle
import numpy as np
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
from visualization import *
from connectorBehavior import *

# Import profiles database and profile metadata
profiles_file = open("profiles.pkl",'rb')
profiles = pickle.load(profiles_file)
profiles_file.close()

profiles_file = open("meta.pkl",'rb')
profiles_meta = pickle.load(profiles_file)
profiles_file.close()

## 1st Phase: Buckling analysis

#for i in range(profiles.shape[0]):
#	for j in range(profiles.shape[1]):
#		for k in range(profiles.shape[2]):
#			for l in range(profiles_meta.shape[3])
for i in range(1):
	for j in range(1):
		for k in range(1):
			for l in range(1):
				# Variables holding information of the current profile
				current_model = str(i+1)+'-'+str(j+1)+'-'+str(k+1)+'-'+str(l+1)
				current_d = float(profiles_meta[i][j][k][l][0][0])
				current_t = float(profiles_meta[i][j][k][l][1][0])
				current_tg = float(profiles_meta[i][j][k][l][2][0])
				current_fy = float(profiles_meta[i][j][k][l][3][0])
				current_l = float(profiles_meta[i][j][k][l][7][0])
				current_llip = sqrt((profiles[i][j][k][0][0]-profiles[i][j][k][0][2])**2+(profiles[i][j][k][1][0]-profiles[i][j][k][1][2])**2)
				
				# Create model
				mdb.Model(modelType=STANDARD_EXPLICIT, name=current_model)
	
				# Material
				mdb.models[current_model].Material(name='pure-elastic')
				mdb.models[current_model].materials['pure-elastic'].Elastic(table=((210000.0, 0.3), ))
	
				# Create sections
				# -for sector
				mdb.models[current_model].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
					integrationRule=SIMPSON, material='pure-elastic', name='sector', numIntPts=
					5, poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
					thickness=current_t, thicknessField='', thicknessModulus=None, thicknessType=
					UNIFORM, useDensity=OFF)
				# -for gusset
				mdb.models[current_model].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
					integrationRule=SIMPSON, material='pure-elastic', name='gusset', numIntPts=
					5, poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
					thickness=current_tg, thicknessField='', thicknessModulus=None, thicknessType=
					UNIFORM, useDensity=OFF)
	
				# Create Parts
				
				# Sector
				
				# Loop through the different bolt spacings (temporary b to change)
				b = [0.5, 1, 1.5]
				for m in range(1):
					# -Profile sketch for sector
					mdb.models[current_model].ConstrainedSketch(name='__profile__', sheetSize=1200.0)
					
					# -Sketch sector lines
					for n in range(profiles[i][j][k].shape[1]-1):
						mdb.models[current_model].sketches['__profile__'].Line(point1=(profiles[i][j][k][0][n], profiles[i][j][k][1][n]), 
							point2=(profiles[i][j][k][0][n+1], profiles[i][j][k][1][n+1]))
					
					# -Extrude sector part
					l_tot = 2*current_l + 3*current_d
					mdb.models[current_model].Part(dimensionality=THREE_D, name='sector', type=
						DEFORMABLE_BODY)
					mdb.models[current_model].parts['sector'].BaseShellExtrude(depth=l_tot, sketch=
						mdb.models[current_model].sketches['__profile__'])
					del mdb.models[current_model].sketches['__profile__']
					
					# Calculate bolt positions
					
					# -Distance on the width
					bolts_w = current_llip/2
					
					# -Distances on the length
					current_b =  b[m]
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
										
					# Make holes
					for o in range(int(bolts_z.shape[0])):
#					for o in range(5):
						mdb.models['1-1-1-1'].parts['sector'].HoleBlindFromEdges(depth=1.0, diameter=d_washer
							, distance1=bolts_z[o], distance2=bolts_w, edge1=
							mdb.models['1-1-1-1'].parts['sector'].edges.getClosest(coordinates=((profiles[i][j][k][0][1], profiles[i][j][k][1][1], 0),))[0][0], edge2=
							mdb.models['1-1-1-1'].parts['sector'].edges.getClosest(coordinates=((profiles[i][j][k][0][0], profiles[i][j][k][1][0], 1),))[0][0], plane=
							mdb.models['1-1-1-1'].parts['sector'].faces.getClosest(coordinates=((profiles[i][j][k][0][0], profiles[i][j][k][1][0], 0),))[0][0], planeSide=SIDE1)
						
						mdb.models['1-1-1-1'].parts['sector'].HoleBlindFromEdges(depth=1.0, diameter=d_washer
							, distance1=bolts_z[o], distance2=bolts_w, edge1=
							mdb.models['1-1-1-1'].parts['sector'].edges.getClosest(coordinates=((profiles[i][j][k][0][-2], profiles[i][j][k][1][-2], 0),))[0][0], edge2=
							mdb.models['1-1-1-1'].parts['sector'].edges.getClosest(coordinates=((profiles[i][j][k][0][-1], profiles[i][j][k][1][-1], 1),))[0][0], plane=
							mdb.models['1-1-1-1'].parts['sector'].faces.getClosest(coordinates=((profiles[i][j][k][0][-1], profiles[i][j][k][1][-1], 0),))[0][0], planeSide=SIDE1)
						
						# Create datum planes to be used for partitioning the sector
						mdb.models['1-1-1-1'].parts['sector'].DatumPlaneByPrincipalPlane(offset=bolts_z[o]-bolts_w, 
							principalPlane=XYPLANE)
						mdb.models['1-1-1-1'].parts['sector'].DatumPlaneByPrincipalPlane(offset=bolts_z[o]+bolts_w, 
							principalPlane=XYPLANE)
					
					# Partition the sector
					
					# -Number of datum planes
#					n_dat = int(len(mdb.models['1-1-1-1'].parts['sector'].datums))
					
					# cut all the faces using the datum planes
#					for o in range((n_dat-2)/3):
					for o in range(1):						
						mdb.models['1-1-1-1'].parts['sector'].PartitionFaceByDatumPlane(datumPlane=
							mdb.models['1-1-1-1'].parts['sector'].datums.items()[o+1][1], faces=
							mdb.models['1-1-1-1'].parts['sector'].faces[:])
					
					# Gusset
					
					# -Profile sketch for gusset
					mdb.models[current_model].ConstrainedSketch(name='__profile__', sheetSize=1200.0)
					
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
					mdb.models[current_model].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(gp1[0], gp1[1]))
					mdb.models[current_model].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(gp2[0], gp2[1]))
					mdb.models[current_model].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(gp3[0], gp3[1]))
		
					# -Extrude gusset part
					mdb.models[current_model].Part(dimensionality=THREE_D, name='gusset', type=
						DEFORMABLE_BODY)
					mdb.models[current_model].parts['gusset'].BaseShellExtrude(depth=current_d, sketch=
						mdb.models[current_model].sketches['__profile__'])
					del mdb.models[current_model].sketches['__profile__']
					
					# -Holes
					for o in range(int(bolts_z1.shape[0])):
						mdb.models['1-1-1-1'].parts['gusset'].HoleBlindFromEdges(depth=1.0, diameter=d_washer
							, distance1=bolts_z1[o], distance2=bolts_w, edge1=
							mdb.models['1-1-1-1'].parts['gusset'].edges.getClosest(coordinates=((gp1[0]/2, gp1[1]/2, 0),))[0][0], edge2=
							mdb.models['1-1-1-1'].parts['gusset'].edges.getClosest(coordinates=((gp1[0], gp1[1], 1),))[0][0], plane=
							mdb.models['1-1-1-1'].parts['gusset'].faces.getClosest(coordinates=((gp1[0], gp1[1], 0),))[0][0], planeSide=SIDE1)
						
						mdb.models['1-1-1-1'].parts['gusset'].HoleBlindFromEdges(depth=1.0, diameter=d_washer
							, distance1=bolts_z1[o], distance2=bolts_w, edge1=
							mdb.models['1-1-1-1'].parts['gusset'].edges.getClosest(coordinates=((gp2[0]/2, gp2[1]/2, 0),))[0][0], edge2=
							mdb.models['1-1-1-1'].parts['gusset'].edges.getClosest(coordinates=((gp2[0], gp2[1], 1),))[0][0], plane=
							mdb.models['1-1-1-1'].parts['gusset'].faces.getClosest(coordinates=((gp2[0], gp2[1], 0),))[0][0], planeSide=SIDE1)
						
						mdb.models['1-1-1-1'].parts['gusset'].HoleBlindFromEdges(
							depth=1.0,
							diameter=d_washer,
							distance1=bolts_z1[o],
							distance2=bolts_w,
							edge1=mdb.models['1-1-1-1'].parts['gusset'].edges.getClosest(coordinates=((gp3[0]/2, gp3[1]/2, 0),))[0][0],
							edge2=mdb.models['1-1-1-1'].parts['gusset'].edges.getClosest(coordinates=((gp3[0], gp3[1], 1),))[0][0],
							plane=mdb.models['1-1-1-1'].parts['gusset'].faces.getClosest(coordinates=((gp3[0], gp3[1], 0),))[0][0],
							planeSide=SIDE1
							)
					
					# Partition gusset
					
					mdb.models['1-1-1-1'].parts['gusset'].DatumPointByCoordinate((gp1[0]-current_llip*cos(5*pi/6), gp1[1]-current_llip*sin(5*pi/6), 0),)
					mdb.models['1-1-1-1'].parts['gusset'].DatumPointByCoordinate((gp1[0]-current_llip*cos(5*pi/6), gp1[1]-current_llip*sin(5*pi/6), current_d),)
					
					mdb.models['1-1-1-1'].parts['gusset'].PartitionFaceByShortestPath(
						faces=mdb.models['1-1-1-1'].parts['gusset'].faces.getClosest(coordinates=((gp1[0], gp1[1], 0),))[0][0],
						point1=mdb.models['1-1-1-1'].parts['gusset'].datum.items()[0][1],
						point2=mdb.models['1-1-1-1'].parts['gusset'].datum.items()[1][1],
						)
					
					mdb.models['1-1-1-1'].parts['gusset'].DatumPointByCoordinate((gp2[0]-current_llip*cos(-pi/2), gp2[1]-current_llip*sin(-pi/2), 0),)
					mdb.models['1-1-1-1'].parts['gusset'].DatumPointByCoordinate((gp2[0]-current_llip*cos(-pi/2), gp2[1]-current_llip*sin(-pi/2), current_d),)
					
					mdb.models['1-1-1-1'].parts['gusset'].PartitionFaceByShortestPath(
						faces=mdb.models['1-1-1-1'].parts['gusset'].faces.getClosest(coordinates=((gp2[0], gp2[1], 0),))[0][0],
						point1=mdb.models['1-1-1-1'].parts['gusset'].datum.items()[2][1],
						point2=mdb.models['1-1-1-1'].parts['gusset'].datum.items()[3][1],
						)
					
					mdb.models['1-1-1-1'].parts['gusset'].DatumPointByCoordinate((gp3[0]-current_llip*cos(pi/6), gp3[1]-current_llip*sin(pi/6), 0),)
					mdb.models['1-1-1-1'].parts['gusset'].DatumPointByCoordinate((gp3[0]-current_llip*cos(pi/6), gp3[1]-current_llip*sin(pi/6), current_d),)
					
					mdb.models['1-1-1-1'].parts['gusset'].PartitionFaceByShortestPath(
						faces=mdb.models['1-1-1-1'].parts['gusset'].faces.getClosest(coordinates=((gp3[0], gp3[1], 0),))[0][0],
						point1=mdb.models['1-1-1-1'].parts['gusset'].datum.items()[4][1],
						point2=mdb.models['1-1-1-1'].parts['gusset'].datum.items()[5][1],
						)
					
					# Assign sections
					
					# -for sector
					mdb.models[current_model].parts['sector'].Set(faces=
						mdb.models[current_model].parts['sector'].faces[:], name='AllSectorFaces')
					mdb.models[current_model].parts['sector'].SectionAssignment(offset=0.0, offsetField=''
						, offsetType=MIDDLE_SURFACE, region=
						mdb.models[current_model].parts['sector'].sets['AllSectorFaces'], sectionName='sector', 
						thicknessAssignment=FROM_SECTION)
		
					# -for gusset
					mdb.models[current_model].parts['gusset'].Set(faces=
						mdb.models[current_model].parts['gusset'].faces[:]
						, name='AllGussetFaces')
					mdb.models[current_model].parts['gusset'].SectionAssignment(offset=0.0, offsetField=''
						, offsetType=MIDDLE_SURFACE, region=
						mdb.models[current_model].parts['gusset'].sets['AllGussetFaces'], sectionName='gusset', 
						thicknessAssignment=FROM_SECTION)
					
					# Meshing
					
					# Global seeding in mm
					seedsize = 30
					
					# -Sector
					mdb.models['1-1-1-1'].parts['sector'].setMeshControls(
						algorithm=MEDIAL_AXIS,
						elemShape=QUAD, 
						regions=mdb.models['1-1-1-1'].parts['sector'].faces[:])
					mdb.models['1-1-1-1'].parts['sector'].seedPart(deviationFactor=0.1, 
						minSizeFactor=0.1, size=seedsize)
					mdb.models['1-1-1-1'].parts['sector'].generateMesh()
					
					# -Gusset
					mdb.models['1-1-1-1'].parts['gusset'].setMeshControls(
						algorithm=MEDIAL_AXIS,
						elemShape=QUAD,
						regions=mdb.models['1-1-1-1'].parts['gusset'].faces[:])
					mdb.models['1-1-1-1'].parts['gusset'].seedPart(deviationFactor=0.1, 
						minSizeFactor=0.1, size=seedsize)
					mdb.models['1-1-1-1'].parts['gusset'].generateMesh()
					
					# Create assembly
					
					mdb.models[current_model].rootAssembly.DatumCsysByDefault(CARTESIAN)
					
					# -Sectors
					mdb.models[current_model].rootAssembly.Instance(dependent=ON, name='sector-1', part=
						mdb.models[current_model].parts['sector'])
					mdb.models[current_model].rootAssembly.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
					mdb.models[current_model].rootAssembly.RadialInstancePattern(axis=(0.0, 0.0, 1.0), 
						instanceList=('sector-1', ), number=3, point=(0.0, 0.0, 0.0), totalAngle=
						240.0)
		
					# -Gusset plate
					
					# --Create the instances
					mdb.models[current_model].rootAssembly.Instance(dependent=ON, name='gusset-1', part=
						mdb.models[current_model].parts['gusset'])
					mdb.models[current_model].rootAssembly.Instance(dependent=ON, name='gusset-2', part=
						mdb.models[current_model].parts['gusset'])
					mdb.models[current_model].rootAssembly.Instance(dependent=ON, name='gusset-3', part=
						mdb.models[current_model].parts['gusset'])
					
					# --Translate them to the right position
					mdb.models[current_model].rootAssembly.instances['gusset-2'].translate(vector=(
						0.0, 0.0, (current_l + current_d)))
					mdb.models[current_model].rootAssembly.instances['gusset-3'].translate(vector=(
						0.0, 0.0, 2*(current_l + current_d)))
					
					# Create reference points for BCs/loads.
					
					# -RPs for the faces at the two ends of the columns
					mdb.models[current_model].rootAssembly.ReferencePoint(point=(0.0, 0.0, 0.0))
					mdb.models[current_model].rootAssembly.ReferencePoint(point=(0.0, 0.0, (2*current_l + 3*current_d)))
					
					# - RP at the middle
					mdb.models[current_model].rootAssembly.ReferencePoint(point=(0.0, 0.0, (current_l + 1.5*current_d)))
					
					# Meshing
					
#					# -Sector
#					mdb.models['1-1-1-1'].parts['sector'].setMeshControls(elemShape=QUAD, regions=
#						mdb.models['1-1-1-1'].parts['sector'].faces[:])
#					mdb.models['1-1-1-1'].parts['sector'].seedPart(deviationFactor=0.1, 
#						minSizeFactor=0.1, size=6.3)
#					mdb.models['1-1-1-1'].parts['sector'].generateMesh()
#					
#					# -Gusset
#					mdb.models['1-1-1-1'].parts['gusset'].setMeshControls(elemShape=QUAD, regions=
#						mdb.models['1-1-1-1'].parts['gusset'].faces[:])
#					mdb.models['1-1-1-1'].parts['gusset'].seedPart(deviationFactor=0.1, 
#						minSizeFactor=0.1, size=6.3)
#					mdb.models['1-1-1-1'].parts['gusset'].generateMesh()
					
					
		
#		
#					# Create buckling step
#					mdb.models[current_model].BuckleStep(maxIterations=300, name='Step-1', numEigen=10, 
#						previous='Initial', vectors=18)
#		
#					# Create face couplings for BCs
#					# -Face 1
#					mdb.models[current_model].rootAssembly.Set(name='m_Set-1', referencePoints=(
#						mdb.models[current_model].rootAssembly.referencePoints[17], ))
#					mdb.models[current_model].rootAssembly.Set(edges=
#						mdb.models[current_model].rootAssembly.instances['sector-1-rad-2'].edges.getSequenceFromMask(
#						mask=('[#49249244 #92492492 #924 ]', ), )+\
#						mdb.models[current_model].rootAssembly.instances['sector-1-rad-2-rad-2'].edges.getSequenceFromMask(
#						mask=('[#49249244 #92492492 #924 ]', ), )+\
#						mdb.models[current_model].rootAssembly.instances['sector-1'].edges.getSequenceFromMask(
#						mask=('[#49249244 #92492492 #924 ]', ), ), name='s_Set-1')
#					mdb.models[current_model].Coupling(controlPoint=
#						mdb.models[current_model].rootAssembly.sets['m_Set-1'], couplingType=KINEMATIC, 
#						influenceRadius=WHOLE_SURFACE, localCsys=None, name='Constraint-1', 
#						surface=mdb.models[current_model].rootAssembly.sets['s_Set-1'], u1=ON, u2=ON, u3=
#						ON, ur1=ON, ur2=ON, ur3=ON)
#						
#					# -Face 2
#					mdb.models[current_model].rootAssembly.Set(name='m_Set-3', referencePoints=(
#						mdb.models[current_model].rootAssembly.referencePoints[18], ))
#					mdb.models[current_model].rootAssembly.Set(edges=
#						mdb.models[current_model].rootAssembly.instances['sector-1-rad-2'].edges.getSequenceFromMask(
#						mask=('[#92492491 #24924924 #249 ]', ), )+\
#						mdb.models[current_model].rootAssembly.instances['sector-1'].edges.getSequenceFromMask(
#						mask=('[#92492491 #24924924 #249 ]', ), )+\
#						mdb.models[current_model].rootAssembly.instances['sector-1-rad-2-rad-2'].edges.getSequenceFromMask(
#						mask=('[#92492491 #24924924 #249 ]', ), ), name='s_Set-3')
#					mdb.models[current_model].Coupling(controlPoint=
#						mdb.models[current_model].rootAssembly.sets['m_Set-3'], couplingType=KINEMATIC, 
#						influenceRadius=WHOLE_SURFACE, localCsys=None, name='Constraint-2', 
#						surface=mdb.models[current_model].rootAssembly.sets['s_Set-3'], u1=ON, u2=ON, u3=
#						ON, ur1=ON, ur2=ON, ur3=ON)
#						
#					# Fasteners
#					# -Create datum points
#					mdb.models[current_model].rootAssembly.DatumPointByOffset(point=
#						mdb.models[current_model].rootAssembly.instances['sector-1'].InterestingPoint(
#						mdb.models[current_model].rootAssembly.instances['sector-1'].edges[2], MIDDLE), 
#						vector=(0.0, 0.0, 50.0))
#					# Boundary conditions
#					mdb.models[current_model].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
#						distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1', 
#						region=mdb.models[current_model].rootAssembly.sets['m_Set-1'], u1=SET, u2=SET, u3=
#						UNSET, ur1=UNSET, ur2=UNSET, ur3=SET)
#					mdb.models[current_model].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
#						distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-2', 
#						region=mdb.models[current_model].rootAssembly.sets['m_Set-3'], u1=SET, u2=SET, u3=
#						SET, ur1=UNSET, ur2=UNSET, ur3=SET)
#						
#					# Apply load
#					mdb.models[current_model].ConcentratedForce(cf3=1000.0, createStepName='Step-1', 
#						distributionType=UNIFORM, field='', localCsys=None, name='Load-1', region=
#						mdb.models[current_model].rootAssembly.sets['m_Set-1'])
#		
#					# Meshing
#					mdb.models[current_model].parts['gusset'].seedPart(deviationFactor=0.1, minSizeFactor=
#						0.1, size=25.0)
#					mdb.models[current_model].parts['sector'].seedPart(deviationFactor=0.1, minSizeFactor=
#						0.1, size=25.0)
#					mdb.models[current_model].parts['sector'].generateMesh()
#					mdb.models[current_model].parts['gusset'].generateMesh()
#					mdb.models[current_model].rootAssembly.regenerate()
#					
##					# Modify keyword for nodefile
##					mdb.models[current_model].keywordBlock.synchVersions(storeNodesAndElements=False)
##					mdb.models[current_model].keywordBlock.replace(101, '\n*Output, field, variable=PRESELECT\n*NODEFILE\nU')
#					
##					2nd Phase: Convert model to riks analysis
#				
#					riks_model = 'RIKS-'+str(i+1)+'-'+str(j+1)+'-'+str(k+1)
#			
#					# copy model from buckling analysis
#					mdb.Model(name=riks_model, objectToCopy=mdb.models[current_model])
#			
##					# Delete keyword nodefile
##					mdb.models[riks_model].keywordBlock.synchVersions(storeNodesAndElements=False)
##					mdb.models[riks_model].keywordBlock.replace(102, '\n')
#		
#					# Change material model, plasticity added (more points needed in plasticity table)
#					mdb.models[riks_model].materials['pure-elastic'].Plastic(table=((355.0, 0.0), ))
#					mdb.models[riks_model].StaticRiksStep(maintainAttributes=True, name='Step-1', 
#					nlgeom=ON, previous='Initial')
#			
#					# Load supressed
#					mdb.models[riks_model].loads['Load-1'].suppress()
#		
#					# Change boundary conditions(recheck/delete the first part)
#					mdb.models[riks_model].DisplacementBC(amplitude=UNSET, createStepName='Step-1'
#					, distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
#					'BC-3', region=mdb.models[riks_model].rootAssembly.sets['m_Set-1'], u1=
#					UNSET, u2=UNSET, u3=1.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
#					mdb.models[riks_model].boundaryConditions['BC-1'].setValuesInStep(stepName=
#					'Step-1', u3=1.0)
#					del mdb.models[riks_model].boundaryConditions['BC-3']
#		
##					# Change keywords to include initial imperfections file (filename was given wrong initially and corrected later)
##					mdb.models[riks_model].keywordBlock.synchVersions(storeNodesAndElements=False)
##					mdb.models[riks_model].keywordBlock.replace(88, 
##					'\n** ----------------------------------------------------------------\n** \n**********GEOMETRICAL IMPERFECTIONS\n*IMPERFECTION,FILE=current_model,STEP=1\n1,4\n\n** STEP: Step-1\n**')

# Delete initial model
del mdb.models['Model-1']

