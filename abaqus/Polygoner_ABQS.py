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
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
# Import profiles database and profile metadata --------------------------------------------------------------

# Import pickle to load the .pkl database
import pickle

# Open and read the database
profiles_file = open("profiles.pkl",'rb')
profiles = pickle.load(profiles_file)
profiles_file.close()

profiles_file = open("meta.pkl",'rb')
profiles_meta = pickle.load(profiles_file)
profiles_file.close()


#for i in range(profiles.shape[0]):
#	for j in range(profiles.shape[1]):
#		for k in range(profiles.shape[2]):
#			for l in range(profiles_meta.shape[3]):
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
				
				# Create model ----------------------------------------------------------------------------------------------------------------------
				mdb.Model(modelType=STANDARD_EXPLICIT, name=current_model)
				c_model = mdb.models[current_model]
				
				# Create Parts ----------------------------------------------------------------------------------------------------------------------
				
				# Sector
				
				# Loop through the different bolt spacings (temporary b to change)
				b = [0.5, 1, 1.5]
				for m in range(1):
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
#					for o in range(2):
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
					
					c_model.Material(name='pure-elastic')
					c_model.materials['pure-elastic'].Elastic(table=((210000.0, 0.3), ))
		
					# Create sections ---------------------------------------------------------------------------------------------------------------
					
					# -for sector
					c_model.HomogeneousShellSection(
						idealization=NO_IDEALIZATION, 
						integrationRule=SIMPSON,
						material='pure-elastic',
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
						material='pure-elastic',
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
						
					# Step -----------------------------------------------------------------------------------
					
					load_step = c_model.StaticStep(
						name='Load', 
						previous='Initial'
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
					# Modify the BC end-2 and apply a displacement
					
					c_model.boundaryConditions['fix-end2'].setValuesInStep(
						stepName='Load',
						u3=-5.0)
					
					# Create the job
					
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
						model=current_model,
						modelPrint=OFF, 
					    multiprocessingMode=DEFAULT,
						name='Job-'+current_model,
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

# Delete initial model
del mdb.models['Model-1']

