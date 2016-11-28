# -*- coding: mbcs -*-
import pickle
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
for i in range(1):
	for j in range(1):
		for k in range(1):
		    # Variables holding information of the current profile
			current_model = str(i+1)+'-'+str(j+1)+'-'+str(k+1)
			current_d = float(profiles_meta[i][j][k][0][0])
			current_t = float(profiles_meta[i][j][k][1][0])
			current_tg = float(profiles_meta[i][j][k][2][0])
			current_fy = float(profiles_meta[i][j][k][3][0])
			
			# Create model
			mdb.Model(modelType=STANDARD_EXPLICIT, name=current_model)

			# Create Parts
			# -Profile sketch for sector
			mdb.models[current_model].ConstrainedSketch(name='__profile__', sheetSize=1200.0)
			
			# -Sketch sector lines
			for n in range(profiles[i][j][k].shape[1]-1):
				mdb.models[current_model].sketches['__profile__'].Line(point1=(profiles[i][j][k][0][n], profiles[i][j][k][1][n]), 
					point2=(profiles[i][j][k][0][n+1], profiles[i][j][k][1][n+1]))
			
			# -Extrude sector part
			mdb.models[current_model].Part(dimensionality=THREE_D, name='sector', type=
				DEFORMABLE_BODY)
			mdb.models[current_model].parts['sector'].BaseShellExtrude(depth=4000.0, sketch=
				mdb.models[current_model].sketches['__profile__'])
			del mdb.models[current_model].sketches['__profile__']
			
			# -Profile sketch for gusset
			mdb.models[current_model].ConstrainedSketch(name='__profile__', sheetSize=1200.0)
			
			# -Sketch gusset lines
			mdb.models[current_model].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(0.0, 
				-180.0))
			mdb.models[current_model].sketches['__profile__'].VerticalConstraint(addUndoState=
				False, entity=mdb.models[current_model].sketches['__profile__'].geometry[2])
			mdb.models[current_model].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
				155.8845727, 90.0))
			mdb.models[current_model].sketches['__profile__'].Line(point1=(0.0, 0.0), point2=(
				-155.8845727, 90.0))

			# -Extrude gusset part
			mdb.models[current_model].Part(dimensionality=THREE_D, name='gusset', type=
				DEFORMABLE_BODY)
			mdb.models[current_model].parts['gusset'].BaseShellExtrude(depth=300.0, sketch=
				mdb.models[current_model].sketches['__profile__'])
			del mdb.models[current_model].sketches['__profile__']

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

			# Assign sections
			# -for sector
			mdb.models[current_model].parts['sector'].Set(faces=
				mdb.models[current_model].parts['sector'].faces.getSequenceFromMask((
				'[#1ffffff ]', ), ), name='Set-1')
			mdb.models[current_model].parts['sector'].SectionAssignment(offset=0.0, offsetField=''
				, offsetType=MIDDLE_SURFACE, region=
				mdb.models[current_model].parts['sector'].sets['Set-1'], sectionName='sector', 
				thicknessAssignment=FROM_SECTION)

			# -for gusset
			mdb.models[current_model].parts['gusset'].Set(faces=
				mdb.models[current_model].parts['gusset'].faces.getSequenceFromMask(('[#7 ]', ), )
				, name='Set-1')
			mdb.models[current_model].parts['gusset'].SectionAssignment(offset=0.0, offsetField=''
				, offsetType=MIDDLE_SURFACE, region=
				mdb.models[current_model].parts['gusset'].sets['Set-1'], sectionName='gusset', 
				thicknessAssignment=FROM_SECTION)

			# Create assembly
			mdb.models[current_model].rootAssembly.DatumCsysByDefault(CARTESIAN)
			# -Sectors (make dependent=ON)
			mdb.models[current_model].rootAssembly.Instance(dependent=ON, name='sector-1', part=
				mdb.models[current_model].parts['sector'])
			mdb.models[current_model].rootAssembly.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
			mdb.models[current_model].rootAssembly.RadialInstancePattern(axis=(0.0, 0.0, 1.0), 
				instanceList=('sector-1', ), number=2, point=(0.0, 0.0, 0.0), totalAngle=
				120.0)
			mdb.models[current_model].rootAssembly.RadialInstancePattern(axis=(0.0, 0.0, 1.0), 
				instanceList=('sector-1-rad-2', ), number=2, point=(0.0, 0.0, 0.0), 
				totalAngle=120.0)

			# -Gusset plate (Translation to be as constraint)
			mdb.models[current_model].rootAssembly.Instance(dependent=ON, name='gusset-1', part=
				mdb.models[current_model].parts['gusset'])
			mdb.models[current_model].rootAssembly.Instance(dependent=ON, name='gusset-2', part=
				mdb.models[current_model].parts['gusset'])
			mdb.models[current_model].rootAssembly.instances['gusset-2'].translate(vector=(
				385.995645402491, 0.0, 0.0))
			mdb.models[current_model].rootAssembly.Instance(dependent=ON, name='gusset-3', part=
				mdb.models[current_model].parts['gusset'])
			mdb.models[current_model].rootAssembly.instances['gusset-3'].translate(vector=(
				728.941705342491, 0.0, 0.0))
			mdb.models[current_model].rootAssembly.DatumPointByCoordinate(coords=(0.0, 0.0, 
				2000.0))
			mdb.models[current_model].rootAssembly.CoincidentPoint(fixedPoint=
				mdb.models[current_model].rootAssembly.datums[15], movablePoint=
				mdb.models[current_model].rootAssembly.instances['gusset-2'].InterestingPoint(
				mdb.models[current_model].rootAssembly.instances['gusset-2'].edges[1], MIDDLE))
			mdb.models[current_model].rootAssembly.ReferencePoint(point=(0.0, 0.0, 0.0))
			mdb.models[current_model].rootAssembly.ReferencePoint(point=(0.0, 0.0, 4000.0))
			mdb.models[current_model].rootAssembly.CoincidentPoint(fixedPoint=
				mdb.models[current_model].rootAssembly.referencePoints[18], movablePoint=
				mdb.models[current_model].rootAssembly.instances['gusset-3'].vertices[1])

			# Create buckling step
			mdb.models[current_model].BuckleStep(maxIterations=300, name='Step-1', numEigen=10, 
				previous='Initial', vectors=18)

			# Create face couplings for BCs
			# -Face 1
			mdb.models[current_model].rootAssembly.Set(name='m_Set-1', referencePoints=(
				mdb.models[current_model].rootAssembly.referencePoints[17], ))
			mdb.models[current_model].rootAssembly.Set(edges=
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2'].edges.getSequenceFromMask(
				mask=('[#49249244 #92492492 #924 ]', ), )+\
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2-rad-2'].edges.getSequenceFromMask(
				mask=('[#49249244 #92492492 #924 ]', ), )+\
				mdb.models[current_model].rootAssembly.instances['sector-1'].edges.getSequenceFromMask(
				mask=('[#49249244 #92492492 #924 ]', ), ), name='s_Set-1')
			mdb.models[current_model].Coupling(controlPoint=
				mdb.models[current_model].rootAssembly.sets['m_Set-1'], couplingType=KINEMATIC, 
				influenceRadius=WHOLE_SURFACE, localCsys=None, name='Constraint-1', 
				surface=mdb.models[current_model].rootAssembly.sets['s_Set-1'], u1=ON, u2=ON, u3=
				ON, ur1=ON, ur2=ON, ur3=ON)
				
			# -Face 2
			mdb.models[current_model].rootAssembly.Set(name='m_Set-3', referencePoints=(
				mdb.models[current_model].rootAssembly.referencePoints[18], ))
			mdb.models[current_model].rootAssembly.Set(edges=
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2'].edges.getSequenceFromMask(
				mask=('[#92492491 #24924924 #249 ]', ), )+\
				mdb.models[current_model].rootAssembly.instances['sector-1'].edges.getSequenceFromMask(
				mask=('[#92492491 #24924924 #249 ]', ), )+\
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2-rad-2'].edges.getSequenceFromMask(
				mask=('[#92492491 #24924924 #249 ]', ), ), name='s_Set-3')
			mdb.models[current_model].Coupling(controlPoint=
				mdb.models[current_model].rootAssembly.sets['m_Set-3'], couplingType=KINEMATIC, 
				influenceRadius=WHOLE_SURFACE, localCsys=None, name='Constraint-2', 
				surface=mdb.models[current_model].rootAssembly.sets['s_Set-3'], u1=ON, u2=ON, u3=
				ON, ur1=ON, ur2=ON, ur3=ON)
				
			# Fasteners
			# -Create datum points
			mdb.models[current_model].rootAssembly.DatumPointByOffset(point=
				mdb.models[current_model].rootAssembly.instances['sector-1'].InterestingPoint(
				mdb.models[current_model].rootAssembly.instances['sector-1'].edges[2], MIDDLE), 
				vector=(0.0, 0.0, 50.0))
			# Boundary conditions
			mdb.models[current_model].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
				distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1', 
				region=mdb.models[current_model].rootAssembly.sets['m_Set-1'], u1=SET, u2=SET, u3=
				UNSET, ur1=UNSET, ur2=UNSET, ur3=SET)
			mdb.models[current_model].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
				distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-2', 
				region=mdb.models[current_model].rootAssembly.sets['m_Set-3'], u1=SET, u2=SET, u3=
				SET, ur1=UNSET, ur2=UNSET, ur3=SET)
				
			# Apply load
			mdb.models[current_model].ConcentratedForce(cf3=1000.0, createStepName='Step-1', 
				distributionType=UNIFORM, field='', localCsys=None, name='Load-1', region=
				mdb.models[current_model].rootAssembly.sets['m_Set-1'])

			# Meshing
			mdb.models[current_model].parts['gusset'].seedPart(deviationFactor=0.1, minSizeFactor=
				0.1, size=25.0)
			mdb.models[current_model].parts['sector'].seedPart(deviationFactor=0.1, minSizeFactor=
				0.1, size=25.0)
			mdb.models[current_model].parts['sector'].generateMesh()
			mdb.models[current_model].parts['gusset'].generateMesh()
			mdb.models[current_model].rootAssembly.regenerate()
			
#			# Modify keyword for nodefile
#			mdb.models[current_model].keywordBlock.synchVersions(storeNodesAndElements=False)
#			mdb.models[current_model].keywordBlock.replace(101, '\n*Output, field, variable=PRESELECT\n*NODEFILE\nU')
			
		# ## 2nd Phase: Convert model to riks analysis
		
			riks_model = 'RIKS-'+str(i+1)+'-'+str(j+1)+'-'+str(k+1)
	
			# copy model from buckling analysis
			mdb.Model(name=riks_model, objectToCopy=mdb.models[current_model])
	
#			# Delete keyword nodefile
#			mdb.models[riks_model].keywordBlock.synchVersions(storeNodesAndElements=False)
#			mdb.models[riks_model].keywordBlock.replace(102, '\n')

			# Change material model, plasticity added (more points needed in plasticity table)
			mdb.models[riks_model].materials['pure-elastic'].Plastic(table=((355.0, 0.0), ))
			mdb.models[riks_model].StaticRiksStep(maintainAttributes=True, name='Step-1', 
			nlgeom=ON, previous='Initial')
    
			# Load supressed
			mdb.models[riks_model].loads['Load-1'].suppress()

			# Change boundary conditions(recheck/delete the first part)
			mdb.models[riks_model].DisplacementBC(amplitude=UNSET, createStepName='Step-1'
			, distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
			'BC-3', region=mdb.models[riks_model].rootAssembly.sets['m_Set-1'], u1=
			UNSET, u2=UNSET, u3=1.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
			mdb.models[riks_model].boundaryConditions['BC-1'].setValuesInStep(stepName=
			'Step-1', u3=1.0)
			del mdb.models[riks_model].boundaryConditions['BC-3']
   
#			# Change keywords to include initial imperfections file (filename was given wrong initially and corrected later)
#			mdb.models[riks_model].keywordBlock.synchVersions(storeNodesAndElements=False)
#			mdb.models[riks_model].keywordBlock.replace(88, 
#			'\n** ----------------------------------------------------------------\n** \n**********GEOMETRICAL IMPERFECTIONS\n*IMPERFECTION,FILE=current_model,STEP=1\n1,4\n\n** STEP: Step-1\n**')

