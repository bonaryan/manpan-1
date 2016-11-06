# Import profiles database
import pickle
profiles_file = open("profiles.pkl",'rb')
profiles = pickle.load(profiles_file)
profiles_file.close()

# -*- coding: mbcs -*-
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

## 1st Phase: Buckling analysis

for i in range(profiles.shape[0]):
	for j in range(profiles.shape[1]):
		for k in range(profiles.shape[2]):
			current_model = str(i+1)+'-'+str(j+1)+'-'+str(k+1)
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
				thickness=5.0, thicknessField='', thicknessModulus=None, thicknessType=
				UNIFORM, useDensity=OFF)
			# -for gusset
			mdb.models[current_model].HomogeneousShellSection(idealization=NO_IDEALIZATION, 
				integrationRule=SIMPSON, material='pure-elastic', name='gusset', numIntPts=
				5, poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
				thickness=10.0, thicknessField='', thicknessModulus=None, thicknessType=
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
			mdb.models[current_model].rootAssembly.Instance(dependent=OFF, name='sector-1', part=
				mdb.models[current_model].parts['sector'])
			mdb.models[current_model].rootAssembly.DatumAxisByPrincipalAxis(principalAxis=ZAXIS)
			mdb.models[current_model].rootAssembly.RadialInstancePattern(axis=(0.0, 0.0, 1.0), 
				instanceList=('sector-1', ), number=2, point=(0.0, 0.0, 0.0), totalAngle=
				120.0)
			mdb.models[current_model].rootAssembly.RadialInstancePattern(axis=(0.0, 0.0, 1.0), 
				instanceList=('sector-1-rad-2', ), number=2, point=(0.0, 0.0, 0.0), 
				totalAngle=120.0)

			# -Gusset plate (Translation to be as constraint)
			mdb.models[current_model].rootAssembly.Instance(dependent=OFF, name='gusset-1', part=
				mdb.models[current_model].parts['gusset'])
			mdb.models[current_model].rootAssembly.Instance(dependent=OFF, name='gusset-2', part=
				mdb.models[current_model].parts['gusset'])
			mdb.models[current_model].rootAssembly.instances['gusset-2'].translate(vector=(
				385.995645402491, 0.0, 0.0))
			mdb.models[current_model].rootAssembly.Instance(dependent=OFF, name='gusset-3', part=
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
				vector=(0.0, 0.0, -50.0))
			del mdb.models[current_model].rootAssembly.features['Datum pt-2']
			mdb.models[current_model].rootAssembly.DatumPointByOffset(point=
				mdb.models[current_model].rootAssembly.instances['sector-1'].InterestingPoint(
				mdb.models[current_model].rootAssembly.instances['sector-1'].edges[2], MIDDLE), 
				vector=(0.0, 0.0, 50.0))
			mdb.models[current_model].rootAssembly.DatumPointByOffset(point=
				mdb.models[current_model].rootAssembly.instances['sector-1'].InterestingPoint(
				mdb.models[current_model].rootAssembly.instances['sector-1'].edges[2], MIDDLE), 
				vector=(0.0, 0.0, 3950.0))
			mdb.models[current_model].rootAssembly.DatumPointByOffset(point=
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2-rad-2'].InterestingPoint(
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2-rad-2'].edges[2], 
				MIDDLE), vector=(0.0, 0.0, 50.0))
			mdb.models[current_model].rootAssembly.DatumPointByOffset(point=
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2-rad-2'].InterestingPoint(
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2-rad-2'].edges[2], 
				MIDDLE), vector=(0.0, 0.0, 3950.0))
			mdb.models[current_model].rootAssembly.DatumPointByOffset(point=
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2'].InterestingPoint(
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2'].edges[2], 
				MIDDLE), vector=(0.0, 0.0, 50.0))
			mdb.models[current_model].rootAssembly.DatumPointByOffset(point=
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2'].InterestingPoint(
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2'].edges[2], 
				MIDDLE), vector=(0.0, 0.0, 3950.0))
			mdb.models[current_model].rootAssembly.AttachmentPointsAlongDirection(endPoint=
				mdb.models[current_model].rootAssembly.datums[30], name='Attachment Points-1', 
				pointCreationMethod=AUTO_FIT, setName='Attachment Points-1-Set-1', spacing=
				100.0, startPoint=mdb.models[current_model].rootAssembly.datums[29])
			mdb.models[current_model].rootAssembly.AttachmentPointsAlongDirection(endPoint=
				mdb.models[current_model].rootAssembly.datums[28], name='Attachment Points-2', 
				pointCreationMethod=AUTO_FIT, setName='Attachment Points-2-Set-1', spacing=
				100.0, startPoint=mdb.models[current_model].rootAssembly.datums[27])
			mdb.models[current_model].rootAssembly.AttachmentPointsAlongDirection(endPoint=
				mdb.models[current_model].rootAssembly.datums[26], name='Attachment Points-3', 
				pointCreationMethod=AUTO_FIT, setName='Attachment Points-3-Set-1', spacing=
				100.0, startPoint=mdb.models[current_model].rootAssembly.datums[25])

			# -Create connector section
			mdb.models[current_model].ConnectorSection(assembledType=BEAM, name='ConnSect-1')
			mdb.models[current_model].rootAssembly.engineeringFeatures.PointFastener(name=
				'Fasteners-1', physicalRadius=8.0, region=
				mdb.models[current_model].rootAssembly.sets['Attachment Points-1-Set-1'], 
				sectionName='ConnSect-1')
			mdb.models[current_model].rootAssembly.engineeringFeatures.PointFastener(name=
				'Fasteners-2', physicalRadius=8.0, region=
				mdb.models[current_model].rootAssembly.sets['Attachment Points-2-Set-1'], 
				sectionName='ConnSect-1')
			mdb.models[current_model].rootAssembly.engineeringFeatures.PointFastener(name=
				'Fasteners-3', physicalRadius=8.0, region=
				mdb.models[current_model].rootAssembly.sets['Attachment Points-3-Set-1'], 
				sectionName='ConnSect-1')

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
			mdb.models[current_model].rootAssembly.seedPartInstance(deviationFactor=0.1, 
				minSizeFactor=0.1, regions=(
				mdb.models[current_model].rootAssembly.instances['sector-1'], 
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2'], 
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2-rad-2'], 
				mdb.models[current_model].rootAssembly.instances['gusset-1'], 
				mdb.models[current_model].rootAssembly.instances['gusset-2'], 
				mdb.models[current_model].rootAssembly.instances['gusset-3']), size=25.0)
			mdb.models[current_model].rootAssembly.generateMesh(regions=(
				mdb.models[current_model].rootAssembly.instances['sector-1'], 
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2'], 
				mdb.models[current_model].rootAssembly.instances['sector-1-rad-2-rad-2'], 
				mdb.models[current_model].rootAssembly.instances['gusset-1'], 
				mdb.models[current_model].rootAssembly.instances['gusset-2'], 
				mdb.models[current_model].rootAssembly.instances['gusset-3']))

# # Create jobs
# mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    # explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    # memory=90, memoryUnits=PERCENTAGE, model=current_model, modelPrint=OFF, 
    # multiprocessingMode=DEFAULT, name='poly-buckle', nodalOutputPrecision=
    # SINGLE, numCpus=4, numDomains=4, numGPUs=0, queue=None, resultsFormat=ODB, 
    # scratch='', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
# mdb.models[current_model].keywordBlock.synchVersions(storeNodesAndElements=False)
# mdb.models[current_model].keywordBlock.replace(125, 
    # '\n*Output, field, variable=PRESELECT\n*NODEFILE\nU')

# # Submit buckling for analysis to solver
# mdb.jobs['poly-buckle'].submit(consistencyChecking=OFF)

# Message returns from solver
# mdb.jobs['poly-buckle']._Message(STARTED, {'phase': BATCHPRE_PHASE, 
    # 'clientHost': 'Lenovo-PC', 'handle': 0, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(WARNING, {'phase': BATCHPRE_PHASE, 
    # 'message': 'BUCKLE OPTION IS NOT SUPPORTED FOR ELEMENT LOOP PARALLELIZATION. IF YOU HAVE SPECIFIED ELEMENT LOOP PARALLELIZATION, IT WILL BE TURNED OFF FOR THIS ANALYSIS.', 
    # 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(WARNING, {'phase': BATCHPRE_PHASE, 
    # 'message': 'WHENEVER A TRANSLATION (ROTATION) DOF AT A NODE IS CONSTRAINED BY A KINEMATIC COUPLING DEFINITION THE TRANSLATION (ROTATION) DOFS FOR THAT NODE CANNOT BE INCLUDED IN ANY OTHER CONSTRAINT INCLUDING MPCS, RIGID BODIES, ETC.', 
    # 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(WARNING, {'phase': BATCHPRE_PHASE, 
    # 'message': 'MPCS (EXTERNAL or INTERNAL, including those generated from rigid body definitions), KINEMATIC COUPLINGS, AND/OR EQUATIONS WILL ACTIVATE ADDITIONAL DEGREES OF FREEDOM', 
    # 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(ODB_FILE, {'phase': BATCHPRE_PHASE, 
    # 'file': 'C:\\Temp\\poly-buckle.odb', 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(COMPLETED, {'phase': BATCHPRE_PHASE, 
    # 'message': 'Analysis phase complete', 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(STARTED, {'phase': STANDARD_PHASE, 
    # 'clientHost': 'Lenovo-PC', 'handle': 176, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(STEP, {'phase': STANDARD_PHASE, 'stepId': 1, 
    # 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 
    # 'step': 0, 'frame': 0, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(STATUS, {'totalTime': 0.0, 'attempts': 1, 
    # 'timeIncrement': 1e-36, 'increment': 0, 'stepTime': 0.0, 'step': 1, 
    # 'jobName': 'poly-buckle', 'severe': 0, 'iterations': 0, 
    # 'phase': STANDARD_PHASE, 'equilibrium': 0})
# mdb.jobs['poly-buckle']._Message(MEMORY_ESTIMATE, {'phase': STANDARD_PHASE, 
    # 'jobName': 'poly-buckle', 'memory': 383.0})
# mdb.jobs['poly-buckle']._Message(MEMORY_ESTIMATE, {'phase': STANDARD_PHASE, 
    # 'jobName': 'poly-buckle', 'memory': 382.0})
# mdb.jobs['poly-buckle']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 
    # 'step': 0, 'frame': 1, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 
    # 'step': 0, 'frame': 2, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 
    # 'step': 0, 'frame': 3, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 
    # 'step': 0, 'frame': 4, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 
    # 'step': 0, 'frame': 5, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 
    # 'step': 0, 'frame': 6, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 
    # 'step': 0, 'frame': 7, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 
    # 'step': 0, 'frame': 8, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 
    # 'step': 0, 'frame': 9, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 
    # 'step': 0, 'frame': 10, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(STATUS, {'totalTime': 0.0, 'attempts': 1, 
    # 'timeIncrement': 1e-36, 'increment': 1, 'stepTime': 1e-36, 'step': 1, 
    # 'jobName': 'poly-buckle', 'severe': 0, 'iterations': 0, 
    # 'phase': STANDARD_PHASE, 'equilibrium': 0})
# mdb.jobs['poly-buckle']._Message(END_STEP, {'phase': STANDARD_PHASE, 
    # 'stepId': 1, 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(COMPLETED, {'phase': STANDARD_PHASE, 
    # 'message': 'Analysis phase complete', 'jobName': 'poly-buckle'})
# mdb.jobs['poly-buckle']._Message(JOB_COMPLETED, {
    # 'time': 'Fri Nov 04 23:14:14 2016', 'jobName': 'poly-buckle'})




# ## 2nd Phase: Convert model to riks analysis


# # Change material model, plasticity added (more points needed in plasticity table)
# mdb.models['p111'].materials['pure-elastic'].Plastic(table=((355.0, 0.0), ))
# mdb.models['p111'].keywordBlock.synchVersions(storeNodesAndElements=False)
# mdb.models['p111'].keywordBlock.replace(127, '\n')
# mdb.models['p111'].StaticRiksStep(maintainAttributes=True, name='Step-1', 
    # nlgeom=ON, previous='Initial')

# # Load deleted
# del mdb.models['p111'].loads['Load-1']

# # Change boundary conditions
# mdb.models['p111'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
    # distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    # 'BC-3', region=mdb.models['p111'].rootAssembly.sets['m_Set-1'], u1=UNSET, 
    # u2=UNSET, u3=1.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)

# # Change keywords to include initial imperfections file (filename was given wrong initially and corrected later)
# mdb.models['p111'].keywordBlock.synchVersions(storeNodesAndElements=False)
# mdb.models['p111'].keywordBlock.replace(112, 
    # '\n** ----------------------------------------------------------------\n** \n****************GEOMETRICAL IMPERFECTIONS\n*IMPERFECTION,FILE=polygon-buckle,STEP=1\n1,4\n\n** STEP: Step-1\n**')
# mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    # explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    # memory=90, memoryUnits=PERCENTAGE, model='p111', modelPrint=OFF, 
    # multiprocessingMode=DEFAULT, name='poly-riks', nodalOutputPrecision=SINGLE, 
    # numCpus=4, numDomains=4, numGPUs=0, queue=None, resultsFormat=ODB, scratch=
    # '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
# mdb.jobs['poly-riks'].submit(consistencyChecking=OFF)
# #* The following results file(s) could not be located: polygon-buckle.fil
# mdb.models['p111'].keywordBlock.synchVersions(storeNodesAndElements=False)
# mdb.models['p111'].keywordBlock.replace(113, 
    # '\n*IMPERFECTION,FILE=poly-buckle,STEP=1\n1,4')

# # Create history output request
# mdb.models['p111'].HistoryOutputRequest(createStepName='Step-1', name=
    # 'H-Output-2', rebar=EXCLUDE, region=
    # mdb.models['p111'].rootAssembly.sets['m_Set-1'], sectionPoints=DEFAULT, 
    # variables=('RF3', ))
# mdb.models['p111'].HistoryOutputRequest(createStepName='Step-1', name=
    # 'H-Output-3', rebar=EXCLUDE, region=
    # mdb.models['p111'].rootAssembly.sets['m_Set-1'], sectionPoints=DEFAULT, 
    # variables=('U3', ))

# # Submit riks job to solver
# mdb.jobs['poly-riks'].submit(consistencyChecking=OFF)
# mdb.jobs['poly-riks']._Message(STARTED, {'phase': BATCHPRE_PHASE, 
    # 'clientHost': 'Lenovo-PC', 'handle': 0, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STARTED, {'phase': BATCHPRE_PHASE, 
    # 'clientHost': 'Lenovo-PC', 'handle': 0, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(WARNING, {'phase': BATCHPRE_PHASE, 
    # 'message': 'WHENEVER A TRANSLATION (ROTATION) DOF AT A NODE IS CONSTRAINED BY A KINEMATIC COUPLING DEFINITION THE TRANSLATION (ROTATION) DOFS FOR THAT NODE CANNOT BE INCLUDED IN ANY OTHER CONSTRAINT INCLUDING MPCS, RIGID BODIES, ETC.', 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(WARNING, {'phase': BATCHPRE_PHASE, 
    # 'message': 'WHENEVER A TRANSLATION (ROTATION) DOF AT A NODE IS CONSTRAINED BY A KINEMATIC COUPLING DEFINITION THE TRANSLATION (ROTATION) DOFS FOR THAT NODE CANNOT BE INCLUDED IN ANY OTHER CONSTRAINT INCLUDING MPCS, RIGID BODIES, ETC.', 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(WARNING, {'phase': BATCHPRE_PHASE, 
    # 'message': 'MPCS (EXTERNAL or INTERNAL, including those generated from rigid body definitions), KINEMATIC COUPLINGS, AND/OR EQUATIONS WILL ACTIVATE ADDITIONAL DEGREES OF FREEDOM', 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(WARNING, {'phase': BATCHPRE_PHASE, 
    # 'message': 'MPCS (EXTERNAL or INTERNAL, including those generated from rigid body definitions), KINEMATIC COUPLINGS, AND/OR EQUATIONS WILL ACTIVATE ADDITIONAL DEGREES OF FREEDOM', 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FILE, {'phase': BATCHPRE_PHASE, 
    # 'file': 'C:\\Temp\\poly-riks.odb', 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FILE, {'phase': BATCHPRE_PHASE, 
    # 'file': 'C:\\Temp\\poly-riks.odb', 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(COMPLETED, {'phase': BATCHPRE_PHASE, 
    # 'message': 'Analysis phase complete', 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(COMPLETED, {'phase': BATCHPRE_PHASE, 
    # 'message': 'Analysis phase complete', 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STARTED, {'phase': STANDARD_PHASE, 
    # 'clientHost': 'Lenovo-PC', 'handle': 11416, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STARTED, {'phase': STANDARD_PHASE, 
    # 'clientHost': 'Lenovo-PC', 'handle': 11416, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STEP, {'phase': STANDARD_PHASE, 'stepId': 1, 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STEP, {'phase': STANDARD_PHASE, 'stepId': 1, 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 0, 'incrementLPF': 1.0, 
    # 'increment': 0, 'step': 1, 'jobName': 'poly-riks', 'severe': 0, 
    # 'iterations': 0, 'phase': STANDARD_PHASE, 'lpf': 0.0, 'equilibrium': 0})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 0, 'incrementLPF': 1.0, 
    # 'increment': 0, 'step': 1, 'jobName': 'poly-riks', 'severe': 0, 
    # 'iterations': 0, 'phase': STANDARD_PHASE, 'lpf': 0.0, 'equilibrium': 0})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 0, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 0, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(MEMORY_ESTIMATE, {'phase': STANDARD_PHASE, 
    # 'jobName': 'poly-riks', 'memory': 334.0})
# mdb.jobs['poly-riks']._Message(MEMORY_ESTIMATE, {'phase': STANDARD_PHASE, 
    # 'jobName': 'poly-riks', 'memory': 334.0})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.998706098041639, 'increment': 1, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 0.998706098041639, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.998706098041639, 'increment': 1, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 0.998706098041639, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 1, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 1, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.995703192175198, 'increment': 2, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 1.99440929021684, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.995703192175198, 'increment': 2, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 1.99440929021684, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 2, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 2, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 1.48662427640729, 'increment': 3, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 3.48103356662412, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 1.48662427640729, 'increment': 3, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 3.48103356662412, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 3, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 3, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.1537712558842, 'increment': 4, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 5.63480482250832, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.1537712558842, 'increment': 4, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 5.63480482250832, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 4, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 4, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(WARNING, {'phase': STANDARD_PHASE, 
    # 'message': 'THE STRAIN INCREMENT HAS EXCEEDED FIFTY TIMES THE STRAIN TO CAUSE FIRST YIELD AT 58215 POINTS', 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(WARNING, {'phase': STANDARD_PHASE, 
    # 'message': 'THE STRAIN INCREMENT HAS EXCEEDED FIFTY TIMES THE STRAIN TO CAUSE FIRST YIELD AT 58215 POINTS', 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(WARNING, {'phase': STANDARD_PHASE, 
    # 'message': 'THE STRAIN INCREMENT IS SO LARGE THAT THE PROGRAM WILL NOT ATTEMPT THE PLASTICITY CALCULATION AT 20627 POINTS', 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(WARNING, {'phase': STANDARD_PHASE, 
    # 'message': 'THE STRAIN INCREMENT IS SO LARGE THAT THE PROGRAM WILL NOT ATTEMPT THE PLASTICITY CALCULATION AT 20627 POINTS', 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': ' 1U', 
    # 'incrementLPF': -16.6978880207524, 'increment': 5, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 5.63480482250832, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': ' 1U', 
    # 'incrementLPF': -16.6978880207524, 'increment': 5, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 5.63480482250832, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 2, 
    # 'incrementLPF': 0.342124948346917, 'increment': 5, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 5.97692977085524, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 2, 
    # 'incrementLPF': 0.342124948346917, 'increment': 5, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 5.97692977085524, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 5, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 5, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.19810828577914, 'increment': 6, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.17503805663438, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.19810828577914, 'increment': 6, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.17503805663438, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 6, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 6, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0871364719004648, 'increment': 7, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.26217452853484, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0871364719004648, 'increment': 7, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.26217452853484, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 7, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 7, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0504515643797298, 'increment': 8, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.31262609291457, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0504515643797298, 'increment': 8, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.31262609291457, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 8, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 8, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0410839762790534, 'increment': 9, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 7, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.35371006919362, 'equilibrium': 7})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0410839762790534, 'increment': 9, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 7, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.35371006919362, 'equilibrium': 7})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 9, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 9, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0367130179214757, 'increment': 10, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.3904230871151, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0367130179214757, 'increment': 10, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.3904230871151, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 10, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 10, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0348163403912681, 'increment': 11, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.42523942750637, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0348163403912681, 'increment': 11, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.42523942750637, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 11, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 11, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0335043890486987, 'increment': 12, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.45874381655507, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0335043890486987, 'increment': 12, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.45874381655507, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 12, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 12, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0328667555687136, 'increment': 13, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.49161057212378, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0328667555687136, 'increment': 13, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.49161057212378, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 13, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 13, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0489799379170349, 'increment': 14, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.54059051004081, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0489799379170349, 'increment': 14, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.54059051004081, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 14, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 14, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0744908409368068, 'increment': 15, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.61508135097762, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0744908409368068, 'increment': 15, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.61508135097762, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 15, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 15, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0764106014715805, 'increment': 16, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.6914919524492, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0764106014715805, 'increment': 16, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.6914919524492, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 16, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 16, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0783378585603091, 'increment': 17, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.76982981100951, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.0783378585603091, 'increment': 17, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.76982981100951, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 17, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 17, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.120818181469894, 'increment': 18, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.89064799247941, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.120818181469894, 'increment': 18, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 6.89064799247941, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 18, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 18, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.188098490929076, 'increment': 19, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 7.07874648340848, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.188098490929076, 'increment': 19, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 7.07874648340848, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 19, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 19, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.293289761228218, 'increment': 20, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 7.3720362446367, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.293289761228218, 'increment': 20, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 7.3720362446367, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 20, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 20, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.448338939363209, 'increment': 21, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 7.82037518399991, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.448338939363209, 'increment': 21, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 7.82037518399991, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 21, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 21, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.479487013952232, 'increment': 22, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 8.29986219795214, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.479487013952232, 'increment': 22, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 8.29986219795214, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 22, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 22, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.549880385242431, 'increment': 23, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 8.84974258319457, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.549880385242431, 'increment': 23, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 8.84974258319457, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 23, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 23, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.638089296154483, 'increment': 24, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 9.48783187934906, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.638089296154483, 'increment': 24, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 9.48783187934906, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 24, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 24, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.727980050397557, 'increment': 25, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 10.2158119297466, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.727980050397557, 'increment': 25, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 10.2158119297466, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 25, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 25, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.810526105227293, 'increment': 26, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 11.0263380349739, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.810526105227293, 'increment': 26, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 11.0263380349739, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 26, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 26, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.883238783658779, 'increment': 27, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 11.9095768186327, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 0.883238783658779, 'increment': 27, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 11.9095768186327, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 27, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 27, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 1.39991527712032, 'increment': 28, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 13.309492095753, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 1.39991527712032, 'increment': 28, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 13.309492095753, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 28, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 28, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(WARNING, {'phase': STANDARD_PHASE, 
    # 'message': 'THE STRAIN INCREMENT HAS EXCEEDED FIFTY TIMES THE STRAIN TO CAUSE FIRST YIELD AT 68 POINTS', 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(WARNING, {'phase': STANDARD_PHASE, 
    # 'message': 'THE STRAIN INCREMENT HAS EXCEEDED FIFTY TIMES THE STRAIN TO CAUSE FIRST YIELD AT 68 POINTS', 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': ' 1U', 
    # 'incrementLPF': 2.4299240750218, 'increment': 29, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 13.309492095753, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': ' 1U', 
    # 'incrementLPF': 2.4299240750218, 'increment': 29, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 13.309492095753, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 2, 
    # 'incrementLPF': 0.817318333277877, 'increment': 29, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 14.1268104290309, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 2, 
    # 'incrementLPF': 0.817318333277877, 'increment': 29, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 14.1268104290309, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 29, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 29, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 1.28404498944272, 'increment': 30, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 15.4108554184736, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 1.28404498944272, 'increment': 30, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 15.4108554184736, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 30, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 30, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 1.32720187376532, 'increment': 31, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 16.7380572922389, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 1.32720187376532, 'increment': 31, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 16.7380572922389, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 31, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 31, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 1.39602196396671, 'increment': 32, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 18.1340792562056, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 1.39602196396671, 'increment': 32, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 18.1340792562056, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 32, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 32, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.23714780028839, 'increment': 33, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 20.371227056494, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.23714780028839, 'increment': 33, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 20.371227056494, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 33, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 33, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.39669294905761, 'increment': 34, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 22.7679200055516, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.39669294905761, 'increment': 34, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 22.7679200055516, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 34, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 34, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.55639691232532, 'increment': 35, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 25.324316917877, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.55639691232532, 'increment': 35, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 25.324316917877, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 35, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 35, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.68801279181373, 'increment': 36, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 28.0123297096907, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.68801279181373, 'increment': 36, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 28.0123297096907, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 36, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 36, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.82172784954582, 'increment': 37, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 30.8340575592365, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.82172784954582, 'increment': 37, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 30.8340575592365, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 37, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 37, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.94535033273028, 'increment': 38, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 33.7794078919668, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.94535033273028, 'increment': 38, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 33.7794078919668, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 38, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 38, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 3.0665937407253, 'increment': 39, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 36.8460016326921, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 3.0665937407253, 'increment': 39, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 36.8460016326921, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 39, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 39, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 4.83592776876619, 'increment': 40, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 41.6819294014583, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 4.83592776876619, 'increment': 40, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 41.6819294014583, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 40, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 40, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 5.12882033666434, 'increment': 41, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 46.8107497381226, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 5.12882033666434, 'increment': 41, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 46.8107497381226, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 41, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 41, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 5.42489306659213, 'increment': 42, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 7, 
    # 'phase': STANDARD_PHASE, 'lpf': 52.2356428047148, 'equilibrium': 7})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 5.42489306659213, 'increment': 42, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 7, 
    # 'phase': STANDARD_PHASE, 'lpf': 52.2356428047148, 'equilibrium': 7})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 42, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 42, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 5.71289779518435, 'increment': 43, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 57.9485405998991, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 5.71289779518435, 'increment': 43, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 57.9485405998991, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 43, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 43, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 5.9951032661697, 'increment': 44, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 63.9436438660688, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 5.9951032661697, 'increment': 44, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 63.9436438660688, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 44, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 44, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.27781196246904, 'increment': 45, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 70.2214558285378, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.27781196246904, 'increment': 45, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 70.2214558285378, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 45, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 45, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': ' 1U', 
    # 'incrementLPF': 6.55224708133305, 'increment': 46, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 70.2214558285378, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': ' 1U', 
    # 'incrementLPF': 6.55224708133305, 'increment': 46, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 70.2214558285378, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 2, 
    # 'incrementLPF': 1.6151251489464, 'increment': 46, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 71.8365809774842, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 2, 
    # 'incrementLPF': 1.6151251489464, 'increment': 46, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 71.8365809774842, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 46, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 46, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 1.62995832267568, 'increment': 47, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 73.4665393001599, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 1.62995832267568, 'increment': 47, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 73.4665393001599, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 47, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 47, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.47599250406385, 'increment': 48, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 75.9425318042238, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 2.47599250406385, 'increment': 48, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 3, 
    # 'phase': STANDARD_PHASE, 'lpf': 75.9425318042238, 'equilibrium': 3})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 48, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 48, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 3.78478812112267, 'increment': 49, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 79.7273199253465, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 3.78478812112267, 'increment': 49, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 79.7273199253465, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 49, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 49, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 5.83385550208664, 'increment': 50, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 85.5611754274331, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 5.83385550208664, 'increment': 50, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 85.5611754274331, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 50, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 50, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.01970712189285, 'increment': 51, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 91.5808825493259, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.01970712189285, 'increment': 51, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 91.5808825493259, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 51, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 51, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.20166265330222, 'increment': 52, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 97.7825452026282, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.20166265330222, 'increment': 52, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 97.7825452026282, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 52, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 52, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.37932978631551, 'increment': 53, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 104.161874988944, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.37932978631551, 'increment': 53, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 104.161874988944, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 53, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 53, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.55205376911452, 'increment': 54, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 110.713928758058, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.55205376911452, 'increment': 54, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 110.713928758058, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 54, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 54, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.72053561390797, 'increment': 55, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 117.434464371966, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.72053561390797, 'increment': 55, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 117.434464371966, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 55, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 55, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.88613870165196, 'increment': 56, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 124.320603073618, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 6.88613870165196, 'increment': 56, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 124.320603073618, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 56, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 56, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.05051413290862, 'increment': 57, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 131.371117206527, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.05051413290862, 'increment': 57, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 131.371117206527, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 57, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 57, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.21022056811365, 'increment': 58, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 138.58133777464, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.21022056811365, 'increment': 58, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 138.58133777464, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 58, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 58, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.36720688303465, 'increment': 59, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 145.948544657675, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.36720688303465, 'increment': 59, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 145.948544657675, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 59, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 59, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.52084991774662, 'increment': 60, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 153.469394575422, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.52084991774662, 'increment': 60, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 153.469394575422, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 60, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 60, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.67106432284142, 'increment': 61, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 161.140458898263, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.67106432284142, 'increment': 61, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 161.140458898263, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 61, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 61, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.81945766191926, 'increment': 62, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 168.959916560182, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.81945766191926, 'increment': 62, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 168.959916560182, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 62, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 62, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.96403467804308, 'increment': 63, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 176.923951238225, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 7.96403467804308, 'increment': 63, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 176.923951238225, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 63, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 63, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.10556333174781, 'increment': 64, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 185.029514569973, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.10556333174781, 'increment': 64, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 185.029514569973, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 64, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 64, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.24353786293044, 'increment': 65, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 193.273052432904, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.24353786293044, 'increment': 65, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 193.273052432904, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 65, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 65, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.37864113620333, 'increment': 66, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 201.651693569107, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.37864113620333, 'increment': 66, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 201.651693569107, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 66, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 66, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.51078755010732, 'increment': 67, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 210.162481119214, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.51078755010732, 'increment': 67, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 210.162481119214, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 67, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 67, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.64051570390474, 'increment': 68, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 218.802996823119, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.64051570390474, 'increment': 68, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 218.802996823119, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 68, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 68, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.76850118352529, 'increment': 69, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 227.571498006644, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.76850118352529, 'increment': 69, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 227.571498006644, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 69, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 69, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.89374989019844, 'increment': 70, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 236.465247896843, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 8.89374989019844, 'increment': 70, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 236.465247896843, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 70, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 70, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.01662139250328, 'increment': 71, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 245.481869289346, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.01662139250328, 'increment': 71, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 245.481869289346, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 71, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 71, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.13670416806264, 'increment': 72, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 254.618573457409, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.13670416806264, 'increment': 72, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 254.618573457409, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 72, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 72, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.25479724974754, 'increment': 73, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 263.873370707156, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.25479724974754, 'increment': 73, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 263.873370707156, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 73, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 73, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.37073508611639, 'increment': 74, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 273.244105793273, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.37073508611639, 'increment': 74, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 273.244105793273, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 74, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 74, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.48475822629464, 'increment': 75, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 282.728864019567, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.48475822629464, 'increment': 75, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 282.728864019567, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 75, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 75, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.59635833630395, 'increment': 76, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 292.325222355871, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.59635833630395, 'increment': 76, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 292.325222355871, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 76, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 76, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.70624770030058, 'increment': 77, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 302.031470056172, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.70624770030058, 'increment': 77, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 302.031470056172, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 77, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 77, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.81394275842722, 'increment': 78, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 311.845412814599, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.81394275842722, 'increment': 78, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 311.845412814599, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 78, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 78, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.91932568967703, 'increment': 79, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 321.764738504276, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 9.91932568967703, 'increment': 79, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 321.764738504276, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 79, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 79, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.0227856858904, 'increment': 80, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 331.787524190166, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.0227856858904, 'increment': 80, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 331.787524190166, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 80, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 80, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.1240918376564, 'increment': 81, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 341.911616027823, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.1240918376564, 'increment': 81, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 341.911616027823, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 81, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 81, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.2243490450178, 'increment': 82, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 352.135965072841, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.2243490450178, 'increment': 82, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 352.135965072841, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 82, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 82, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.3206938101824, 'increment': 83, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 362.456658883023, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.3206938101824, 'increment': 83, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 362.456658883023, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 83, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 83, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.4169521438206, 'increment': 84, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 372.873611026844, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.4169521438206, 'increment': 84, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 372.873611026844, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 84, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 84, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.5109588927755, 'increment': 85, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 383.384569919619, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.5109588927755, 'increment': 85, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 383.384569919619, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 85, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 85, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.6035831982631, 'increment': 86, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 393.988153117882, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.6035831982631, 'increment': 86, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 393.988153117882, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 86, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 86, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.693899004742, 'increment': 87, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 404.682052122624, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.693899004742, 'increment': 87, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 404.682052122624, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 87, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 87, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.7834622763314, 'increment': 88, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 415.465514398956, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.7834622763314, 'increment': 88, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 415.465514398956, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 88, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 88, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.8710024506654, 'increment': 89, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 426.336516849621, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.8710024506654, 'increment': 89, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 426.336516849621, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 89, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 89, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.9597377842003, 'increment': 90, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 437.296254633821, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 10.9597377842003, 'increment': 90, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 437.296254633821, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 90, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 90, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.0442116698608, 'increment': 91, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 448.340466303682, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.0442116698608, 'increment': 91, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 6, 
    # 'phase': STANDARD_PHASE, 'lpf': 448.340466303682, 'equilibrium': 6})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 91, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 91, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.1270574162477, 'increment': 92, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 7, 
    # 'phase': STANDARD_PHASE, 'lpf': 459.46752371993, 'equilibrium': 7})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.1270574162477, 'increment': 92, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 7, 
    # 'phase': STANDARD_PHASE, 'lpf': 459.46752371993, 'equilibrium': 7})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 92, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 92, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.2067475855585, 'increment': 93, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 8, 
    # 'phase': STANDARD_PHASE, 'lpf': 470.674271305488, 'equilibrium': 8})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.2067475855585, 'increment': 93, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 8, 
    # 'phase': STANDARD_PHASE, 'lpf': 470.674271305488, 'equilibrium': 8})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 93, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 93, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.2883743761859, 'increment': 94, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 8, 
    # 'phase': STANDARD_PHASE, 'lpf': 481.962645681674, 'equilibrium': 8})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.2883743761859, 'increment': 94, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 8, 
    # 'phase': STANDARD_PHASE, 'lpf': 481.962645681674, 'equilibrium': 8})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 94, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 94, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.3655450282219, 'increment': 95, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 8, 
    # 'phase': STANDARD_PHASE, 'lpf': 493.328190709896, 'equilibrium': 8})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.3655450282219, 'increment': 95, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 8, 
    # 'phase': STANDARD_PHASE, 'lpf': 493.328190709896, 'equilibrium': 8})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 95, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 95, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.4429723466631, 'increment': 96, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 8, 
    # 'phase': STANDARD_PHASE, 'lpf': 504.771163056559, 'equilibrium': 8})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.4429723466631, 'increment': 96, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 8, 
    # 'phase': STANDARD_PHASE, 'lpf': 504.771163056559, 'equilibrium': 8})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 96, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 96, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.5090903803119, 'increment': 97, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 9, 
    # 'phase': STANDARD_PHASE, 'lpf': 516.280253436871, 'equilibrium': 9})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.5090903803119, 'increment': 97, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 9, 
    # 'phase': STANDARD_PHASE, 'lpf': 516.280253436871, 'equilibrium': 9})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 97, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 97, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.5874836424827, 'increment': 98, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 9, 
    # 'phase': STANDARD_PHASE, 'lpf': 527.867737079354, 'equilibrium': 9})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.5874836424827, 'increment': 98, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 9, 
    # 'phase': STANDARD_PHASE, 'lpf': 527.867737079354, 'equilibrium': 9})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 98, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 98, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(WARNING, {'phase': STANDARD_PHASE, 
    # 'message': 'MOMENT equilibrium accepted using the alternate tolerance.', 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(WARNING, {'phase': STANDARD_PHASE, 
    # 'message': 'MOMENT equilibrium accepted using the alternate tolerance.', 
    # 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.6644919062404, 'increment': 99, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 10, 
    # 'phase': STANDARD_PHASE, 'lpf': 539.532228985594, 'equilibrium': 10})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 1, 
    # 'incrementLPF': 11.6644919062404, 'increment': 99, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 10, 
    # 'phase': STANDARD_PHASE, 'lpf': 539.532228985594, 'equilibrium': 10})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 99, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 99, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': ' 1U', 
    # 'incrementLPF': 11.8077614726334, 'increment': 100, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 539.532228985594, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': ' 1U', 
    # 'incrementLPF': 11.8077614726334, 'increment': 100, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 5, 
    # 'phase': STANDARD_PHASE, 'lpf': 539.532228985594, 'equilibrium': 5})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 2, 
    # 'incrementLPF': 2.93014138045372, 'increment': 100, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 542.462370366048, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(STATUS, {'attempts': 2, 
    # 'incrementLPF': 2.93014138045372, 'increment': 100, 'step': 1, 
    # 'jobName': 'poly-riks', 'severe': 0, 'iterations': 4, 
    # 'phase': STANDARD_PHASE, 'lpf': 542.462370366048, 'equilibrium': 4})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 100, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(ODB_FRAME, {'phase': STANDARD_PHASE, 'step': 0, 
    # 'frame': 100, 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(COMPLETED, {'phase': STANDARD_PHASE, 
    # 'message': 'Analysis phase complete', 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(COMPLETED, {'phase': STANDARD_PHASE, 
    # 'message': 'Analysis phase complete', 'jobName': 'poly-riks'})
# mdb.jobs['poly-riks']._Message(JOB_COMPLETED, {
    # 'time': 'Sat Nov 05 00:01:44 2016', 'jobName': 'poly-riks'})
