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

# Additional modules
import numpy as np
import abq_toolset as xtr
import os
import sys
import odbAccess
from shutil import copyfile
from EN_tools import N_pl_Rd as single_plate_Rd
from EN_tools import sigma_x_Rd as shell_buckling_stress
from closed_polygons import cs_calculator


def polygon_2_cylinder(n_sides, r_circle, p_classification, f_yield):
    ## MODEL ##
    
    # Calculate the cs geometric and resistance propertiesusinc cs_calculator
    # and break the return tuple to separate variables
    properties = cs_calculator(
        n_sides = n_sides, 
        r_circle = r_circle, 
        p_classification = p_classification, 
        f_yield = f_yield
        )
    
    shell_thickness = properties[1]
    N_b_sh = properties[3]
    
    # Return values
    return shell_thickness


def modeler(
        r_circle = None,
        shell_thickness = None,
        column_length = None,
        f_yield = None,
        n_eigen = None,
        n_of_waves = None,
        m_of_waves = None,
        u_max = None,
        IDstring = None, 
        proj_name = None, 
        ):
    
    # Radius
    if r_circle is None:
        r_circle = 250.
    else:
        r_circle = float(r_circle)
    
    # Thickness
    if shell_thickness is None:
        shell_thickness = 1.
    else:
        shell_thickness = float(shell_thickness)
    
    # Length
    if column_length is None:
        column_length = 2 * pi * r_circle
    else:
        column_length = float(column_length)
    
    # Yield stress
    if f_yield is None:
        f_yield = 381.
    else:
        f_yield = float(f_yield)
    
    # Number of eigenvalues
    if n_eigen is None:
        n_eigen = 1
    else:
        n_eigen = int(n_eigen)
    
    # Number of circumverential waves for imperfection
    if n_of_waves is None:
        n_of_waves = 2
    else:
        n_of_waves = int(n_of_waves)
    
    # Number of meridional waves for imperfection
    if m_of_waves is None:
        m_of_waves = 2
    else:
        m_of_waves = int(m_of_waves)
    
    # Amplitude of imperfection
    if u_max is None:
        u_max = 0.016
    else:
        u_max = float(u_max)
    
    # ID of the current job
    if IDstring is None:
        IDstring = 'NA'
    
    # Name of the parametric project
    if proj_name is None:
        proj_name = 'NA'
    
    # Calculate the cylinder's resistance to meridional compression
    # Using the functions of EN_tools.py
    N_shell_Rd = 2 * pi * r_circle * shell_thickness * shell_buckling_stress(
        shell_thickness,
        r_circle,
        column_length,
        f_yield,
        3,
        1.
        )
    
    # Create a new model database. This will also start a new journal file for the current session.
    Mdb()
    ## MODEL ##
    
    # Epsilon for the material
    epsilon = sqrt(235. / f_yield)
    
    # Diameter of the circumferential circle
    diameter = 2 * r_circle
    
    # Buckling model name
    bckl_model_name = 'bckl_model'
    
    # Change the pre-existing model name
    mdb.models.changeKey(
        fromName = 'Model-1',
        toName = bckl_model_name
        )
    
    #Create a variable for the model
    bckl_model = mdb.models.values()[0]
    
    # Create sketch
    cs_sketch = bckl_model.ConstrainedSketch(
        name = 'cs_sketch',
        sheetSize = 2 * r_circle
        )
    
    # Draw the cs circle
    cs_sketch.CircleByCenterPerimeter(
        center=(0.0, 0.0),
        point1=(r_circle, 0.0)
        )

    # Create the part
    p_part = bckl_model.Part(
        dimensionality = THREE_D,
        name = 'short_column',
        type = DEFORMABLE_BODY
        )
    
    p_part.BaseShellExtrude(
        depth = column_length,
        sketch = cs_sketch
        )
    
    # Create material
    el_material = bckl_model.Material(name='elastic')
    
    el_material.Elastic(table=((210000.0, 0.3), ))
    
    # Create shell section
    shell_section = bckl_model.HomogeneousShellSection(
        idealization = NO_IDEALIZATION, 
        integrationRule = SIMPSON,
        material = 'elastic',
        name = 'shell_section', 
        numIntPts = 5,
        poissonDefinition = DEFAULT,
        preIntegrate = OFF,
        temperature = GRADIENT,
        thickness = shell_thickness,
        thicknessField = '',
        thicknessModulus = None, 
        thicknessType = UNIFORM,
        useDensity = OFF
        )
    
    # Create a set with all the faces
    all_faces_set = p_part.Set(
        faces = p_part.faces[:],
        name = 'all_faces'
        )
    
    # Assign the section to the set
    p_part.SectionAssignment(
        offset = 0.0, 
        offsetField = '',
        offsetType = MIDDLE_SURFACE,
        region = all_faces_set,
        sectionName = 'shell_section',
        thicknessAssignment = FROM_SECTION
        )
    
    # Seed the part according to the shell diameter
    p_part.seedPart(
        deviationFactor = 0.1, 
        minSizeFactor = 0.1,
        size = diameter / 65
        )
    
    # Mesh the part
    p_part.generateMesh()
    
    # Create variable and coordinate system for the assembly
    r_assembly = bckl_model.rootAssembly
    r_assembly.DatumCsysByDefault(CARTESIAN)
    
    # Create instance
    column_instance = r_assembly.Instance(
        dependent=ON,
        name='short_column', 
        part=p_part
        )
    
    # Create reference points at the ends of the column for BC couplings
    r_assembly.ReferencePoint(point=(0.0, 0.0, 0.0))
    r_assembly.ReferencePoint(point=(0.0, 0.0, column_length))
    
    rp_base = r_assembly.referencePoints.items()[1][1]
    rp_head = r_assembly.referencePoints.items()[0][1]
    
    # Create sets for the two reference points
    rp_base_set = r_assembly.Set(
        name='base_rp',
        referencePoints = (rp_base, ))
    
    rp_head_set = r_assembly.Set(
        name='head_rp',
        referencePoints = (rp_head, ))
    
    # Create sets for the base and the head of the column
    base_edges_set = r_assembly.Set(
        edges = column_instance.edges.getByBoundingBox(-diameter, -diameter, 0, diameter, diameter, 0),
        name = 'base_edges'
        )
    
    head_edges_set = r_assembly.Set( 
        edges = column_instance.edges.getByBoundingBox(-diameter, -diameter, column_length, diameter, diameter, column_length),
        name = 'head_edges'
        )
    
    # Create column end couplings
    # Current coupling settings restrain the shell membrain rotation
    # For free edge shell rotation, change to: ur1 = OFF, ur2 = OFF
    base_coupling = bckl_model.Coupling(
        controlPoint = rp_base_set,
        couplingType = KINEMATIC,
        influenceRadius = WHOLE_SURFACE,
        localCsys = None,
        name = 'base_coupling',
        surface = base_edges_set,
        u1 = ON, u2 = ON, u3 = ON, ur1 = ON, ur2 = ON, ur3 = ON
        )
    
    head_coupling = bckl_model.Coupling(
        controlPoint = rp_head_set,
        couplingType = KINEMATIC,
        influenceRadius = WHOLE_SURFACE,
        localCsys = None,
        name = 'head_coupling',
        surface = head_edges_set,
        u1 = ON, u2 = ON, u3 = ON, ur1 = ON, ur2 = ON, ur3 = ON
        )
    
    # Create buckling analysis step
    bckl_model.BuckleStep(
        name = 'bckl',
        numEigen = n_eigen,
        previous = 'Initial', 
        vectors = 8,
        maxIterations = 1000
        )
    
    # Apply concentrated load
    bckl_model.ConcentratedForce(
        cf3 = -1,
        createStepName = 'bckl', 
        distributionType = UNIFORM,
        field = '',
        localCsys = None,
        name = 'compression', 
        region = rp_head_set
        )
    
    # Fix column base
    bckl_model.DisplacementBC(
        amplitude = UNSET,
        createStepName = 'Initial', 
        distributionType = UNIFORM,
        fieldName = '',
        localCsys = None,
        name = 'fix_base', 
        region = rp_base_set,
        u1 = SET, u2 = SET, u3 = SET, ur1 = SET, ur2 = SET, ur3 = SET
        )
    
    # Fix column head
    bckl_model.DisplacementBC(
        amplitude = UNSET,
        createStepName = 'Initial', 
        distributionType = UNIFORM,
        fieldName = '',
        localCsys = None,
        name = 'fix_head', 
        region = rp_head_set,
        u1 = SET, u2 = SET, u3 = UNSET, ur1 = SET, ur2 = SET, ur3 = SET
        )
    
    # Set field output requests
    bckl_model.fieldOutputRequests['F-Output-1'].setValues(
        variables = ('U',)
        )
    
    # Apply imperfections
    # Edit the coordinates of nodes to get the imperfect shape
    
    # Circumferencial half wavelength.
    # perimeter divided by the number of waves
    l_gx_circum = pi * diameter / (2 * n_of_waves)
    
    # Meridional half wavelength
    l_gx_meridi = column_length / (2 * m_of_waves)
    
    # l_gx is calculated as the mean value of the two wavelengths
    # (This requires justification)
    l_gx = min(l_gx_circum, l_gx_meridi)
    
    #theta_yoshi = (m_of_waves * pi * 2 * r_circle) / (n_of_waves * column_length)
    
    #n_twists = column_length * tan(theta_yoshi) / (2 * pi * r_circle)
    
    for j in range(len(column_instance.nodes)):
        xi = column_instance.nodes[j].coordinates[0]
        yi = column_instance.nodes[j].coordinates[1]
        zi = column_instance.nodes[j].coordinates[2]
        if xi > 1e-14:
            crrnt_pt_angle = atan(yi / xi)
        elif xi < -1e-14:
            crrnt_pt_angle = pi + atan(yi / xi)
        else:
            crrnt_pt_angle = pi
        
        circum_wave = sin(n_of_waves * crrnt_pt_angle)
        meridi_wave = sin(m_of_waves * 2 * pi * zi / column_length)
        
        r_assembly.editNode(
            nodes = column_instance.nodes[j],
            offset1 = u_max * l_gx * (circum_wave * meridi_wave) * cos(crrnt_pt_angle),
            offset2 = u_max * l_gx * (circum_wave * meridi_wave) * sin(crrnt_pt_angle)
            )
    
    ## Save the model
    mdb.saveAs(pathName=os.getcwd()+'/'+IDstring+'.cae')
    
    #### RIKS MODEL ####
    
    # Riks model name
    riks_model_name = 'riks_model'
    # Copy the buckling model
    riks_mdl = mdb.Model(
        name = riks_model_name,
        objectToCopy = bckl_model
        )
    
    # Delete buckling step
    del riks_mdl.steps['bckl']
    
    # Create RIKS step
    step_name = 'RIKS'
    riks_mdl.StaticRiksStep(
        name=step_name,
        previous='Initial',
        nlgeom=ON,
        maxNumInc=1,
        extrapolation=LINEAR,
        initialArcInc=0.1,
        minArcInc=1e-07,
        totalArcLength=2
        )
    
    # Rename the material
    riks_mdl.materials.changeKey(
        fromName='elastic',
        toName='optim355')
    
    # Change to plastic material, optim355
    riks_mdl.materials['optim355'].Plastic(
        table=(
            (381.1, 0.0),
            (391.2, 0.0053),
            (404.8, 0.0197),
            (418.0, 0.0228),
            (444.2, 0.0310),
            (499.8, 0.0503),
            (539.1, 0.0764),
            (562.1, 0.1009),
            (584.6, 0.1221),
            (594.4, 0.1394)
            )
        )
    
    # Change the section material name accordingly
    riks_mdl.sections['shell_section'].setValues(
        material='optim355',
        )
    
    # Create a set to act as a control point for the solver killer subroutine
    assmbl_riks = riks_mdl.rootAssembly
    
    assmbl_riks.Set(
        name='RIKS_NODE', 
        nodes=(assmbl_riks.instances['short_column'].nodes[1:2], )
        )
    
    # Set history output request for displacement
    disp_history = riks_mdl.HistoryOutputRequest(
        createStepName = step_name,
        name = 'disp', 
        rebar = EXCLUDE,
        region = rp_head_set, 
        sectionPoints = DEFAULT,
        variables = ('U3', )
        )
    
    # Set history output request for load
    load_history = riks_mdl.HistoryOutputRequest(
        createStepName = step_name,
        name = 'load', 
        rebar = EXCLUDE,
        region = rp_base_set, 
        sectionPoints = DEFAULT,
        variables = ('RF3', )
        )
    
    # Delete pre-existing history request: H-Output-1
    riks_mdl.historyOutputRequests.delete(['H-Output-1'])
    
    # Apply concentrated load
    riks_mdl.ConcentratedForce(
        cf3 = -N_shell_Rd,
        createStepName = step_name, 
        distributionType = UNIFORM,
        field = '',
        localCsys = None,
        name = 'compression', 
        region = rp_head_set
        )
    
    ###### END RIKS MODEL ######
    
    ###### BCKL JOB ######
    
    ## Edit the keywords for the buckling model to write 'U' on file
    #bckl_model.keywordBlock.synchVersions(storeNodesAndElements = False)
    #bckl_model.keywordBlock.insert(xtr.GetBlockPosition(bckl_model,'*End Step')-1, '*NODE FILE\nU')
    #
    ## Create the job
    #bckl_job = mdb.Job(
    #    atTime = None,
    #    contactPrint = OFF,
    #    description = '',
    #    echoPrint = OFF, 
    #    explicitPrecision = SINGLE,
    #    getMemoryFromAnalysis = True,
    #    historyPrint = OFF, 
    #    memory = 90,
    #    memoryUnits = PERCENTAGE,
    #    model = 'bckl_model',
    #    modelPrint = OFF, 
    #    multiprocessingMode = DEFAULT,
    #    name = 'BCKL-'+IDstring,
    #    nodalOutputPrecision = SINGLE, 
    #    numCpus = 1,
    #    numGPUs = 0,
    #    queue = None,
    #    resultsFormat = ODB,
    #    scratch = '',
    #    type = ANALYSIS,
    #    userSubroutine = '',
    #    waitHours = 0,
    #    waitMinutes = 0
    #    )
    #
    ## Submit buckling job
    #bckl_job.submit(consistencyChecking=OFF)
    #bckl_job.waitForCompletion()
    
    ###### END BCKL JOB ######
    
    ###### RIKS JOB #######
    
    ## find the maximum displacement from the buckling analysis
    #bckl_odb = xtr.open_odb('BCKL-'+IDstring+'.odb')
    #Umax = xtr.max_result(bckl_odb, ['U', 'Magnitude'])
    #odbAccess.closeOdb(bckl_odb)
    #
    ## Imperfection amplitude
    #a_imp = w_side / (200 * Umax)
    #
    ## Edit the keywords for the compression riks model to include imperfections from buckling analysis and to output the RIKS_NODE for GN_killer to work
    riks_mdl.keywordBlock.synchVersions(storeNodesAndElements=False)
    #riks_mdl.keywordBlock.replace(xtr.GetBlockPosition(riks_mdl, '*step')-1, 
    #'\n** ----------------------------------------------------------------\n** \n**********GEOMETRICAL IMPERFECTIONS\n*IMPERFECTION,FILE=BCKL'+IDstring+',STEP=1\n1,'+str(a_imp)+'\n**')
    riks_mdl.keywordBlock.insert(xtr.GetBlockPosition(riks_mdl,'*End Step')-1, '\n*NODE FILE, GLOBAL=YES, NSET=RIKS_NODE\nU')
    
    # Create the RIKS job
    riks_job = mdb.Job(
        atTime = None,
        contactPrint = OFF,
        description = '',
        echoPrint = OFF, 
        explicitPrecision = SINGLE,
        getMemoryFromAnalysis = True,
        historyPrint = OFF, 
        memory = 90,
        memoryUnits = PERCENTAGE,
        model = riks_model_name,
        modelPrint = OFF, 
        multiprocessingMode = DEFAULT,
        name = 'RIKS-'+IDstring,
        nodalOutputPrecision = SINGLE, 
        numCpus = 1,
        numGPUs = 0,
        queue = None,
        resultsFormat = ODB,
        scratch = '',
        type = ANALYSIS,
        #userSubroutine='../GN_Riks_killer.f',
        waitHours = 0,
        waitMinutes = 0
        )
    
    # Submit RIKS job
    riks_job.submit(consistencyChecking = OFF)
    riks_job.waitForCompletion()
    
    # Save the model
    mdb.saveAs(pathName=os.getcwd()+'/'+IDstring+'.cae')
    
    ##### END RIKS JOB ##########
    
    ##### POST-PROCESSING #####\
    
    ## BCKL post processing
    ## Open the buckling step odb file
    #bckl_odb = odbAccess.openOdb(path='BCKL-'+IDstring+'.odb')
    #bckl_step = bckl_odb.steps['bckl']
    #
    ## Gather the eigenvalues
    #eigenvalues = ()
    #eigen_string = ""
    #for J_eigenvalues in range(1, n_eigen + 1):
    #    current_eigen = float(bckl_step.frames[J_eigenvalues].description[-11:])
    #    eigenvalues = eigenvalues + (current_eigen, )
    #    eigen_string = eigen_string + "%.3E "%(current_eigen)
    #
    ## RIKS post processing
    ## Find max LPF
    LPF = xtr.history_max('RIKS-'+IDstring, step_name)
    
    # Create and populate an output text with model information
    
    out_file = open('./'+IDstring+'_info.dat', 'w')
    out_file.write('\n-GEOMETRIC CHARACTERISTICS\n')
    out_file.write('Diameter of equal perimeter circle:...................................... '+str(2 * r_circle)+' [mm]\n')
    out_file.write('Total length:............................................................ '+str(column_length/1000)+' [m]\n')
    out_file.write('Profile thickness:....................................................... '+str(shell_thickness)+' [mm]\n')
    out_file.write('\n-STRUCTURAL CHARACTERISTICS'+'\n')
    out_file.write('Yield strength:.......................................................... '+str(f_yield)+' [MPa]\n')
    out_file.write('Resistance, N_shell_Rd:............................................. '+str(N_shell_Rd/1000.)+' [kN]\n')
    out_file.write('\n-RESULTS\n')
    out_file.write('Load proportionality factor:............................................. '+str(N_shell_Rd/1000.)+' [kN]\n')
    out_file.close()
    
    # Compile a string with the model parameters and results
    return_string = ("%03d %03d %07.3f %07.3f %03d %.3E %05.3f "
        %(
            m_of_waves,
            n_of_waves,
            diameter,
            shell_thickness,
            f_yield,
            N_shell_Rd,
            LPF
            )
        )
    
    # Return the results
    return return_string


def fabclass_2_umax(fab_class):
    # Assign imperfection amplitude, u_max acc. to the fabrication class
    if fab_class is 'fcA':
        u_max = 0.006
    elif fab_class is 'fcB':
        u_max = 0.010
    else:
        u_max = 0.016
    
    # Return values
    return u_max