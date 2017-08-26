# Abaqus python method creating a model of a closed polygon with initial imperfections
# Current version #
# Version 8
#
# Version history #
#
# Version 7
# 15/08/2017
# Corrected sigma_cr_shell. Some EN1993-1-6 formulas were written wrong
#
# Version 6
# 13/07/2017
# Added a second function aiding the design of the test specimens
# The new function accepts n_sides, thickness, class, fy and designs a specimen
#
# Version 5
# 05/06/2017
# Modified to run a GNIA
#
# Version4
# The imperfection amplitude is calculated as dimple imperfections
# acc. to EC3-1-6
#
# Version 3
# 31/05/2017
# Accepts number of meridional half-waves instead of the inperfection angle.
#
# Version 2
# 30/05/2017
# This version accepts theta_yoshi instead of "number_of_twists", which is the angle
# of the Yoshimura pattern to the meridians, 0 <theta_yosi > pi/2

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
import EN_tools as en


def modeler(
        n_sides = None, 
        r_circle = None, 
        p_classification = None, 
        column_length = None,
        n_of_waves = None, 
        m_of_waves = None, 
        fab_class = None, 
        f_yield = None, 
        nominal_fy = None,
        n_eigen = None, 
        IDstring = None, 
        proj_name = None, 
        ):
    
    # Number of polygon sides
    if n_sides is None:
        n_sides = 6
    else:
        n_sides = int(n_sides)
    
    # Radius of the circle of perimeter equal to the polygon perimeter
    if r_circle is None:
        r_circle = 250.
    else:
        r_circle = float(r_circle)
    
    # C / (t * epsilon) ratio
    if p_classification is None:
        p_classification = 42
    else:
        p_classification = float(p_classification)
    
    # Column length
    if column_length is None:
        column_length = 2 * pi * r_circle
    else:
        column_length = float(column_length)
    
    # If no waves are given, a buckling model is ran and used for GMNIA imperfections
    if n_of_waves is None and m_of_waves is None:
        bckl_flag = True
    else:
        bckl_flag = False
    
    # Number of circumverential waves for imperfection
    if n_of_waves is None:
        n_of_waves = int(n_sides / 2)
    else:
        n_of_waves = int(n_of_waves)
    
    # Number of meridional waves for imperfection
    if m_of_waves is None:
        m_of_waves = int(n_of_waves)
    else:
        m_of_waves = int(m_of_waves)
    
    # Amplitude of imperfection
    if fab_class is None:
        fab_class = 'fcA'
    elif fab_class is not ('fcA' or 'fcB' or 'fcB'):
        print('Invalid fabrication class input. Choose between \'fcA\', \'fcB\' and \'fcC\' ')
    
    # Yield stress
    if f_yield is None:
        f_yield = 381.
    else:
        f_yield = float(f_yield)
    
    # Nominal yield strength (string)
    if nominal_fy is None:
        nominal_fy = 'S355'
    else:
        nominal_fy = str(nominal_fy)
    
    # Number of eigenvalues requested
    if n_eigen is None:
        n_eigen = 1
    else:
        n_eigen = int(n_eigen)
    
    # ID of the current job
    if IDstring is None:
        IDstring = 'NA'
    
    # Name of the parametric project
    if proj_name is None:
        proj_name = 'NA'
    
    # Create a new model database. This will also start a new journal file for the current session.
    Mdb()
    ## MODEL ##
    
    # Calculate the cs geometric and resistance propertiesusinc cs_calculator
    # and break the return tuple to separate variables
    properties = cs_calculator(
        n_sides = n_sides, 
        r_circle = r_circle, 
        p_classification = p_classification, 
        column_length = column_length,
        f_yield = f_yield,
        fab_class = fab_class,
        )
    
    r_circum = properties[0]
    shell_thickness = properties[1]
    N_pl_Rd = properties[2]
    N_b_sh = properties[3]
    x_corners = properties[4]
    y_corners = properties[5]
    
    # Epsilon for the material
    epsilon = sqrt(235. / f_yield)
    
    # Diameter of the circumferential circle
    diameter = 2 * r_circum
    
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
        sheetSize = 2 * r_circum
        )
    
    # Draw lines sides on the sketch for the polygon
    for current_corner in range(n_sides):
        cs_sketch.Line(
            point1 = (x_corners[current_corner], y_corners[current_corner]),
            point2 = (x_corners[current_corner - 1], y_corners[current_corner - 1])
            )
    
    # Fillet corners
    for current_side in range(n_sides):
        cs_sketch.FilletByRadius(
            curve1 = cs_sketch.geometry.items()[current_side-1][1],
            curve2 = cs_sketch.geometry.items()[current_side][1],
            nearPoint1 = (0, 0),
            nearPoint2 = (0, 0),
            radius = 3 * shell_thickness
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
    
    # Seed the part according to the shell thickness
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
        maxNumInc=3,
        extrapolation=LINEAR,
        initialArcInc=0.1,
        minArcInc=1e-07,
        totalArcLength=2
        )
    
    # Rename the material
    riks_mdl.materials.changeKey(
        fromName='elastic',
        toName=nominal_fy)
    
    # Change to plastic material, optim355
    riks_mdl.materials[nominal_fy].Plastic(
        table= xtr.plastic_table(nominal = nominal_fy)
        )
    
    # Change the section material name accordingly
    riks_mdl.sections['shell_section'].setValues(
        material=nominal_fy,
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
        cf3 = -N_pl_Rd,
        createStepName = step_name, 
        distributionType = UNIFORM,
        field = '',
        localCsys = None,
        name = 'compression', 
        region = rp_head_set
        )
    
    ###### END RIKS MODEL ######
    
    ###### BCKL JOB ######
    
    if bckl_flag is True:
        # Edit the keywords for the buckling model to write 'U' on file
        bckl_model.keywordBlock.synchVersions(storeNodesAndElements = False)
        bckl_model.keywordBlock.insert(xtr.GetBlockPosition(bckl_model,'*End Step')-1, '*NODE FILE\nU')
        
        # Create the job
        bckl_job = mdb.Job(
            atTime = None,
            contactPrint = OFF,
            description = '',
            echoPrint = OFF, 
            explicitPrecision = SINGLE,
            getMemoryFromAnalysis = True,
            historyPrint = OFF, 
            memory = 90,
            memoryUnits = PERCENTAGE,
            model = 'bckl_model',
            modelPrint = OFF, 
            multiprocessingMode = DEFAULT,
            name = 'BCKL-'+IDstring,
            nodalOutputPrecision = SINGLE, 
            numCpus = 1,
            numGPUs = 0,
            queue = None,
            resultsFormat = ODB,
            scratch = '',
            type = ANALYSIS,
            userSubroutine = '',
            waitHours = 0,
            waitMinutes = 0
            )
        
        # Submit buckling job
        bckl_job.submit(consistencyChecking=OFF)
        bckl_job.waitForCompletion()
    
    ###### END BCKL JOB ######
    
    ###### BCKL post processing ######
    # Open the buckling step odb file
    bckl_odb = odbAccess.openOdb(path='BCKL-'+IDstring+'.odb')
    bckl_step = bckl_odb.steps['bckl']
    
    # Gather the eigenvalues
    eigenvalues = ()
    eigen_string = ""
    for J_eigenvalues in range(1, n_eigen + 1):
        current_eigen = float(bckl_step.frames[J_eigenvalues].description[-11:])
        eigenvalues = eigenvalues + (current_eigen, )
        eigen_string = eigen_string + "%.3E "%(current_eigen)
    
    
    ###### IMPERFECTIONS #####
    
    # Edit the coordinates of nodes to get the imperfect shape.
    # In case specific number of waves are given as input, the imperfections
    # are applied directly to the mesh by the wave functions
    # In case no imperfections are given, imperfections are introduced 
    # using the first eigenmode of a buckling analysis. The amplitude of
    # this imperfection is regulated according to the 
    
    if bckl_flag is True:
        # Calculate the proportionality of the two imperfection types
        # based on the deviation of the 1st eigenvalue to the plate critical load
        #
        # Plate critical stress
        sigma_cr = en.sigma_cr_plate(shell_thickness, (pi * diameter / n_sides))
        N_cr_pl = pi * diameter * shell_thickness * sigma_cr
        N_cr_sh = en.N_cr_shell(shell_thickness, r_circle, column_length)
        prop_a = abs(eigenvalues[0] - N_cr_pl)
        prop_b = abs(eigenvalues[0] - N_cr_sh)
        
        # find the maximum displacement from the buckling analysis
        bckl_odb = xtr.open_odb('BCKL-'+IDstring+'.odb')
        Umax = xtr.field_max(bckl_odb, ['U', 'Magnitude'])
        odbAccess.closeOdb(bckl_odb)
        
        # Plate-like imperfection half wavelength
        # perimeter divided by the number of waves
        l_g_I = pi * diameter / (n_sides)
        
        # Plate-like imperfection half wavelength
        l_g_II = pi * diameter / (2 * floor(n_sides / 4))
        
        # l_g is formed from l_g_I and l_g_II.
        # the contribution of each wavelength is based on the deviation of 
        # the eigenvalue analysis to the EC N_cr
        l_g = l_g_I * (prop_b / (prop_a + prop_b)) + l_g_II * (prop_a / (prop_a + prop_b))
        
        # Calculate target maximum imperfection displacement (fabrication class)
        # for a length lg calculated between 1 and 2 sides (plate and spillover waves)
        u_tot = l_g * fabclass_2_umax(fab_class)
        
        # Imperfection amplitude
        a_imp = u_tot / Umax
        
        # Edit the keywords for the compression riks model to include imperfections from buckling analysis and to output the RIKS_NODE for GN_killer to work
        riks_mdl.keywordBlock.synchVersions(storeNodesAndElements=False)
        riks_mdl.keywordBlock.replace(xtr.GetBlockPosition(riks_mdl, '*step')-1, 
        '\n** ----------------------------------------------------------------\n** \n**********GEOMETRICAL IMPERFECTIONS\n*IMPERFECTION,FILE=BCKL-'+IDstring+',STEP=1\n1,'+str(a_imp)+'\n**')
        riks_mdl.keywordBlock.insert(xtr.GetBlockPosition(riks_mdl,'*End Step')-1, '\n*NODE FILE, GLOBAL=YES, NSET=RIKS_NODE\nU')
    
    else:
        # Circumferencial half wavelength.
        # perimeter divided by the number of waves
        l_g_circum = pi * diameter / (2 * n_of_waves)
        
        # Meridional half wavelength
        l_g_meridi = column_length / (2 * m_of_waves)
        
        # l_gx is calculated as the mean value of the two wavelengths
        # (This requires justification)
        l_g = min(l_g_circum, l_g_meridi)
        
        # Loop through the nodes of the mesh and apply displacements
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
                offset1 = fabclass_2_umax(fab_class) * l_g * (circum_wave * meridi_wave) * cos(crrnt_pt_angle),
                offset2 = fabclass_2_umax(fab_class) * l_g * (circum_wave * meridi_wave) * sin(crrnt_pt_angle)
                )
    
    ######END IMPERFECTIONS ####
    
    ###### RIKS JOB #######
        
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

    ## RIKS post processing
    ## Find max LPF
    LPF = xtr.history_max('RIKS-'+IDstring, step_name)
    #
    ## Create and populate an output text with model information
    #
    #out_file = open('./'+IDstring+'_info.dat', 'w')
    #out_file.write('\n-GEOMETRIC CHARACTERISTICS\n')
    #out_file.write('Number of corners:....................................................... '+str(n_sides)+'\n')
    #out_file.write('Diameter of equal perimeter circle:...................................... '+str(2 * r_circle)+' [mm]\n')
    #out_file.write('Diameter of the circumscribed circle:.................................... '+str(2 * r_circum)+' [mm]\n')
    #out_file.write('Width of each side:...................................................... '+str(w_side)+' [mm]\n')
    #out_file.write('Perimeter:............................................................... '+str(perimeter)+' [mm]\n')
    #out_file.write('Total length:............................................................ '+str(column_length/1000)+' [m]\n')
    #out_file.write('Profile thickness:....................................................... '+str(shell_thickness)+' [mm]\n')
    #out_file.write('Bending radius for the creases of the polygon (midline):................. '+str(3 * shell_thickness)+' [mm]\n')
    #out_file.write('Cross-sectional Area:.................................................... '+str(Area)+' [mm^2]\n')
    #out_file.write('\n-STRUCTURAL CHARACTERISTICS'+'\n')
    #out_file.write('Yield strength:.......................................................... '+str(f_yield)+' [MPa]\n')
    #out_file.write('epsilon:................................................................. '+str(epsilon)+'\n')
    #out_file.write('Cross-section classification (as plate):................................. '+str(p_classification)+'\n')
    #out_file.write('Cross-section classification (as tube):.................................. '+str(t_classification)+'\n')
    #out_file.write('Class 4 reduction factor:................................................ '+str(rho)+'\n')
    #out_file.write('Effective area:.......................................................... '+str(A_eff)+' [mm^2]\n')
    #out_file.write('Plate critical stress:................................................... '+str(sigma_cr_plate)+' [MPa]\n')
    #out_file.write('Shell critical stress:................................................... '+str(sigma_x_Rcr)+' [MPa]\n')
    #out_file.write('Plate critical load:..................................................... '+str(N_cr_plate/1000.)+' [kN]\n')
    #out_file.write('Shell critical load:..................................................... '+str(N_cr_shell/1000.)+' [kN]\n')
    #out_file.write('Plastic resistance, N_pl_rd:............................................. '+str(N_pl_Rd/1000.)+' [kN]\n')
    #out_file.write('\n-RESULTS\n')
    ##for J_eigenvalues in range(1, n_eigen + 1):
    ##    out_file.write('Eigen value '+"%02d"%J_eigenvalues+':......................................................... '+str(eigenvalues[J_eigenvalues-1]/1000)+' [kN]\n')
    #
    #out_file.close()
    
    # Compile a string with the model parameters and results
    # This string is returned
    
    
    #out_file = open('./' + proj_name +'.dat', 'a')
    return_string = ("%03d %07.3f %07.3f %07.3f %03d %05.3f %05.1f %.3E %.3E %05.3f "
        %(
            n_sides,
            2 * r_circle,
            2 * r_circum,
            shell_thickness,
            f_yield,
            epsilon,
            p_classification,
            N_pl_Rd,
            N_b_sh,
            LPF
            )
        )
    return return_string


def class_2_radius(
        n_sides = None, 
        p_classification = None, 
        thickness = None, 
        f_yield = None, 
        ):
    
    # Number of polygon sides
    if n_sides is None:
        n_sides = 6
    else:
        n_sides = int(n_sides)
    
    # C / (t * epsilon) ratio
    if p_classification is None:
        p_classification = 42
    else:
        p_classification = float(p_classification)
    
    # Thickness
    if thickness is None:
        thickness = 2 * pi * r_circle
    else:
        thickness = float(thickness)
    
    # Yield stress
    if f_yield is None:
        f_yield = 381.
    else:
        f_yield = float(f_yield)
    
    # Epsilon for the material
    epsilon = sqrt(235. / f_yield)
    
    # Radius of the equal perimeter cylinder
    r_circle = n_sides * thickness * epsilon * p_classification / (2 * pi)
    
    # Return radius of the cylinder of equal perimeter
    return r_circle


def cs_calculator(
        n_sides = None, 
        r_circle = None, 
        p_classification = None, 
        column_length = None,
        f_yield = None,
        fab_class = None
        ):
    
    # Number of polygon sides
    if n_sides is None:
        n_sides = 6
    else:
        n_sides = int(n_sides)
    
    # Radius of the equal perimeter cylinder
    if r_circle is None:
        r_circle = 250.
    else:
        r_circle = float(r_circle)
    
    # C / (t * epsilon) ratio
    if p_classification is None:
        p_classification = 42.
    else:
        p_classification = float(p_classification)
    
    # Column length
    if column_length is None:
        column_length = 2 * pi * r_circle
    else:
        column_length = float(column_length)
    
    # Yield stress
    if f_yield is None:
        f_yield = 381.
    else:
        f_yield = float(f_yield)

    # Fabrication class
    if fab_class is None:
        fab_class = 'fcA'

    ## END INPUT ##
    
    ## GEOMETRY ##
    
    # Epsilon for the material
    epsilon = sqrt(235. / f_yield)
    
    # Radius of the polygon's circumscribed circle
    r_circum = (pi * r_circle) / (n_sides * sin(pi / n_sides))
    
    # Diameter
    diameter = 2 * r_circum
    
    # Central angles
    theta = 2 * pi / n_sides
    
    # Width of each side
    w_side = diameter * sin(pi / n_sides)
    
    # Thickness the given class (classification as plated, not tube)
    thickness = (diameter * sin(theta / 2)) / (p_classification * epsilon)
    
    # Polar coordinate of ths polygon vertices on the cross-section plane
    phii = []
    for i_index in range(n_sides):
        phii.append(i_index * theta)
    
    # Polygon corners coordinates
    x_corners = r_circum * np.cos(phii)
    y_corners = r_circum * np.sin(phii)
    
    # Axial compression resistance , Npl
    N_pl_Rd = n_sides * en.N_pl_Rd(thickness, w_side, f_yield)
    
    # Compression resistance of equivalent cylindrical shell
    N_b_Rd_shell = 2 * pi * r_circle * thickness * en.sigma_x_Rd(
        thickness, 
        r_circle, 
        column_length, 
        f_yield, 
        fab_quality = fab_class, 
        gamma_M1 = 1.
		)
        
    ## END GEOMETRY ##
    
    # Return values
    return r_circum, thickness, N_pl_Rd, N_b_Rd_shell, x_corners, y_corners


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