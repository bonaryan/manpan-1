# TODO: Wrong comment or wrong command. Decide and adjust. (Line 580 for lg_x)
# Abaqus modules
from material import *
from part import *
from section import *
from interaction import *
from step import *
from job import *
from mesh import *
import odbAccess

# Additional modules
import numpy as np
import abq_tools as at
import steel_design as sd
import os

def class_2_radius(n_sides,
                   thickness,
                   p_classification,
                   f_yield,
                   ):
    """ Calculate radius of a polygon for given thickness and classification."""
    # Epsilon for the material
    epsilon = sqrt(235. / f_yield)
    
    # Calculate and return the radius of the equal perimeter cylinder
    return (n_sides * thickness * epsilon * p_classification / (2 * pi))


def class_2_thickness(n_sides,
                      r_circle,
                      p_classification,
                      f_yield,
                      ):
    """ Calculate the thickness of a polygon for given radius and classification."""
    # Epsilon for the material
    epsilon = np.sqrt(235. / f_yield)
    
    # Calculate and return the thickness
    return (2 * np.pi * r_circle / (n_sides * epsilon * p_classification))


# def apply_local_imperfections(model,
#                               r_circle,
#                               column_length,
#                               circum_waves,
#                               meridi_waves,
#                               fab_class,
#                               percentage,
#                               windowing=False):
#     # Create the items for the riks assembly and instance.
#     assembly = model.rootAssembly
#     instance = assembly.instances.values()[0]
#
#     # Circumferencial half wavelength.
#     # perimeter divided by the number of waves
#     l_g_circum = pi * 2  * r_circle / (2 * circum_waves)
#
#     # Meridional half wavelength
#     l_g_meridi = column_length / (2 * meridi_waves)
#     # l_gx is calculated as the mean value of the two wavelengths
#     # (This requires justification)
#     l_g = min(l_g_circum, l_g_meridi)
#
#     # Loop through the nodes of the mesh and apply displacements
#     for node in instance.nodes:
#         xi = node.coordinates[0]
#         yi = node.coordinates[1]
#         zi = node.coordinates[2]
#         if xi > 1e-14:
#             crrnt_pt_angle = atan(yi / xi)
#         elif xi < -1e-14:
#             crrnt_pt_angle = pi + atan(yi / xi)
#         else:
#             crrnt_pt_angle = pi
#
#         circum_wave = sin(circum_waves * crrnt_pt_angle)
#         meridi_wave = cos(meridi_waves * 2 * pi * (zi - column_length/2.) / column_length)
#
#         if windowing:
#             circum_window = blackman_percentage((crrnt_pt_angle + pi / 2.)/(2*pi))
#             meridi_window = blackman_percentage(zi/(column_length+0.1))
#         else:
#             circum_window = 1
#             meridi_window = 1
#
#         assembly.editNode(
#             nodes=node,
#             offset1=percentage * sd.fabclass_2_umax(fab_class) * l_g * (circum_wave * meridi_wave) * cos(crrnt_pt_angle) * (circum_window * meridi_window),
#             offset2=percentage * sd.fabclass_2_umax(fab_class) * l_g * (circum_wave * meridi_wave) * sin(crrnt_pt_angle) * (circum_window * meridi_window)
#         )
#
#
# def apply_flexural_imperfections(model,
#                                  column_length,
#                                  flex_imp):
#     # Create the items for the riks assembly and instance.
#     assembly = model.rootAssembly
#     instance = assembly.instances.values()[0]
#     flex_imperfection_amp = column_length / flex_imp
#
#     for node in instance.nodes:
#         xi = node.coordinates[0]
#         yi = node.coordinates[1]
#         zi = node.coordinates[2]
#         glob_bow = -flex_imperfection_amp * sin(pi * zi / column_length) * sin(pi / 4)
#         assembly.editNode(
#             nodes=node,
#             offset1=glob_bow,
#             offset2=glob_bow,
#         )


def calc_flex_imp(z,
                  column_length,
                  flex_imp,
                  axis_angle):
    if flex_imp==0:
        return 0, 0
    
    # Create the items for the riks assembly and instance.
    flex_imperfection_amp = column_length / flex_imp
    glob_bow = -flex_imperfection_amp * sin(pi * z / column_length) * sin(pi / 4)
    return sin(axis_angle) * glob_bow, cos(axis_angle) * glob_bow

def calc_local_imp(coords,
                   r_circle,
                   column_length,
                   circum_waves,
                   meridi_waves,
                   fab_class,
                   percentage,
                   windowing=False):
    
    # Circumferencial half wavelength.
    # perimeter divided by the number of waves
    l_g_circum = pi * 2  * r_circle / (2 * circum_waves)
    
    # Meridional half wavelength
    l_g_meridi = column_length / (2 * meridi_waves)
    # l_gx is calculated as the mean value of the two wavelengths
    # (This requires justification)
    l_g = min(l_g_circum, l_g_meridi)
    
    # Loop through the nodes of the mesh and apply displacements
    xi = coords[0]
    yi = coords[1]
    zi = coords[2]
    if xi > 1e-14:
        crrnt_pt_angle = atan(yi / xi)
    elif xi < -1e-14:
        crrnt_pt_angle = pi + atan(yi / xi)
    else:
        crrnt_pt_angle = pi
    
    circum_wave = sin(circum_waves * crrnt_pt_angle)
    meridi_wave = cos(meridi_waves * 2 * pi * (zi - column_length/2.) / column_length)
    
    if windowing:
        circum_window = blackman_percentage((crrnt_pt_angle + pi / 2.)/(2*pi))
        meridi_window = blackman_percentage(zi/(column_length+0.1))
    else:
        circum_window = 1
        meridi_window = 1
    
    offset1 = percentage * sd.fabclass_2_umax(fab_class) * l_g * (circum_wave * meridi_wave) * cos(crrnt_pt_angle) * (circum_window * meridi_window)
    offset2 = percentage * sd.fabclass_2_umax(fab_class) * l_g * (circum_wave * meridi_wave) * sin(crrnt_pt_angle) * (circum_window * meridi_window)
    
    return offset1, offset2

def blackman_percentage(percent):
    """
    Return a Blackman windowing value for a given angle.
    Blackman is applied around a full circle.

    Parameters
    ----------
    angle : float
        An angle in radians. Must be between 0 and 2*pi.

    """
    bck = np.blackman(10000)
    crrnt_value = bck[int(np.ceil(10000 * percent))-1]
    return(crrnt_value)



def divisors(n):
    """
    Divisors of an integer.

    Return all the possible divisors for a given integer.

    Parameters
    ----------
    n : int

    Returns
    -------
    int

    """
    large_divisors = []
    for i in range(1, int(sqrt(n) + 1)):
        if n % i == 0:
            yield i
            if i * i != n:
                large_divisors.append(n / i)
    for divisor in reversed(large_divisors):
        yield divisor


def en_calcs(
        n_sides,
        r_circle,
        thickness,
        f_yield,
        lambda_flex,
        fab_class,
        arc_to_thickness
        ):
    """
    EN calculations for a polygon

    Parameters
    ----------

    """
    # Bending radius
    r_bend = arc_to_thickness * thickness
    
    # Radius of the polygon's circumscribed circle.
    r_circum = (pi * r_circle + r_bend * (n_sides * tan(np.pi/n_sides) - np.pi)) / (n_sides * sin(pi/n_sides))
    #r_circum = (pi * r_circle) / (n_sides * sin(pi / n_sides))
    
    # Width of each side.
    w_side = 2 * r_circum * sin(pi / n_sides)
    
    # Width of the corner bend half arc projection on the plane of the facet
    arc_width = r_bend * tan(np.pi / n_sides)
    
    # Flat width of each facet (excluding the bended arcs)
    facet_flat_width = w_side - 2 * arc_width
    
    # Total cross-sectional area of the bended corners
    corner_area = 2 * pi * r_bend * thickness
    
    # Central angles
    theta = 2 * pi / n_sides
    
    # Polar coordinate of ths polygon vertices on the cross-section plane.
    phii = []
    for i_index in range(n_sides):
        phii.append(i_index * theta)
    
    # Coordinates of the polygon vertices.
    x_corners = r_circum * np.cos(phii)
    y_corners = r_circum * np.sin(phii)
    
    # Cross-sectional properties
    nodes = [x_corners, y_corners]
    elem = [
        list(range(0, len(x_corners))),
        list(range(1, len(x_corners))) + [0],
        len(x_corners) * [thickness]
    ]
    
    # Cross-sectional properties
    cs_props = sd.CsProps.from_cs_sketch(sd.CsSketch(nodes, elem))
    
    # Critical stress acc. to plate theory.
    sigma_cr_plate = sd.sigma_cr_plate(thickness, facet_flat_width)
    
    # Critical load acc. to plate theory.
    n_cr_plate = pi * 2 * r_circle * thickness * sigma_cr_plate
    
    # Effective cross section area, A_eff
    a_eff = n_sides * sd.a_eff(thickness, facet_flat_width, f_yield) + corner_area
    
    # Calculate column length for the given flexural slenderness.
    column_length = lambda_flex * pi * sqrt(210000. * cs_props.moi_1 / (a_eff * f_yield))
    
    # Buckling load
    n_b_rd_plate = sd.n_b_rd(column_length, a_eff, cs_props.moi_1, f_yield, "d")
    
    # Compression resistance acc. to plate theory, EC3-1-5.
    n_pl_rd_plate = n_sides * sd.n_pl_rd(thickness, facet_flat_width, f_yield) + corner_area * f_yield
    
    # Critical stress acc. to shell theory.
    sigma_cr_shell = sd.sigma_x_rcr(thickness, r_circle, column_length)
    
    # Critical load acc. to shell theory.
    n_cr_shell = sd.n_cr_shell(thickness, r_circle, column_length)
    
    # Compression resistance acc. to shell theory, EC3-1-6.
    n_b_rd_shell = 2 * pi * r_circle * thickness * sd.sigma_x_rd(
        thickness,
        r_circle,
        column_length,
        f_yield,
        fab_quality=fab_class,
        gamma_m1=1.
    )
    
    # Material epsilon.
    epsilon = sqrt(235. / f_yield)
    
    # Classification as plate.
    p_classification = facet_flat_width / (epsilon * thickness)
    
    # Classification as tube.
    t_classification = 2 * r_circle / (epsilon ** 2 * thickness)
    
    return_dict = {
        "r_circum":r_circum,
        "w_side":w_side,
        "facet_flat_width":facet_flat_width,
        "sigma_cr_plate":sigma_cr_plate,
        "n_cr_plate":n_cr_plate,
        "n_pl_rd_plate":n_pl_rd_plate,
        "sigma_cr_shell":sigma_cr_shell,
        "n_cr_shell":n_cr_shell,
        "n_b_rd_shell":n_b_rd_shell,
        "epsilon":epsilon,
        "p_classification":p_classification,
        "t_classification":t_classification,
        "x_corners":x_corners,
        "y_corners":y_corners,
        "column_length":column_length,
        "n_b_rd_plate":n_b_rd_plate
    }
    
    return return_dict

#TODO: Update docstring for modeler
def modeler(n_sides,
            r_circle,
            thickness,
            f_yield,
            arc_to_thickness=3.,
            lambda_flex=None,
            flex_imp=0,
            imperfections=None,
            windowing=True,
            fab_class=None,
            radius_to_elsize=None,
            biased_mesh=5,
            n_eigen=None,
            submit=False,
            IDstring=None
            ):
    """
    Polygonal column modeler

    The function is responsible for modelling a polygonal column in Abaqus. Eigenvalue and RIKS analyses are performed
    depending on the given inputs.

    Parameters
    ----------
    n_sides : int
        Number of sides of the polygon cross-section.
    r_circle : float
        Radius of the equivalent circle. As equivalent is defined the circle having equal perimeter to the polygon.
    thickness : float
        Shell thickness.
    column_length : float, optional
        Length of the column. Default value is 2 times the perimeter.
    imperfections : {tuple of floats, tuple of tuples}, optional
        The imperfection variable defines the type or combination of types of the applied initial imperfections. It
        either take 2 forms:  i) a combination of explicitly defined imperfection modes. A tuple of tuples of 3 numbers
                                 should be given, where each sub-tuple is one imperfection mode and the 3 numbers are
                                 the number of circumferencial waves (integer), the number of meridional waves (integer)
                                 and the scale factor (float) respectively,
                             ii) a combination of plate-like and spillover imperfections. A tuple of 2 floats should be
                                 given, where the 2 floats are the scale factors for the plate-like and the spillover
                                 imperfection respectively.
        By default, a None value is assigned, which triggers an eigenvalue analysis, from which the first eigenmode is
        used for initial imperfections.
    windowing : bool
        Apply windowing on circumferencial imperfections.
    fab_class : {"fcA", "fcB", "fcC"}, optional
        The fabrication class used to scale the imperfections according to EC3-1-6.
    radius_to_elsize : float, optional
        Mesh density defined as the ratio of the column radius to the element size.
        Reasonable values are between 10 and 100. Default value is 20.
    f_yield : float
        Yield stress.
    n_eigen : int, optional
        The number of requested eigenvalues. The default value is 1.
    submit : bool, optional
        Flag to either submit the buckling analysis or not. It does not affect the submission of the eigenvalue
        analysis, i.e if no imperfections are defined (see "imperfections" parameter above), the eigenvalue analysis
        will be executed, calculating the requested number of eigenvalues (see "n_eigen" parameter above).
    IDstring : str
         Identification string for the model. Used in various instances as filename and folder name for the execution.

    """
    #### INPUT DEFAULTS ####
    session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
    
    # Column length
    if lambda_flex is None:
        lambda_flex = 1.
    else:
        lambda_flex = float(lambda_flex)
    
    # Amplitude of imperfection
    if fab_class is None:
        fab_class = 'fcA'
    elif not((fab_class is 'fcA') or (fab_class is 'fcB') or (fab_class is 'fcC')):
        print('Invalid fabrication class input. Choose between \'fcA\', \'fcB\' and \'fcC\' ')
    
    # Element size.
    if radius_to_elsize is None:
        radius_to_elsize = 20.
    else:
        radius_to_elsize = float(radius_to_elsize)
    
    # Number of eigenvalues requested
    if n_eigen is None:
        n_eigen = 1
    else:
        n_eigen = int(n_eigen)
    
    # ID of the current job
    if IDstring is None:
        IDstring = 'NA'
    
    #### END DEFAULT INPUTS ####
    
    #### GENERAL CALCULATIONS ####
    
    props = en_calcs(
        n_sides,
        r_circle,
        thickness,
        f_yield,
        lambda_flex,
        fab_class,
        arc_to_thickness)
    
    # Radius of the polygon's circumscribed circle.
    r_circum = props["r_circum"]
    
    # Full width of each side.
    w_side = props["w_side"]
    
    # Reduced width of each side (bended part excluded, as requested on EC31-1 and 1-5)
    facet_flat_width = props["facet_flat_width"]
    
    # Critical stress acc. to plate theory.
    sigma_cr_plate = props["sigma_cr_plate"]
    
    # Critical load acc. to plate theory.
    n_cr_plate = props["n_cr_plate"]
    
    # Compression resistance acc. to plate theory, EC3-1-5.
    n_pl_rd_plate = props["n_pl_rd_plate"]
    
    # Critical stress acc. to shell theory.
    sigma_cr_shell = props["sigma_cr_shell"]
    
    # Critical load acc. to shell theory.
    n_cr_shell = props["n_cr_shell"]
    
    # Compression resistance acc. to shell theory, EC3-1-6.
    n_b_rd_shell = props["n_b_rd_shell"]
    
    # Material epsilon.
    epsilon = props["epsilon"]
    
    # Classification as plate.
    p_classification = props["p_classification"]
    
    # Classification as tube.
    t_classification = props["t_classification"]
    
    # Polygon vertex coordinates (for the cross-section)
    x_corners = props["x_corners"]
    y_corners = props["y_corners"]
    
    # Column length
    column_length = props["column_length"]
    
    # Flexural buckling resistance
    n_b_rd_plate = props["n_b_rd_plate"]
    
    # Calculate imperfections.
    # If no waves are given, a buckling model is ran and used for GMNIA imperfections
    if imperfections is None:
        imperf_from_eigen = True
    elif isinstance(imperfections, tuple):
        imperf_from_eigen = False
        if all([isinstance(i, float) for i in imperfections]) and len(imperfections) == 2:
            imperfections = (
                (n_sides / 2., n_sides * column_length / (2 * 2 * pi * r_circle), imperfections[0]),
                (n_sides / 4., n_sides * column_length / (4 * 2 * pi * r_circle), imperfections[1])
            )
        elif all([isinstance(i, tuple) for i in imperfections]):
            for j in imperfections:
                if len(j) == 3:
                    if all([isinstance(j[0], int), isinstance(j[1], int), isinstance(j[2], float)]):
                        pass
    else:
        print("Wrong input on imperfections. Must be a tuple of 2 floats or a tuple of tuples of 2 integers and "
              "1 float. See documentation.")
    
    #### START BCKL MODEL ####
    
    # Create a new model database. This will also start a new journal file for the current session.
    Mdb()
    
    # Change the pre-existing model name
    bckl_model_name = 'bckl_model'
    mdb.models.changeKey(
        fromName = 'Model-1',
        toName = bckl_model_name
        )
    
    # Create a variable for the model
    bckl_model = mdb.models[bckl_model_name]
    
    # Create sketch
    skecth_name = 'cs_sketch'
    cs_sketch = bckl_model.ConstrainedSketch(
        name = skecth_name,
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
            radius = 3 * thickness
            )
    
    # Create the part
    part_name = 'short_column'
    p_part = bckl_model.Part(
        dimensionality = THREE_D,
        name = part_name,
        type = DEFORMABLE_BODY
        )
    
    p_part.BaseShellExtrude(
        depth = column_length,
        sketch = cs_sketch
        )
    
    # Create material
    young = 210000.
    poisson = 0.3
    el_material_name = 'elastic'
    el_material = bckl_model.Material(name=el_material_name)
    el_material.Elastic(table=((young, poisson), ))
    
    # Create shell section
    shell_section_name = 'shell_section'
    shell_section = bckl_model.HomogeneousShellSection(
        idealization = NO_IDEALIZATION, 
        integrationRule = SIMPSON,
        material = el_material_name,
        name = shell_section_name,
        numIntPts = 5,
        poissonDefinition = DEFAULT,
        preIntegrate = OFF,
        temperature = GRADIENT,
        thickness = thickness,
        thicknessField = '',
        thicknessModulus = None, 
        thicknessType = UNIFORM,
        useDensity = OFF
        )
    
    # Create a set with all the faces
    faces_set_name = 'all_faces'
    all_faces_set = p_part.Set(
        faces = p_part.faces[:],
        name = faces_set_name
        )
    
    # Assign the section to the set
    p_part.SectionAssignment(
        offset = 0.0,
        offsetField = '',
        offsetType = MIDDLE_SURFACE,
        region = all_faces_set,
        sectionName = shell_section_name,
        thicknessAssignment = FROM_SECTION
        )
    
    # Global seeding of the part according to the shell thickness
    elem_size = r_circum / radius_to_elsize
    p_part.seedPart(
        deviationFactor = 0.1, 
        minSizeFactor = 0.1,
        size = elem_size
        )
    
    if biased_mesh:
        mid_datum = p_part.DatumPlaneByPrincipalPlane(
            offset=column_length/2.,
            principalPlane=XYPLANE
            )
        p_part.PartitionFaceByDatumPlane(
            datumPlane=p_part.datums[mid_datum.id],
            faces=p_part.faces[:]
            )
        
        merid_1_ind, merid_2_ind, end_1_ind, mid_circ_ind, end_2_ind = [], [], [], [], []
        for i in p_part.edges:
            z1 = p_part.vertices[i.getVertices()[0]].pointOn[0][2]
            z2 = p_part.vertices[i.getVertices()[1]].pointOn[0][2]
            if z1 != z2 and (abs(z1)<1e-4 or abs(z1)<1e-4):
                merid_1_ind.append(i.index)
            elif z1 != z2 and(abs(z1-column_length)<1e-4 or abs(z2-column_length)<1e-4):
                merid_2_ind.append(i.index)
            elif abs(z1) < 1e-4 and abs(z2) < 1e-4:
                end_1_ind.append(i.index)
            elif abs(z1 - column_length/2.)<1e-4 and abs(z2 - column_length/2.)<1e-4:
                mid_circ_ind.append(i.index)
            elif abs(z1 - column_length)<1e-4 and abs(z2 - column_length)<1e-4:
                end_2_ind.append(i.index)
            else:
                pass
        
        merid_1 = p_part.edges[merid_1_ind[0]:merid_1_ind[0] + 1]
        merid_2 = p_part.edges[merid_2_ind[0]:merid_2_ind[0] + 1]
        end_1_circ = p_part.edges[end_1_ind[0]:end_1_ind[0] + 1]
        mid_circ = p_part.edges[mid_circ_ind[0]:mid_circ_ind[0] + 1]
        end_2_circ = p_part.edges[end_2_ind[0]:end_2_ind[0] + 1]
        
        
        for i in merid_1_ind[1:]:
            merid_1 = merid_1 + p_part.edges[i:i + 1]
        
        for i in merid_2_ind[1:]:
            merid_2 = merid_2 + p_part.edges[i:i + 1]
        
        for i in end_1_ind[1:]:
            end_1_circ = end_1_circ + p_part.edges[i:i + 1]
        
        for i in mid_circ_ind[1:]:
            mid_circ = mid_circ + p_part.edges[i:i + 1]
        
        for i in end_2_ind[1:]:
            end_2_circ = end_2_circ + p_part.edges[i:i + 1]
        
        merid_set = p_part.Set(edges=merid_1, name="merid_1_edges")
        merid_set = p_part.Set(edges=merid_2, name="merid_2_edges")
        end_1_circ_set = p_part.Set(edges=end_1_circ, name="end_1")
        mid_circ_set = p_part.Set(edges=mid_circ, name="mid_circ")
        end_2_set = p_part.Set(edges=end_2_circ, name="end_2")
        
        p_part.seedEdgeByBias(
            biasMethod=SINGLE,
            constraint=FINER,
            end2Edges=merid_1,
            maxSize=biased_mesh*elem_size,
            minSize=elem_size
            )
        p_part.seedEdgeByBias(
            biasMethod=SINGLE,
            constraint=FINER,
            end1Edges=merid_2,
            maxSize=biased_mesh*elem_size,
            minSize=elem_size
            )
        p_part.seedEdgeBySize(
            constraint=FINER, 
            deviationFactor=0.1,
            edges=end_1_circ, 
            minSizeFactor=0.1,
            size=biased_mesh*elem_size
            )
        p_part.seedEdgeBySize(
            constraint=FINER, 
            deviationFactor=0.1,
            edges=mid_circ, 
            minSizeFactor=0.1,
            size=elem_size
            )
        p_part.seedEdgeBySize(
            constraint=FINER, 
            deviationFactor=0.1,
            edges=end_2_circ, 
            minSizeFactor=0.1,
            size=biased_mesh*elem_size
            )
        p_part.setMeshControls(
            elemShape=QUAD, 
            regions=p_part.faces[:]
            )
    
    # Mesh the part
    p_part.generateMesh()
    
    # Create variable and coordinate system for the assembly
    b_assembly = bckl_model.rootAssembly
    b_assembly.DatumCsysByDefault(CARTESIAN)
    
    # Create instance
    column_instance = b_assembly.Instance(
        dependent=ON,
        name=part_name,
        part=p_part
        )
    
    # Create reference points at the ends of the column for BC couplings
    base_rp_feature = b_assembly.ReferencePoint(point=(0.0, 0.0, 0.0))
    head_rp_feature = b_assembly.ReferencePoint(point=(0.0, 0.0, column_length))

    rp_base = b_assembly.referencePoints[base_rp_feature.id]
    rp_head = b_assembly.referencePoints[head_rp_feature.id]

    # Create sets for the two reference points
    base_rp_set_name = 'base_rp'
    rp_base_set = b_assembly.Set(
        name=base_rp_set_name,
        referencePoints = (rp_base, ))

    head_rp_set_name = 'head_rp'
    rp_head_set = b_assembly.Set(
        name=head_rp_set_name,
        referencePoints = (rp_head, ))

    # Create sets for the base and the head of the column
    base_edges_set_name = 'base_edges'
    base_edges_set = b_assembly.Set(
        edges = column_instance.edges.getByBoundingBox(
            -2 * r_circum,
            -2 * r_circum,
            0,
            2 * r_circum,
            2 * r_circum,
            0),
        name = base_edges_set_name
        )

    head_edges_set_name = 'head_edges'
    head_edges_set = b_assembly.Set( 
        edges = column_instance.edges.getByBoundingBox(
            -2 * r_circum,
            -2 * r_circum,
            column_length,
            2 * r_circum,
            2 * r_circum,
            column_length),
        name = head_edges_set_name
        )

    # Create column end couplings
    # Current coupling settings restrain the shell membrain rotation
    # For free edge shell rotation, change to: ur1 = OFF, ur2 = OFF
    base_coupling_name = 'base_coupling'
    base_coupling = bckl_model.Coupling(
        controlPoint = rp_base_set,
        couplingType = KINEMATIC,
        influenceRadius = WHOLE_SURFACE,
        localCsys = None,
        name = base_coupling_name,
        surface = base_edges_set,
        u1 = ON, u2 = ON, u3 = ON, ur1 = ON, ur2 = ON, ur3 = ON
        )

    head_coupling_name = 'head_coupling'
    head_coupling = bckl_model.Coupling(
        controlPoint = rp_head_set,
        couplingType = KINEMATIC,
        influenceRadius = WHOLE_SURFACE,
        localCsys = None,
        name = head_coupling_name,
        surface = head_edges_set,
        u1 = ON, u2 = ON, u3 = ON, ur1 = ON, ur2 = ON, ur3 = ON
        )

    # Create buckling analysis step
    bckl_step_name = 'bckl-step'
    bckl_model.BuckleStep(
        name = bckl_step_name,
        numEigen = n_eigen,
        previous = 'Initial', 
        vectors = 8,
        maxIterations = 1000
        )

    # Apply concentrated load
    load_name = 'compression'
    bckl_model.ConcentratedForce(
        cf3 = -1,
        createStepName = bckl_step_name,
        distributionType = UNIFORM,
        field = '',
        localCsys = None,
        name = load_name,
        region = rp_head_set
        )

    # Hinge column base
    base_BC_name = 'fix_base'
    bckl_model.DisplacementBC(
        amplitude = UNSET,
        createStepName = 'Initial', 
        distributionType = UNIFORM,
        fieldName = '',
        localCsys = None,
        name = base_BC_name,
        region = rp_base_set,
        u1 = SET, u2 = SET, u3 = SET, ur1 = UNSET, ur2 = UNSET, ur3 = SET
        )

    # Hinge column head
    head_BC_name = 'fix_head'
    bckl_model.DisplacementBC(
        amplitude = UNSET,
        createStepName = 'Initial', 
        distributionType = UNIFORM,
        fieldName = '',
        localCsys = None,
        name = head_BC_name,
        region = rp_head_set,
        u1 = SET, u2 = SET, u3 = UNSET, ur1 = UNSET, ur2 = UNSET, ur3 = SET
        )

    # Set field output requests
    bckl_model.fieldOutputRequests['F-Output-1'].setValues(
        variables = ('U',)
        )

    #### END BCKL MODEL ####

    #### START RIKS MODEL ####
    # Riks model name
    riks_model_name = 'riks_model'

    # Copy the buckling model
    riks_mdl = mdb.Model(
        name = riks_model_name,
        objectToCopy = bckl_model
        )

    # Delete buckling step
    del riks_mdl.steps[bckl_step_name]

    # Create RIKS step
    riks_step_name = 'riks-step'
    max_increments = 20
    riks_mdl.StaticRiksStep(
        name=riks_step_name,
        previous='Initial',
        nlgeom=ON,
        maxNumInc=max_increments,
        extrapolation=LINEAR,
        initialArcInc=0.5,
        minArcInc=1e-12,
        totalArcLength=1.
        )

    # Rename the material
    material_name = "S"+"%d"%f_yield
    riks_mdl.materials.changeKey(
        fromName=el_material_name,
        toName=material_name)

    # Change to plastic material
    riks_mdl.materials[material_name].Plastic(
        table = sd.Material.plastic_table(nominal=material_name)
        )

    # Change the section material name accordingly
    riks_mdl.sections[shell_section_name].setValues(
        material=material_name,
        )

    # Create a set to act as a control point for the solver killer subroutine
    assmbl_riks = riks_mdl.rootAssembly
    instance = assmbl_riks.instances[part_name]
    killer_cp_name = 'RIKS_NODE'
    assmbl_riks.Set(
        name=killer_cp_name,
        nodes=(instance.nodes[1:2],)
        )

    # Set history output request for displacement
    disp_hist_name = "disp"
    disp_history = riks_mdl.HistoryOutputRequest(
        createStepName = riks_step_name,
        name = disp_hist_name,
        rebar = EXCLUDE,
        region = rp_head_set, 
        sectionPoints = DEFAULT,
        variables = ('U3', )
        )

    # Set history output request for load
    load_hist_name = "load"
    load_history = riks_mdl.HistoryOutputRequest(
        createStepName = riks_step_name,
        name = load_hist_name,
        rebar = EXCLUDE,
        region = rp_base_set, 
        sectionPoints = DEFAULT,
        variables = ('RF3', )
        )

    # Delete pre-existing history request: H-Output-1
    riks_mdl.historyOutputRequests.delete(['H-Output-1'])

    # Apply concentrated load
    riks_mdl.ConcentratedForce(
        cf3 = -n_b_rd_plate,
        createStepName = riks_step_name,
        distributionType = UNIFORM,
        field = '',
        localCsys = None,
        name = load_name,
        region = rp_head_set
        )

    ###### END RIKS MODEL ######

    ###### IMPERFECTIONS ######

    # Edit the coordinates of nodes to get the imperfect shape.
    # In case specific number of waves are given as input, the imperfections
    # are applied directly to the mesh by the wave functions
    # In case no imperfections are given, imperfections are introduced
    # using the first eigenmode of a buckling analysis. The amplitude of
    # this imperfection is regulated according to the

    if imperf_from_eigen:
        # Edit the keywords for the buckling model to write 'U' on file
        bckl_model.keywordBlock.synchVersions(storeNodesAndElements = False)
        bckl_model.keywordBlock.insert(at.get_block_position(bckl_model, '*End Step') - 1, '*NODE FILE\nU')

        # Create the job
        bckl_job_name = 'BCKL-'+IDstring
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
            model = bckl_model_name,
            modelPrint = OFF, 
            multiprocessingMode = DEFAULT,
            name = bckl_job_name,
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

        # Open the buckling step odb file
        eigen_data = at.fetch_eigenv(bckl_job_name, bckl_step_name, n_eigen)
        eigenvalues = eigen_data[0]
        eigen_string = eigen_data[1]

        # Calculate the proportionality of the two imperfection types
        # based on the deviation of the 1st eigenvalue to the plate critical load
        # Plate critical stress
        diff_I = (eigenvalues[0] - n_cr_plate) ** 2
        diff_II = (eigenvalues[0] - n_cr_shell) ** 2
        prop_I = diff_II / (diff_I + diff_II)
        prop_II = diff_I / (diff_I + diff_II)

        # find the maximum displacement from the buckling analysis
        bckl_odb = at.open_odb(bckl_job_name + '.odb')
        Umax = at.field_max(bckl_odb, ['U', 'Magnitude'])
        odbAccess.closeOdb(bckl_odb)

        # Plate-like imperfection half wavelength
        # perimeter divided by the number of waves
        l_g_I = 2 * pi * r_circle / (n_sides)

        # Plate-like imperfection half wavelength
        l_g_II = 2 * pi * r_circle / (2 * floor(n_sides / 4))

        # l_g is formed from l_g_I and l_g_II.
        # the contribution of each wavelength is based on the deviation of 
        # the eigenvalue analysis to the EC N_cr
        l_g = l_g_I * prop_I + l_g_II * prop_II

        # Calculate target maximum imperfection displacement (fabrication class)
        # for a length lg calculated between 1 and 2 sides (plate and spillover waves)
        u_tot = l_g * sd.fabclass_2_umax(fab_class)

        # Imperfection amplitude
        a_imp = u_tot / Umax

        # Edit the keywords for the compression riks model to include imperfections from buckling analysis
        riks_mdl.keywordBlock.synchVersions(storeNodesAndElements=False)
        riks_mdl.keywordBlock.replace(
            at.get_block_position(riks_mdl, '*step') - 1,
            '\n** ----------------------------------------------------------------\n** \n**********GEOMETRICAL '
            'IMPERFECTIONS\n*IMPERFECTION,FILE=' + bckl_job_name +',STEP=1\n1,' + str(a_imp) +'\n**'
            )
    else:
        eigenvalues = None
        for node in instance.nodes:
            coords = node.coordinates
            offset1, offset2 = calc_flex_imp(coords[2], column_length, flex_imp, 0.)
            for imp_case in imperfections:
                circum_waves, meridi_waves, scale = imp_case
                curr_offsets = calc_local_imp(coords,
                                              r_circle,
                                              column_length,
                                              circum_waves,
                                              meridi_waves,
                                              fab_class,
                                              scale,
                                              windowing=windowing)
                offset1, offset2 = offset1 + curr_offsets[0], offset2 + curr_offsets[1]

            assmbl_riks.editNode(
                nodes=node,
                offset1=offset1,
                offset2=offset2
            )

    ###### END IMPERFECTIONS ####

    ###### RIKS JOB #######

    #  Output the RIKS_NODE for GN_killer to work
    riks_mdl.keywordBlock.synchVersions(storeNodesAndElements=False)
    riks_mdl.keywordBlock.insert(
        at.get_block_position(riks_mdl, '*End Step') - 1,
            '\n*NODE FILE, GLOBAL=YES, NSET=RIKS_NODE\nU'
            )

    # Check what type of system is the modelling executed and adjust the filename of the GN_killler subroutine.
    if os.name == "posix":
        subroutine_name = '../GN_Riks_killer.f'
    else:
        subroutine_name = None

    # Create the RIKS job
    riks_job_name = 'RIKS-'+IDstring
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
        name = riks_job_name,
        nodalOutputPrecision = SINGLE, 
        numCpus = 1,
        numGPUs = 0,
        queue = None,
        resultsFormat = ODB,
        scratch = '',
        type = ANALYSIS,
        userSubroutine=subroutine_name,
        waitHours = 0,
        waitMinutes = 0
        )

    # Save the model
    mdb.saveAs(pathName=os.getcwd()+'/'+IDstring+'.cae')

    ##### END RIKS JOB ##########

    ##### PROCESSING/POST PROCESSING ########

    # Submit RIKS job
    # Assign None values to the variables hosting the analysis results, in case the user asks for no job submission.
    max_lpf, max_load, max_disp = None, None, None
    # Submit the job.
    if submit:
        riks_job.submit(consistencyChecking = OFF)
        riks_job.waitForCompletion()

        # nOpen the results database
        odb_obj = odbAccess.openOdb(path=riks_job_name + ".odb")

        # Find max values from history output
        max_lpf, max_load, max_disp = at.history_max(odb_obj, riks_step_name)

        # close database
        odbAccess.closeOdb(odb_obj)

    ##### END PROCESSING/POST PROCESSING ########

    # Create and populate an output text with model information
    out_file = open('./'+IDstring+'_info.dat', 'w')
    out_file.write('\n-GEOMETRIC CHARACTERISTICS\n')
    out_file.write('Number of corners:........................................ '+str(n_sides)+'\n')
    out_file.write('Diameter of equal perimeter circle:....................... '+str(2 * r_circle)+' [mm]\n')
    out_file.write('Diameter of the circumscribed circle:..................... '+str(2 * r_circum)+' [mm]\n')
    out_file.write('Full width of each side:.................................. '+str(w_side)+' [mm]\n')
    out_file.write('Clear flat width of each side:............................ '+str(facet_flat_width)+' [mm]\n')
    out_file.write('Perimeter:................................................ '+str(2 * pi * r_circle)+' [mm]\n')
    out_file.write('Total length:............................................. '+str(column_length/1000)+' [m]\n')
    out_file.write('Profile thickness:........................................ '+str(thickness)+' [mm]\n')
    out_file.write('Bending radius for the creases of the polygon (midline):.. '+str(3 * thickness)+' [mm]\n')
    out_file.write('Cross-sectional Area:..................................... '+str(2 * pi * r_circle * thickness)+' [mm^2]\n')
    out_file.write('\n-STRUCTURAL CHARACTERISTICS'+'\n')
    out_file.write('Yield strength:........................................... '+str(f_yield)+' [MPa]\n')
    out_file.write('Epsilon:.................................................. '+str(epsilon)+' [MPa]\n')
    out_file.write('Cross-section classification (as plate):.................. '+str(p_classification)+'\n')
    out_file.write('Cross-section classification (as tube):................... '+str(t_classification)+'\n')
    out_file.write('Plate critical stress:.................................... '+str(sigma_cr_plate)+' [MPa]\n')
    out_file.write('Shell critical stress:.................................... '+str(sigma_cr_shell[0])+' [MPa]\n')
    out_file.write('Plate critical load:...................................... '+str(n_cr_plate/1000.)+' [kN]\n')
    out_file.write('Shell critical load:...................................... '+str(n_cr_shell/1000.)+' [kN]\n')
    out_file.write('Plastic resistance, N_pl_rd:.............................. '+str(n_pl_rd_plate/1000.)+' [kN]\n')
    out_file.write('Buckling resistance, N_b_rd:.............................. '+str(n_b_rd_plate/1000.)+' [kN]\n')
    out_file.write('\n-RESULTS\n')
    if imperf_from_eigen:
        for i, e_value in eigenvalues:
            out_file.write('Eigen value '+"%02d"%(i+1)+':............................. '+str(e_value/1000)+' [kN]\n')

    if submit:
        out_file.write('Max LPF:.................................................. '+str(max_lpf)+'\n')
        out_file.write('Max load:................................................. '+str(max_load/1000.)+' [kN]\n')
        out_file.write('Displacement at max load:................................. '+str(max_disp)+' [mm]\n')

    out_file.close()

    return_string = ("%07.3f,%07.3f,%07.3f,%.5E,%.5E,%.5E"
        %(
            r_circum,
            thickness,
            column_length,
            n_pl_rd_plate,
            n_b_rd_plate,
            n_b_rd_shell
            )
        )
    if submit:
        utilization_1_5 = max_load / n_b_rd_plate
        utilization_1_6 = max_load / n_b_rd_shell
        return_string = return_string + ",%06.5f,%.5E,%07.4, %06.5f,%06.5f"%(max_lpf, max_load, max_disp, utilization_1_5, utilization_1_6)

    return return_string

def modeler_classif(n_sides,
            r_circle,
            p_classification,
            f_yield,
            arc_to_thickness=3.,
            lambda_flex=None,
            flex_imp=0,
            imperfections=None,
            windowing=True,
            fab_class=None,
            radius_to_elsize=None,
            biased_mesh=5,
            n_eigen=None,
            submit=False,
            IDstring=None
            ):
    """
    Execution of the modeler function for given classification instead of given thickness.
    See documentation od "modeler". The input parameters are identical, except the replacement of "thickness" with
    "p_classification".

    """
    # Calculate the thickness for a given slenderness.
    thickness = class_2_thickness(n_sides,
                      r_circle,
                      p_classification,
                      f_yield
                      )

    # Execute the modeler function.
    return_string = modeler(n_sides,
            r_circle,
            thickness,
            f_yield,
            arc_to_thickness=arc_to_thickness,
            lambda_flex=lambda_flex,
            flex_imp=flex_imp,
            imperfections=imperfections,
            windowing=windowing,
            fab_class=fab_class,
            radius_to_elsize=radius_to_elsize,
            biased_mesh=biased_mesh,
            n_eigen=n_eigen,
            submit=submit,
            IDstring=IDstring
            )

    return return_string


def results_from_odb(filename=None):
    """
    Get results from polygonal odb.

    Parameters
    ----------
    odb_file : str
        Relative path/filename of the odb file.
    requests : dict, optional
        Results requested from odb. Default is all.
    headers : bool of list of strings, optional
        If given, returns only the headers of the results
    Returns
    -------
    str

    """
    if filename is None:
        filename = at.find_odb_in_cwd()

    # Open odb database
    # odb = at.open_odb(filename)
    odb = odbAccess.openOdb(path=filename, readOnly=True)

    # Find max values from history output
    max_lpf, max_load, max_disp = at.history_max(odb, "riks-step")

    # Close file
    odbAccess.closeOdb(odb)

    # Return
    return ("%08.6f"%(max_lpf)+" %013.3f"%(max_load)+" %08.3f"%(max_disp))
