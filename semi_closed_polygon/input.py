# Input variables for the for the semi-closed polygonal column model
from math import pi
import collections

def polygon_input(*arg):
    # Initiate a named tuple for the input parameters
    Parameters = collections.namedtuple('Parameters', ['n_sides',
                                                       'diameter',
                                                       'classification',
                                                       'slenderness',
                                                       'bending_arc_radius',
                                                       'gusset_thickness_ratio',
                                                       'yield_stress',
                                                       'Youngs_modulus',
                                                       'bolt_spacing',
                                                       'bow_imperfections',
                                                       'imperfections_angle',
                                                       'distortional_imperfections',
                                                       'bolt_diameter',
                                                       'bend_clearence',
                                                       'edge_clearence',
                                                       'max_RIKS_increments',
                                                       'elem_size',
                                                       'classify_as',
                                                       'first_bolt',
                                                       'joint_bolts',
                                                       'joint_length',
                                                       'load_cross_section',
                                                       'hard_contact'
                                                       ])
    #--------------------------------------------------------------------------------------------------------------
    # The first 8 input parameters in order define the filename
    # Number of corners

    n = 9
    
    # Diameter of the circumscribed circle in mm

    d = 500.
    
    # Bolt spacing given as a ratio to the prescribed circle diameter, b=s/d.

    b = 1.
    
    # Classification as tube or plates
    #class_type = 'tube'
    class_type = 'plate'
    
    # Cross-section slenderness defined by the prescribed circle lambda1=(d/(t^2*epsilon))

    cs_slenderness = 42.
    
    # Member slenderness for overall column buckling lambda2= sqrt(A*fy/Ncr)

    mb_slenderness = 0.2
    
    # Yield strength in MPa. Used for epsilon, not for the modelling material properties

    fy = 381.
    
    # Imperfection factor for overall bowing of the column u1=l/impamp
	#  Taken from EC3-1-1 table 5.1 elastic analysis curve c

    bow_imp = 200
    
    # Imperfection factor for distortional imperfections u2=s/dist_imp
	# Takesn from EC3-1-5 annex C table C.2

    dist_imp = 200
    
    #--------------------------------------------------------------------------------------------------------------
    # The following parameters are not included in the filename.
    # Different parent directory should be used to avoid matching filenames
    
    # Radius of the bended corners of the polygon given as a ratio to the thickness r=rcoef*t
    # It regards to the bends of the polygon. The arc radious of the lips' bends is half this value
    rcoef = 6.
    
    # Thickness of the gusset plates given as a ratio to the profile thickness tgusset=t_ratio*t
    # For equal cross sectional area between gusset plate and polygonal, t_ratio = (2. * pi) / 3.
    t_ratio = 1.1 * (2. * pi) / 3.
    
    # Young's modulus
    E_young = 210000.
    
    # Direction angle for the overall bowing, 0 is on global y axis in rads
    theta_bow = pi/2
    
    # Bolt diameter in mm
    M_bolt = 16
    
    # Clearence from the washer to the start of the bending arc in mm
    bend_clearence = 3
     
    # Clearence from the washer to the edge of the lip in mm
    edge_clearence = 15
   
    # Maximum number of increments for the RIKS solver
    max_inc = 40
    
    # Mesh seeding size
    seed_size = d/30

    # Number of bolts on the connection
    joint_bolts = 5
    	
    # Distance of first bolt adjascent to the connection given as a ratio to the cs plate thickness
    first_bolt = 35.

    # Length of the joint as a ratio to the diameter
    joint_length = 0.8

    # Is the load applied on the gusset plate or on the entire cross-cross section? (0 for gusset, 1 for entire section)
    load_cross_section = 1

    # Is hard contact interaction applied between the lips? (0 for no, 1 for yes)
    hard_contact = 1

    # Collect and return the input values in a named tuple
    
    return Parameters(n,
                      d,
                      cs_slenderness,
                      mb_slenderness,
                      rcoef,
                      t_ratio,
                      fy,
                      E_young,
                      b,
                      bow_imp,
                      theta_bow,
                      dist_imp,
                      M_bolt,
                      bend_clearence,
                      edge_clearence,
                      max_inc,
                      seed_size,
                      class_type,
                      first_bolt,
                      joint_bolts,
                      joint_length,
                      load_cross_section,
                      hard_contact
                      )
