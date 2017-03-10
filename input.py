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
                                                       'clearence',
                                                       'max_RIKS_increments',
                                                       'elem_size',
                                                       'classify_as'
                                                       ])
    #--------------------------------------------------------------------------------------------------------------
    # The first 8 input parameters in order define the filename
    # Number of corners
    try:
        n = sys.argv[-8]
    except:
        n = 6
    
    # Diameter of the circumscribed circlein mm
    try:
        d = int(sys.argv[-7])
    except:
        d = 500.
    
    # Bolt spacing given as a ratio to the prescribed circle diameter, b=s/d.
    try:
        b = int(sys.argv[-6])
    except:
        b = 1
    
    # Classification as tube or plates
    #class_type = 'tube'
    class_type = 'plate'
    
    # Cross-section slenderness defined by the prescribed circle lambda1=(d/(t^2*epsilon))
    try:
        cs_slenderness = int(sys.argv[-5])
    except:
        cs_slenderness = 60.
    
    # Member slenderness for overall column buckling lambda2= sqrt(A*fy/Ncr)
    try:
        mb_slenderness = int(sys.argv[-4])/10
    except:
        mb_slenderness = 0.2
    
    # Yield strength in MPa. Used for epsilon, not for the modelling material properties
    try:
        fy = int(sys.argv[-3])
    except:
        fy = 355.
    
    # Imperfection factor for overall bowing of the column u1=l/impamp
    try:
        bow_imp = int(sys.argv[-2])
    except:
        bow_imp = 250
    
    # Imperfection factor for distortional imperfections u2=s/dist_imp
    try:
        dist_imp = int(sys.argv[-1])
    except:
        dist_imp = 300
    
    #--------------------------------------------------------------------------------------------------------------
    # The following parameters are not included in the filename.
    # Different parent directory should be used to avoid matching filenames
    
    # Radius of the bended corners of the polygon given as a ratio to the thickness r=rcoef*t
    # It regards to the bends of the polygon. The arc radious of the lips' bends is half this value
    rcoef = 6.
    
    # Thickness of the gusset plates given as a ratio to the profile thickness tgusset=t_ratio*t
    t_ratio = 2
    
    # Young's modulus
    E_young = 210000.
    
    # Direction angle for the overall bowing, 0 is on global y axis in rads
    theta_bow = pi/2
    
    # Bolt diameter in mm
    M_bolt = 16
    
    # Clearence from the washer to the edge of the lip and the start of the bending arc in mm
    clearence = 3
    
    # Maximum number of increments for the RIKS solver
    max_inc = 100
    
    # Mesh seeding size
    seed_size = d/15
    
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
                      clearence,
                      max_inc,
                      seed_size,
                      class_type
                      )
