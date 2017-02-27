# Input variables for the for the semi-closed polygonal column model
from math import pi

def polygon_input(*arg):
    # Number of corners
    n = 6
    
    # Diameter of the circumscribed circlein mm
    d = 1000.
    
    # Cross-section slenderness defined by the prescribed circle lambda1=(d/(t^2*epsilon))
    cs_slenderness = 120.
    
    # Member slenderness for overall column buckling lambda2= sqrt(A*fy/Ncr)
    mb_slenderness = 1.
    
    # Radius if the bended corners of the polygon given as a ratio to the thickness r=rcoef*t
    # It regards to the bends of the polygon. The arc radious of the lips' bends is half this value
    rcoef = 6.
    
    # Number of elements along the arc length of the bended corners (to be changed: the number of elements have to be according to the global seeding)
    nbend = 3
    
    # Thickness of the gusset plates given as a ratio to the profile thickness tgusset=t_ratio*t
    t_ratio = 1.2
    
    # Yield strength in MPa. Used for epsilon, not for the modelling material properties
    fy = 355.
    
    # Young's modulus
    E_young = 210000.
    
    # Bolt spacing given as a ratio to the prescribed circle diameter, b=s/d.
    b = 3
    
    # Imperfection factor for overall bowing of the column u1=l/impamp
    flx_imp = 100
    
    # Direction angle for the overall bowing, 0 is on global y axis in rads
    theta_bow = pi/2
    
    # Imperfection factor for distortional imperfections u2=s/dist_imp
    dist_imp = 50
    
    # Bolt diameter in mm
    M_bolt = 16
    
    # Clearence from the washer to the edge of the lip and the start of the bending arc in mm
    clearence = 3
	
	# Maximum number of increments for the RIKS solver
	
    return n, d, cs_slenderness, mb_slenderness, rcoef, nbend, t_ratio, fy, E_young, b, flx_imp, theta_bow, dist_imp, M_bolt, clearence
	