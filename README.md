# Polygoner

Script creating sectors of polygonal profiles.

# Use
Just edit the input variables and execute the script in abaqus. The input variables are currently within the first 100 lines of the script Polygoner.py. 


## The input variables are as follows:

 Number of corners
n = 6

 Diameter of prescribed circle in mm
d = 1000.

 Cross-section slenderness defined by the prescribed circle lambda1=(d/(t^2\*epsilon))
cs_slenderness = 120.

 Member slenderness for overall column buckling lambda2= sqrt(A\*fy/Ncr)
mb_slenderness = 1.

 Radius if the bended corners of the polygon given as a ratio to the thickness r=rcoef\*t
 It regards to the bends of the polygon. The arc radious of the lips' bends is half this value
rcoef = 6.

 Number of elements along the arc length of the bended corners
nbend = 3

 length of the lips given as a ratio to the prescribed circle diameter l=d\*l_ratio
l_ratio = 0.14

 Thickness of the gusset plates given as a ratio to the profile thickness tgusset=t_ratio\*t
t_ratio = 1.2

 Yield strength in MPa. Used for epsilon, not for the modelling material properties
fy = 355.

 Young's modulus
E_young = 210000.

 Bolt spacing given as a ratio to the prescribed circle diameter, b=s/d.
b = 1

 Imperfection factor for overall bowing of the column u1=l/impamp
flx_imp = 250

 Direction angle for the overall bowing, 0 is on global y axis in rads
theta = pi/2

 Imperfection factor for distortional imperfections u2=s/dist_imp
dist_imp = 250

 Bolt diameter in mm
M_bolt = 16

 Clearence from the washer to the edge of the lip and the start of the bending arc in mm
clearence = 3

# Acknowledgements
This code was used for masters thesis and publication as follows:
"Finite Element Modelling and Parametric Studies of Semi-Closed Thin-Walled Steel Polygonal Columns For The Application on Steel Lattice Towers" https://doi.org/10.32783/csid-jid.v2i2.65
