### Comments on the present script
# Updated 18/04/2017
# This script describes the modelling of a polygonal semi-closed column.
# The following files are required to be in the same directory:
#  input.py			-> file that lists the input parameters for the model
#  abq_toolset.py	-> additional python tools for abaqus
#  GN_Riks_killer	-> Optional. Fortran code that kills the analysys under certain conditions. To neglect, comment the corresponding line on the riks job.
# For the GN_Riks_killer to operate on windows, the extension must be .for or .obj. To run on linux, it must be .f. Change the filename and the job creation accordingly
'''
'''

# Calculate xy of nodes for a given polygonal profile
# Returns points for the entire profile (1st and 2nd returned values)
# and points for a single sector (3rd and 4th returned values)
def polygon_sector(n, R, t, tg, rbend, nbend, l_lip):
    # Angle corresponding to one face of the polygon
    theta = 2*math.pi/n;
    
    # Angles of radii (measured from x-axis)
    phi=np.linspace(5*math.pi/6, math.pi/6, n/3+1)
    
    # xy coords of the polygon's corners
    x = R*np.cos(phi);
    y = R*np.sin(phi);
    
    ## Bends
    
    # Distance between bending centre and corner
    lc = rbend/np.cos(theta/2)
    
    # Centers of bending arcs
    xc  = x[1:-1] - lc*np.cos(phi[1:-1])
    yc  = y[1:-1] - lc*np.sin(phi[1:-1])
    
    # Bending arc angle
    theta_b = math.pi - theta
    
    # Angles of the edges' midlines (measured from x-axis)
    phi_mids = phi[0:-1] - theta/2 
    
    # xy coords of the arc's points
    xarc = [[0 for j in range(nbend+1)] for i in range(int(n/3 -1))]
    yarc = [[0 for j in range(nbend+1)] for i in range(int(n/3 -1))] 
    for i in range(int(n/3 -1)):
        for j in range(nbend+1):
            xarc[i][j] = xc[i] + rbend*np.cos(phi_mids[i]-(j)*(theta/nbend))
            yarc[i][j] = yc[i] + rbend*np.sin(phi_mids[i]-(j)*(theta/nbend))
    
    ## Start-end extensions
    # Bending radius
    rs = rbend/2
    xcs = [0, 0]
    ycs = [0, 0]
    
    # First bend
    v1 = phi_mids[0]-math.pi/2
    v2 = (phi[0]+phi_mids[0]-math.pi/2)/2
    l1 = (t+tg)/(2*np.cos(phi[0]-phi_mids[0]))
    l2 = rs/np.sin(v2-phi_mids[0]+math.pi/2)
    x1 = x[0]+l1*np.cos(v1)
    y1 = y[0]+l1*np.sin(v1)
    
    # First bend centre coords
    xcs[0] = x1+l2*np.cos(v2)
    ycs[0] = y1+l2*np.sin(v2)
    
    # Last bend
    v1 = phi_mids[-1]+math.pi/2
    v2 = (v1+phi[-1])/2
    l1 = (t+tg)/(2*np.cos(v1-phi[-1]-math.pi/2))
    l2 = rs/np.sin(v2-phi[-1])
    x1 = x[-1]+l1*np.cos(v1)
    y1 = y[-1]+l1*np.sin(v1)
    
    # Last bend centre coords
    xcs[1] = x1+l2*np.cos(v2)
    ycs[1] = y1+l2*np.sin(v2)
    
    # First and last bend arc points coords
    xsarc = [[0 for j in range(nbend+1)] for j in [0,1]]
    ysarc = [[0 for j in range(nbend+1)] for j in [0,1]] 
    for j in range(nbend+1):
        xsarc[0][j] = xcs[0] + rs*np.cos(4*math.pi/3+(j)*((phi_mids[0]-math.pi/3)/nbend))
        ysarc[0][j] = ycs[0] + rs*np.sin(4*math.pi/3+(j)*((phi_mids[0]-math.pi/3)/nbend))
        xsarc[1][j] = xcs[1] + rs*np.cos(phi_mids[-1]+math.pi+(j)*((phi[-1]+math.pi/2-phi_mids[-1])/nbend))
        ysarc[1][j] = ycs[1] + rs*np.sin(phi_mids[-1]+math.pi+(j)*((phi[-1]+math.pi/2-phi_mids[-1])/nbend))
    
    
    ## Points of the lips
    
    # Lip length according to bolt washer diameter
    
    # First lip
    xstart = [xsarc[0][0] + l_lip*np.cos(phi[0]), xsarc[0][0] + l_lip*np.cos(phi[0])/2]
    ystart = [ysarc[0][0] + l_lip*np.sin(phi[0]), ysarc[0][0] + l_lip*np.sin(phi[0])/2]
    
    
    # Last point
    xend = [xsarc[1][-1] + l_lip*np.cos(phi[-1])/2, xsarc[1][-1] + l_lip*np.cos(phi[-1])]
    yend = [ysarc[1][-1] + l_lip*np.sin(phi[-1])/2, ysarc[1][-1] + l_lip*np.sin(phi[-1])]
    
    ## Collect the x, y values in a sorted 2xn array
    xarcs, yarcs=[],[]
    for i in range(len(phi)-2):
        xarcs=xarcs+xarc[i][:]
        yarcs=yarcs+yarc[i][:]
    
    x_sector = xstart+xsarc[0][:]+xarcs[:]+xsarc[1][:]+xend
    y_sector = ystart+ysarc[0][:]+yarcs[:]+ysarc[1][:]+yend
    
    # Copy-rotate the points of the first sector to create the entire CS
    # Rotation matrix
    Rmat = np.array([[math.cos(-2*math.pi/3), -math.sin(-2*math.pi/3)], [math.sin(-2*math.pi/3), math.cos(-2*math.pi/3)]])
    
    # Dot multiply matrices
    coord1 = np.array([x_sector, y_sector])
    coord2 = Rmat.dot(coord1)
    coord3 = Rmat.dot(coord2)
    
    # Concatenate into a single xy array
    x_cs = np.concatenate([coord1[0], coord2[0], coord3[0]])
    y_cs = np.concatenate([coord1[1], coord2[1], coord3[1]])
    
    # Return matrices
    return x_cs, y_cs, x_sector, y_sector


# # Import necessary libraries --------------------------------------------------------------------

import numpy as np
import sys
import os
import meshEdit
import abq_toolset as xtr
import steel_tools as stlt
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
from shutil import copyfile
import input
input = reload(input)
import odbAccess
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

# Fetch the input variables from the file input.py
parameters = input.polygon_input()


# Replace input parameters from the input file with parameters given on the command line
# example of shell arguments:
# -- 9 500 1 42 80 381 200 200
# gives: 9 sided polygon, 500 mm diameter of equal perimeter tube, 1*diameters bolt spacing, lambda = 0.8, fy = 381, flex_imp = l/200, loc_imp = b/200

try:
    parameters = parameters._replace(n_sides = int(sys.argv[-8]))
except:
    pass

# Diameter of the circumscribed circlein mm
try:
    parameters = parameters._replace(diameter = int(sys.argv[-7]))
except:
    pass

# Bolt spacing given as a ratio to the prescribed circle diameter, b=s/d.
try:
    parameters = parameters._replace(bolt_spacing = int(sys.argv[-6]) / 10.)
except:
    pass

# Cross-section slenderness defined by the prescribed circle lambda1=(d/(t^2*epsilon))
try:
    parameters = parameters._replace(classification = int(sys.argv[-5]))
except:
    pass

# Member slenderness for overall column buckling lambda2= sqrt(A*fy/Ncr)
try:
    parameters = parameters._replace(slenderness=int(sys.argv[-4]) / 100.)
except:
    pass

# Yield strength in MPa. Used for epsilon, not for the modelling material properties
try:
    parameters = parameters._replace(yield_stress = int(sys.argv[-3]))
except:
    pass

# Imperfection factor for overall bowing of the column u1=l/impamp
try:
    parameters = parameters._replace(bow_imperfections = int(sys.argv[-2]))
except:
    pass

# Imperfection factor for distortional imperfections u2=s/dist_imp
try:
    parameters = parameters._replace(distortional_imperfections = int(sys.argv[-1]))
except:
    pass




# Model specific ID string. Used for save filename and for the jobnames
# the ID string has the following structure (example given for filename 6-1000-30-120-100-355-250-250):
#
# polygon number of sides | circumcircle diameter | 10*bolt spacing | cs_slenderness | 100*flexural slenderness | fy      | bowing imperfection | distortional imperfection
# 6 (hexagon)             | d = 1000 mm           | Sb = 3*d        | d/(t*e^2) =120 | lambda = 100             | 355 MPa | a = l/250           | a = Sb/250
#
if parameters.classification < 100:
    class_string = '0' + str(int(parameters.classification))
else:
    class_string = str(int(parameters.classification))

if parameters.classify_as == 'plate':
    class_string = 'P' + class_string
else:
    class_string = 'T' + class_string

if 100*parameters.slenderness < 100:
    slend_string = '0' + str(int(100*parameters.slenderness))
else:
    slend_string = str(int(100*parameters.slenderness))

IDstring = str(int(parameters.n_sides))+'-'+\
           str(int(parameters.diameter))+'-'+\
           str(int(parameters.bolt_spacing * 10))+'-'+\
           class_string+'-'+\
           slend_string+'-'+\
           str(int(parameters.yield_stress))+'-'+\
           str(int(parameters.bow_imperfections))+'-'+\
           str(int(parameters.distortional_imperfections))

# Make a new subdirectory for the current session
os.mkdir(IDstring)

# Copy necessary files to the new directory
copyfile('abq_toolset.py', './'+IDstring+'/abq_toolset.py')
copyfile('input.py', './'+IDstring+'/input.py')
copyfile('polygoner.py', './'+IDstring+'/polygoner.py')
copyfile('GN_Riks_killer.f', './'+IDstring+'/GN_Riks_killer.f')

# Change the working directory
os.chdir('./'+IDstring)

# Calculated model characteristics -----------------------------------------------------------------------------------------------------------
# The diameter given as input corresponds to a circle with circumferencial length equal to the perimeter of the polygon
# The polygon is created based on the circumscribed circle.
# Diameter of the circumscribed circle:
d_circumscribed = parameters.diameter
R = d_circumscribed/2

# Width of each side
face_width = d_circumscribed * sin(pi/parameters.n_sides)

# Epsilon value as given in EC3
epsilon = sqrt(235./parameters.yield_stress)

# Thickness of the profile plate
# calculated based on EC3-1-1 for a tube of the same diameter
if parameters.classify_as == 'plate':
    cs_thickness = face_width / (epsilon * parameters.classification)
else:
    cs_thickness = d_circumscribed / (epsilon**2 * parameters.classification)

# Calcilate classification number both as plate and tube
p_classification = face_width / (cs_thickness * epsilon)
if p_classification > 42:
    p_class = 'Class 4'
else:
    p_class = 'Class 1/2/3'

t_classification = d_circumscribed / (cs_thickness * epsilon ** 2)


# Thickness of the gusset plate
gusset_thickness = (parameters.gusset_thickness_ratio * cs_thickness)

# Diameter of the washer for the given bolt
d_washer = stlt.bolt2washer(parameters.bolt_diameter)

# Calculate lip length
l_lip = d_washer + parameters.bend_clearence + parameters.edge_clearence

# Create a new model ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

static_model_name = 'IMP'

stc_mdl = mdb.Model(modelType=STANDARD_EXPLICIT, name=static_model_name)

## Save the model -------------------------------------------------------------------------------------------------------
mdb.saveAs(pathName=os.getcwd()+'/'+IDstring+'.cae')


# Calculate cross section coordinates for one sector --------------------------------------------------------------------------------

# Number of elements along the arc length of the bended corners
r_bend=parameters.bending_arc_radius*cs_thickness
n_bend, null = divmod(pi*(1-2/parameters.n_sides)*r_bend, (parameters.elem_size))
n_bend = int(n_bend)

x_cs, y_cs, x_sector, y_sector = polygon_sector(parameters.n_sides,
                                                    R,
                                                     cs_thickness,
                                                    gusset_thickness,
                                                    r_bend,
                                                    n_bend,
                                                    l_lip)

# Build the node and connectivity matrices
# number of nodes
nnodes = len(x_cs)

# Gather nodes and elements
# Three very thin elements are added to connect the lips and form a closed cs
coord = [x_cs, y_cs]
ends = [range(nnodes), range(1,nnodes)+[0], [cs_thickness]*(nnodes/3-1)+[0.01]+[cs_thickness]*(nnodes/3-1)+[0.01]+[cs_thickness]*(nnodes/3-1)+[0.01]]

# Calculate cross sectional properties
Area, xc, yc, Ix, Iy, Ixy, I1, I2, theta_principal = stlt.cs_prop(coord, ends)

# Bending stiffness E*I
E_I = parameters.Youngs_modulus * I2

# Distance of the first bolt to the end of the gusset plate
end_bolts_dist = parameters.first_bolt * cs_thickness

# Spacing between two middle bolts
bolts_dist = d_circumscribed * parameters.bolt_spacing

# Calculate the minimum column length. This length will be increased to fit a certain number of bolts
min_length = parameters.slenderness*pi*sqrt(parameters.Youngs_modulus * I2 / (Area * parameters.yield_stress))

# The length of the joint region
l_joint = d_circumscribed * parameters.joint_length

# Calculate the number of bolts that can fit between the joints
n_of_bolts = int(ceil((min_length - 2 * (l_joint - end_bolts_dist)) / bolts_dist) - 1)

# Final length, increased to fit a number of bolts with specific bolt spacing
l_tot = 2 * (l_joint + end_bolts_dist) + (n_of_bolts + 1) * bolts_dist

# Design plastic resistance
N_pl_rd = parameters.yield_stress * Area

# Design buckliing resistance
# Euler critical load
N_cr = (pi ** 2 * E_I) / (l_tot ** 2)

# Member slenderness
lambda_final = sqrt(N_pl_rd / N_cr)

# Reduction factor chi
a_imp_fact = 0.49
phi_capital = (1 + a_imp_fact * (lambda_final - 0.2) + lambda_final ** 2) / 2
chi_flex = 1 / (phi_capital + sqrt(phi_capital ** 2 - lambda_final ** 2))

# Aeff is calculated assuming uniform compression on the sectors.
# Not true when the load is applies on the gusset only.

psi = 1
if psi > 0:
    kapa_sigma = 8.2 / (1.05 + psi)
elif psi > -3 and psi < -1:
    kapa_sigma = 5.98 * (1 - psi) ** 2*d_washer
else:
    kapa_sigma = 7.81 - 6.29 * psi + 9.78 * psi ** 2

lambda_p = p_classification / (28.4 * sqrt(kapa_sigma))
if lambda_p > 0.673 and p_classification > 42:
    rho = (lambda_p - 0.055 * (3 + psi)) / lambda_p ** 2
else:
    rho = 1.


# Buckling resistance
N_b_rd = chi_flex * rho * N_pl_rd

# Number of spaces between the inner bolts on the joint
joint_spaces = parameters.joint_bolts - 1

# Create a text file with the model data
out_file = open('./model_info.dat', 'w')
out_file.write('\n-GEOMETRIC CHARACTERISTICS\n')
out_file.write('Number of corners:....................................................... '+str(parameters.n_sides)+'\n')
out_file.write('Diameter of equal perimeter:............................................. '+str(parameters.diameter)+' [mm]\n')
out_file.write('Diameter of the circumscribed circle:.................................... '+str(d_circumscribed)+' [mm]\n')
out_file.write('Total length:............................................................ '+str(l_tot/1000)+' [m]\n')
out_file.write('Profile thickness:....................................................... '+str(cs_thickness)+' [mm]\n')
out_file.write('Gusset thickness:........................................................ '+str(gusset_thickness)+' [mm]\n')
out_file.write('Lip width:............................................................... '+str(l_lip)+' [mm]\n')
out_file.write('Bending radius for the creases of the polygon (midline):................. '+str(r_bend)+' [mm]\n')
out_file.write('Bolt spacing for inner bolts:............................................ '+str(bolts_dist)+' [mm]\n')
out_file.write('    Which is:............................................................ '+str(parameters.bolt_spacing)+' * Diameter\n')
out_file.write('Bolt spacing for bolt adjacent to joint:................................. '+str(end_bolts_dist)+' [mm]\n')
out_file.write('    Which is:............................................................ '+str(parameters.first_bolt)+' * t\n')
out_file.write('Total number of bolts between joints for each pair of lips:.............. '+str(n_of_bolts)+'\n')
out_file.write('Length of the joint:..................................................... '+str(l_joint)+' [mm]\n')
out_file.write('    Which is:............................................................ '+str(parameters.joint_length)+' * Diameter \n')
out_file.write('Number of bolts for each fin of the gusset plate:........................ '+str(parameters.joint_bolts)+'\n')
out_file.write('\n-STRUCTURAL CHARACTERISTICS'+'\n')
out_file.write('Cross-sectional Area:.................................................... '+str(Area)+' [mm^2]\n')
out_file.write('Moment of inertia:....................................................... '+str(I2)+' [mm^4]\n')
out_file.write('Yield strength:.......................................................... '+str(parameters.yield_stress)+' [MPa]\n')
out_file.write('Cross-section classification (as plate):................................. '+str(p_classification)+'\n')
out_file.write('    Which is:............................................................ '+p_class+' * Diameter \n')
out_file.write('Cross-section classification (as tube):.................................. '+str(t_classification)+'\n')
out_file.write('Critical load, N_cr:..................................................... '+str(N_cr/1000)+' [kN]\n')
out_file.write('Plastic resistance, N_pl_rd:............................................. '+str(N_pl_rd/1000)+' [kN]\n')
out_file.write('Buckling resistance, N_b_rd:............................................. '+str(N_b_rd/1000)+' [kN]\n')
out_file.write('Member slenderness, lambda:.............................................. '+str(lambda_final)+'\n')
out_file.write('Flexural buckling reduction factor, chi:................................. '+str(chi_flex)+'\n')
out_file.write('Plate slenderness, lambda_p:............................................. '+str(lambda_p)+'\n')
out_file.write('Effective cross section, rho:............................................ '+str(rho)+'\n')
out_file.write('Input target slenderness:................................................ '+str(parameters.slenderness)+'\n')
out_file.write('\n-MODEL CHARACTERISTICS'+'\n')
out_file.write('Flexural buckling bow imperfections:..................................... l/'+str(parameters.bow_imperfections)+'\n')
out_file.write('Plate imperfections:..................................................... b/'+str(parameters.distortional_imperfections)+'\n')
out_file.write('Target element size:..................................................... '+str(parameters.elem_size)+' [mm]\n')
out_file.write('Load applied on: the entire cs ("True"), the gusset plate ("False"):..... '+str(bool(parameters.load_cross_section))+'\n')
out_file.close()

# Create Parts ----------------------------------------------------------------------------------------------------------------------

# Sector

# -Profile sketch for sector
sector_sketch = stc_mdl.ConstrainedSketch(
    name='sector',
    sheetSize=1200.0
    )

# -Sketch sector lines
for n in range(len(x_sector)-1):
    sector_sketch.Line(
        point1=(x_sector[n], y_sector[n]), 
        point2=(x_sector[n+1], y_sector[n+1])
        )

# -Extrude sector part

sector_part = stc_mdl.Part(
    dimensionality=THREE_D,
    name='sector',
    type=DEFORMABLE_BODY
    )

sector_part.BaseShellExtrude(
    depth=l_tot,
    sketch=sector_sketch
    )

# Calculate bolt positions

# Washer diameter
d_washer = 30

# -Distance on the width
#bolts_w = l_lip/2
bolts_w = d_washer / 2 + parameters.edge_clearence

# -Distances on the length

bolts_z1 = np.concatenate([
                          [bolts_w],
                          bolts_w + ((l_joint - 2 * bolts_w) / joint_spaces) * np.linspace(1, joint_spaces, joint_spaces)
                          ])

bolts_z2 = np.concatenate([
                          [l_joint + end_bolts_dist],
                          l_joint + end_bolts_dist + (bolts_dist * np.linspace(1, n_of_bolts+1, n_of_bolts+1))
                          ])

bolts_z3 = bolts_z1 + (l_tot - l_joint)

bolts_z = np.concatenate([bolts_z1, bolts_z2, bolts_z3])

# Initiate list to store datum planes
datum_p=[]

# Make holes and datum planes
                    
for o in range(int(bolts_z.shape[0])):
                            
    sector_part.HoleBlindFromEdges(
        depth=1.0,
        diameter=d_washer,
        distance1=bolts_z[o],
        distance2=bolts_w,
        edge1=sector_part.edges.getClosest(coordinates=((x_sector[1],y_sector[1], 0),))[0][0],
        edge2=sector_part.edges.getClosest(coordinates=((x_sector[0], y_sector[0], 1),))[0][0],
        plane=sector_part.faces.getClosest(coordinates=((x_sector[0], y_sector[0], 0),))[0][0],
        planeSide=SIDE1
        )
    
    sector_part.HoleBlindFromEdges(
        depth=1.0,
        diameter=d_washer,
        distance1=bolts_z[o],
        distance2=bolts_w,
        edge1=sector_part.edges.getClosest(coordinates=((x_sector[-2], y_sector[-2], 0),))[0][0],
        edge2=sector_part.edges.getClosest(coordinates=((x_sector[-1], y_sector[-1], 1),))[0][0],
        plane=sector_part.faces.getClosest(coordinates=((x_sector[-1], y_sector[-1], 0),))[0][0],
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
n_dat = len(datum_p)

# cut all the faces using the datum planes
                    
for o in range((n_dat-2)):
    sector_part.PartitionFaceByDatumPlane(
        datumPlane=sector_part.datums.items()[o+1][1],
        faces=sector_part.faces[:]
        )

# Gusset

# -Profile sketch for gusset
gusset_sketch=stc_mdl.ConstrainedSketch(name='__profile__', sheetSize=1200.0)

# -Sketch gusset lines
# First point of the first sector
x0 = x_sector[0]
y0 = y_sector[0]

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
gusset_part=stc_mdl.Part(
    dimensionality=THREE_D,
    name='gusset',
    type=DEFORMABLE_BODY
    )

gusset_part.BaseShellExtrude(
    depth=l_joint,
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

# Sketch for the elliptical cut

elliptic_sketch = stc_mdl.ConstrainedSketch(
    name='elliptic_cut_sketch', 
    sheetSize=l_joint,
    transform=gusset_part.MakeSketchTransform(
        sketchPlane=gusset_part.faces.findAt((0.0, -1, 0), ),
        sketchPlaneSide=SIDE1, 
        sketchUpEdge=gusset_part.edges.findAt((0.0, -1, 0), ),
        sketchOrientation=RIGHT,
        origin=(0.0, -gp2[1] / 2, l_joint / 2)
        )
    )

gusset_part.projectReferencesOntoSketch(
    filter=COPLANAR_EDGES,
    sketch=elliptic_sketch
    )

elliptic_sketch.Spline(
    points=(
        (-l_joint / 4, -gp2[1] / 2 + gp2[1]),
        (-l_joint / 4, -gp2[1] / 3 + gp2[1]),
        (l_joint / 2, gp2[1] /2 + l_lip + gp2[1])
        )
    )

elliptic_sketch.Line(
    point1=(-l_joint / 4, -gp2[1] / 2 + gp2[1]),
    point2=(l_joint / 2, -gp2[1] /2 + gp2[1])
    )

elliptic_sketch.Line(
    point1=(l_joint / 2, -gp2[1] /2 + gp2[1]),
    point2=(l_joint / 2, gp2[1] /2 + l_lip + gp2[1])
    )

elliptic_sketch.ConstructionLine(
    angle=0.0,
    point1=(l_joint / 2, -gp2[1] /2 + gp2[1])
    )

gusset_part.CutRevolve(
    angle=360.0,
    flipRevolveDirection=OFF,
    sketch=elliptic_sketch,
    sketchOrientation=RIGHT,
    sketchPlane=gusset_part.faces.findAt((0.0, -1, 1), ),
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=gusset_part.edges.findAt((0.0, -1, 0), )
    )

# Partition gusset

# create assisting datum points to partition the lips
#gusset_part.DatumPointByCoordinate((gp1[0]-l_lip*cos(5*pi/6), gp1[1]-l_lip*sin(5*pi/6), 0),)
#gusset_part.DatumPointByCoordinate((gp1[0]-l_lip*cos(5*pi/6), gp1[1]-l_lip*sin(5*pi/6), diameter),)
#
#gusset_part.DatumPointByCoordinate((gp2[0]-l_lip*cos(-pi/2), gp2[1]-l_lip*sin(-pi/2), 0),)
#gusset_part.DatumPointByCoordinate((gp2[0]-l_lip*cos(-pi/2), gp2[1]-l_lip*sin(-pi/2), diameter),)
#
#gusset_part.DatumPointByCoordinate((gp3[0]-l_lip*cos(pi/6), gp3[1]-l_lip*sin(pi/6), 0),)
#gusset_part.DatumPointByCoordinate((gp3[0]-l_lip*cos(pi/6), gp3[1]-l_lip*sin(pi/6), diameter),)
#
## Partition the 3 gusset plate fins to separate the lip parts
#gusset_part.PartitionFaceByShortestPath(
#    faces=gusset_part.faces.getClosest(coordinates=((gp1[0], gp1[1], 0),))[0][0],
#    point1=gusset_part.datum.items()[0][1],
#    point2=gusset_part.datum.items()[1][1],
#    )
#
#gusset_part.PartitionFaceByShortestPath(
#    faces=gusset_part.faces.getClosest(coordinates=((gp2[0], gp2[1], 0),))[0][0],
#    point1=gusset_part.datum.items()[2][1],
#    point2=gusset_part.datum.items()[3][1],
#    )
#
#gusset_part.PartitionFaceByShortestPath(
#    faces=gusset_part.faces.getClosest(coordinates=((gp3[0], gp3[1], 0),))[0][0],
#    point1=gusset_part.datum.items()[4][1],
#    point2=gusset_part.datum.items()[5][1],
#    )

# Material ----------------------------------------------------------------------------------------------------------------------

stc_mdl.Material(name='elastic')
stc_mdl.materials['elastic'].Elastic(table=((parameters.Youngs_modulus, 0.3), ))

# Create sections ---------------------------------------------------------------------------------------------------------------

# -for sector
stc_mdl.HomogeneousShellSection(
    idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON,
    material='elastic',
    name='sector',
    numIntPts=5,
    poissonDefinition=DEFAULT,
    preIntegrate=OFF,
    temperature=GRADIENT, 
    thickness=cs_thickness,
    thicknessField='',
    thicknessModulus=None,
    thicknessType=UNIFORM,
    useDensity=OFF
    )
    
# -for gusset
stc_mdl.HomogeneousShellSection(
    idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON,
    material='elastic',
    name='gusset',
    numIntPts=5,
    poissonDefinition=DEFAULT,
    preIntegrate=OFF,
    temperature=GRADIENT, 
    thickness=gusset_thickness,
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
#seedsize = pi*d_circumscribed/100

# -Sector
sector_part.setMeshControls(
    algorithm=MEDIAL_AXIS,
    elemShape=QUAD, 
    regions=sector_part.faces[:]
    )
sector_part.seedPart(
    deviationFactor=0.1, 
    minSizeFactor=0.1,
    size=parameters.elem_size
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
    size=parameters.elem_size
    )
gusset_part.generateMesh()

# Create assembly ---------------------------------------------------------------------------------------------------------------

assmbl=stc_mdl.rootAssembly
assmbl.DatumCsysByDefault(CARTESIAN)

# -Sectors
s1_instance=assmbl.Instance(
    dependent=ON,
    name='sector-1',
    part=sector_part
    )
assmbl.DatumAxisByPrincipalAxis(
    principalAxis=ZAXIS
    )
s3_instance=assmbl.RadialInstancePattern(
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
g1_instance=assmbl.Instance(
    dependent=ON,
    name='gusset-1',
    part=gusset_part
    )

g2_instance=assmbl.Instance(
    dependent=ON,
    name='gusset-2',
    part=gusset_part
    )

# --Translate and rotate them to the right position
g2_instance.translate(
    vector=(0.0, 0.0, (l_tot - l_joint))
    )

assmbl.rotate(
    angle=180.0,
    axisDirection=(0.0, 1, 0.0),
    axisPoint=(0.0, 0.0, l_joint / 2),
    instanceList=('gusset-1', )
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

# Position of the holes on the cross-section (x-y coordinates)
sh11 = np.array([x_sector[0] + bolts_w * cos(-pi / 6), y_sector[0] + bolts_w * sin(-pi / 6)])
sh12 = np.array([x_sector[-1] + bolts_w * cos(7 * pi / 6), y_sector[-1] + bolts_w * sin(7 * pi / 6)])

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

# Create reference points for the bolt rigid body couplings

# Create the necessary sets and the tie constraints for all the bolts

# End 1 connection
for oo in (range(3)):
    ii=1
    for o in tuple(bolts_z1):
        
        assmbl.ReferencePoint((gh[oo-3][0], gh[oo-3][1], float(o)))
        
        assmbl.Set(
            edges=s_instance[oo-3].edges.findAt(((sh[oo-3][0][0], sh[oo-3][0][1], float(o)-d_washer/2), ), )+\
            s_instance[oo-2].edges.findAt(((sh[oo-2][1][0], sh[oo-2][1][1], float(o)-d_washer/2), ), )+\
            g1_instance.edges.findAt(((gh[oo-3][0], gh[oo-3][1], float(o)-d_washer/2), ), ),
            name='b'+str(ii)+str(oo)+'set1'
            )
        
        stc_mdl.RigidBody(
            name='b1'+str(ii)+str(oo)+'joint1',
            refPointRegion=Region(referencePoints=(assmbl.referencePoints.findAt((gh[oo-3][0], gh[oo-3][1], float(o))), )),
            tieRegion=assmbl.sets['b'+str(ii)+str(oo)+'set1']
            )
        
        ii+=1
            
# Span

for oo in (range(3)):
    ii=1
    for o in tuple(bolts_z2):
        
        assmbl.ReferencePoint((gh[oo-3][0], gh[oo-3][1], float(o)))
        
        assmbl.Set(
            edges=s_instance[oo-3].edges.findAt(((sh[oo-3][0][0], sh[oo-3][0][1], float(o)-d_washer/2), ), )+\
            s_instance[oo-2].edges.findAt(((sh[oo-2][1][0], sh[oo-2][1][1], float(o)-d_washer/2), ), ),
            name='b'+str(ii)+str(oo)+'-set2'
            )
        
        stc_mdl.RigidBody(
            name='b1'+str(ii)+str(oo)+'span1',
            refPointRegion=Region(referencePoints=(assmbl.referencePoints.findAt((gh[oo-3][0], gh[oo-3][1], float(o))), )),
            tieRegion=assmbl.sets['b'+str(ii)+str(oo)+'-set2']
            )
        
        ii+=1

# End 2 connection

for oo in (range(3)):
    ii=1
    for o in tuple(bolts_z3):
        
        assmbl.ReferencePoint((gh[oo-3][0], gh[oo-3][1], float(o)))
        
        assmbl.Set(
            edges=s_instance[oo-3].edges.findAt(((sh[oo-3][0][0], sh[oo-3][0][1], float(o)-d_washer/2), ), )+\
            s_instance[oo-2].edges.findAt(((sh[oo-2][1][0], sh[oo-2][1][1], float(o)-d_washer/2), ), )+\
            g2_instance.edges.findAt(((gh[oo-3][0], gh[oo-3][1], float(o)-d_washer/2), ), ),
            name='b'+str(ii)+str(oo)+'set3'
            )
        
        stc_mdl.RigidBody(
            name='b1'+str(ii)+str(oo)+'joint2',
            refPointRegion=Region(referencePoints=(assmbl.referencePoints.findAt((gh[oo-3][0], gh[oo-3][1], float(o))), )),
            tieRegion=assmbl.sets['b'+str(ii)+str(oo)+'set3']
            )
        
        ii+=1


# Create reference points for BCs/loads.

# -RPs for the faces at the two ends of the columns
assmbl.ReferencePoint((0.0, 0.0, 0.0))
assmbl.ReferencePoint((0.0, 0.0, l_tot))

rp_end1 = assmbl.referencePoints.findAt((0, 0, 0))
rp_end2 = assmbl.referencePoints.findAt((0, 0, l_tot))

# - End face couplings to reference points

# End 1
rp_end1_set = assmbl.Set(
    name='RP-1-set', 
    referencePoints = (rp_end1, )
    )

if bool(parameters.load_cross_section):
    end1_face_set = assmbl.Set(
        edges=g1_instance.edges.getByBoundingBox(-d_circumscribed,-d_circumscribed,0,d_circumscribed,d_circumscribed,0)+\
        s_instance[0].edges.getByBoundingBox(-d_circumscribed,-d_circumscribed,0,d_circumscribed,d_circumscribed,0)+\
        s_instance[1].edges.getByBoundingBox(-d_circumscribed,-d_circumscribed,0,d_circumscribed,d_circumscribed,0)+\
        s_instance[2].edges.getByBoundingBox(-d_circumscribed,-d_circumscribed,0,d_circumscribed,d_circumscribed,0),
        name='end1-face',
        )
else:
    end1_face_set = assmbl.Set(
        edges=g1_instance.edges.getByBoundingBox(-d_circumscribed,-d_circumscribed,0,d_circumscribed,d_circumscribed,0),
        name='end1-face',
        )

stc_mdl.Coupling(
    controlPoint=rp_end1_set, 
    couplingType=KINEMATIC,
    influenceRadius=WHOLE_SURFACE,
    localCsys=None,
    name='end1-coupling', 
    surface=end1_face_set,
    u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
    )

# End 2

rp_end2_set = assmbl.Set(
    name='RP-2-set',
    referencePoints=(rp_end2, )
    )

if bool(parameters.load_cross_section):
    end2_face_set = assmbl.Set(
        edges=g2_instance.edges.getByBoundingBox(-d_circumscribed,-d_circumscribed, l_tot,d_circumscribed,d_circumscribed,l_tot)+\
        s_instance[0].edges.getByBoundingBox(-d_circumscribed,-d_circumscribed,l_tot,d_circumscribed,d_circumscribed,l_tot)+\
        s_instance[1].edges.getByBoundingBox(-d_circumscribed,-d_circumscribed,l_tot,d_circumscribed,d_circumscribed,l_tot)+\
        s_instance[2].edges.getByBoundingBox(-d_circumscribed,-d_circumscribed,l_tot,d_circumscribed,d_circumscribed,l_tot),
        name='end2-face'
        )
else:
    end2_face_set = assmbl.Set(
        edges=g2_instance.edges.getByBoundingBox(-d_circumscribed, -d_circumscribed, l_tot, d_circumscribed, d_circumscribed, l_tot),
        name='end2-face'
        )

stc_mdl.Coupling(
    controlPoint=rp_end2_set, 
    couplingType=KINEMATIC, influenceRadius=WHOLE_SURFACE,
    localCsys=None,
    name='end2-coupling', 
    surface=end2_face_set,
    u1=ON, u2=ON, u3=ON, ur1=ON, ur2=ON, ur3=ON
    )

# Step -----------------------------------------------------------------------------------
stc_mdl.StaticStep(
    name='IMP',
    previous='Initial')


# Boundary Conditions --------------------------------------------------------------

# BCs                    
end1_BC=stc_mdl.DisplacementBC(
    amplitude=UNSET, 
    createStepName='Initial', 
    distributionType=UNIFORM, 
    fieldName='', 
    localCsys=None, 
    name='fix-end1', 
    region=rp_end1_set, 
    u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=SET
    )

end2_BC=stc_mdl.DisplacementBC(
    amplitude=UNSET,
    createStepName='Initial', 
    distributionType=UNIFORM, 
    fieldName='', 
    localCsys=None, 
    name='fix-end2', 
    region=rp_end2_set, 
    u1=SET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=SET
    )

# Load -------------------------------------------------------------------------------

# Start-end points of the 3 sectors
ss11=np.array([x_sector[0], y_sector[0]])
ss12=np.array([x_sector[-1], y_sector[-1]])

#rotate to get start-end ponts of the other two sectors
ss21=ss11.dot(Rmat)
ss22=ss12.dot(Rmat)
ss31=ss21.dot(Rmat)
ss32=ss22.dot(Rmat)

for i in range(1, len(bolts_z2)):
    assmbl.Surface(
        name='lip_edge_1'+str(i),
        side1Edges=s1_instance.edges.findAt(((ss11[0], ss11[1], bolts_z2[i]-2*d_washer), ), )+\
            s2_instance.edges.findAt(((ss22[0], ss22[1], bolts_z2[i]-2*d_washer), ), )+\
            s3_instance.edges.findAt(((ss31[0], ss31[1], bolts_z2[i]-2*d_washer), ), )+\
            s1_instance.edges.findAt(((ss12[0], ss12[1], bolts_z2[i]-2*d_washer), ), )
        )
    
    assmbl.Surface(
        name='lip_edge_2'+str(i),
        side1Edges=s2_instance.edges.findAt(((ss21[0], ss21[1], bolts_z2[i]-2*d_washer), ), )+\
            s3_instance.edges.findAt(((ss32[0], ss32[1], bolts_z2[i]-2*d_washer), ), )
        )
    
    stc_mdl.ShellEdgeLoad(
        createStepName='IMP',
        distributionType=UNIFORM, 
        field='',
        localCsys=None,
        magnitude=1.0*(-1)**i,
        name='Load-1'+str(i),
        region=assmbl.surfaces['lip_edge_1'+str(i)]
        )
    
    stc_mdl.ShellEdgeLoad(
        createStepName='IMP',
        distributionType=UNIFORM, 
        field='',
        localCsys=None,
        magnitude=1.0*(-1)**(i+1),
        name='Load-2'+str(i),
        region=assmbl.surfaces['lip_edge_2'+str(i)]
        )

# Field output request. Request displacements 'U'

stc_mdl.fieldOutputRequests.changeKey(
    fromName='F-Output-1', 
    toName='fields'
    )

stc_mdl.fieldOutputRequests['fields'].setValues(
    variables=('U',)
    )

## Save the model -------------------------------------------------------------------------------------------------------
mdb.saveAs(pathName=os.getcwd()+'/'+IDstring+'.cae')

# RIKS model ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

riks_model_name = 'RIKS'

# Copy the previous static model
riks_mdl=mdb.Model(
    name=riks_model_name,
    objectToCopy=stc_mdl
    )

# Delete buckling step
del riks_mdl.steps['IMP']

# Create RIKS step
riks_mdl.StaticRiksStep(
    name='RIKS',
    previous='Initial',
    nlgeom=ON,
    maxNumInc=parameters.max_RIKS_increments,
    extrapolation=LINEAR,
    initialArcInc=0.1,
    minArcInc=1e-07,
    totalArcLength=0.5
    )

# Create elastoplastic material
riks_mdl.Material(name='optim355')
riks_mdl.materials['optim355'].Elastic(
    table=((parameters.Youngs_modulus, 0.3), )
    )
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

# Change the section material for the sector to the elastoplastic
riks_mdl.sections['sector'].setValues(
    material='optim355'
    )

# Create a set to act as a control point for the solver killer subroutine
assmbl_riks=riks_mdl.rootAssembly

assmbl_riks.Set(
    name='RIKS_NODE', 
    nodes=assmbl_riks.instances['gusset-2'].nodes.getByBoundingSphere((0, 0, l_tot), parameters.elem_size * 0.1)
    )

# Create variables for the three sector instances
sec1_inst = assmbl_riks.instances['sector-1']
sec2_inst = assmbl_riks.instances['sector-1-rad-2']
sec3_inst = assmbl_riks.instances['sector-1-rad-3']
sectors = (sec1_inst, sec2_inst, sec3_inst)

# Hard contact interaction
if bool(parameters.hard_contact):
    riks_mdl.ContactProperty('contact_prop')
    
    riks_mdl.interactionProperties['contact_prop'].NormalBehavior(
        allowSeparation=ON,
        constraintEnforcementMethod=DEFAULT, 
        pressureOverclosure=HARD
        )
        
    contact_int = riks_mdl.ContactStd(
        createStepName='Initial',
        name='lip_contact'
        )
    
    contact_int.contactPropertyAssignments.appendInStep(
        assignments=((GLOBAL, SELF, 'contact_prop'), ),
        stepName='Initial'
        )
    
    for sec_num in (0, 1, 2):
        assmbl_riks.Surface(
            name='lip0'+str(sec_num)+'1',
            side12Faces=sectors[sec_num].faces.findAt(((sh[sec_num][0][0], sh[sec_num][0][1], bolts_z2[0] - l_lip), ))
            )
            
        assmbl_riks.Surface(
            name='lip0'+str(sec_num)+'2',
            side12Faces=sectors[sec_num].faces.findAt(((sh[sec_num][1][0], sh[sec_num][1][1], bolts_z2[0] - l_lip), ))
            )
    
    for bolt_num in range(len(bolts_z2)):
        for sec_num in (0, 1, 2):
            assmbl_riks.Surface(
                name='lip'+str(bolt_num + 1)+str(sec_num)+'1',
                side12Faces=sectors[sec_num].faces.findAt(((sh[sec_num][0][0], sh[sec_num][0][1], bolts_z2[bolt_num] + l_lip), ))
                )
                
            assmbl_riks.Surface(
                name='lip'+str(bolt_num + 1)+str(sec_num)+'2',
                side12Faces=sectors[sec_num].faces.findAt(((sh[sec_num][1][0], sh[sec_num][1][1], bolts_z2[bolt_num] + l_lip), ))
                )
    
    for bolt_num in range(len(bolts_z2) + 1):
        contact_int.includedPairs.setValuesInStep(
            addPairs=(
                (assmbl_riks.surfaces['lip'+str(bolt_num)+'01'], assmbl_riks.surfaces['lip'+str(bolt_num)+'12']),
                (assmbl_riks.surfaces['lip'+str(bolt_num)+'21'], assmbl_riks.surfaces['lip'+str(bolt_num)+'02']),
                (assmbl_riks.surfaces['lip'+str(bolt_num)+'11'], assmbl_riks.surfaces['lip'+str(bolt_num)+'22']),
                ),
            stepName='Initial', 
            useAllstar=OFF
            )
        
        contact_int.surfaceThicknessAssignments.appendInStep(
            assignments=(
                (assmbl_riks.surfaces['lip' + str(bolt_num)+'01'], ORIGINAL, (1 + parameters.gusset_thickness_ratio)),
                (assmbl_riks.surfaces['lip' + str(bolt_num)+'02'], ORIGINAL, (1 + parameters.gusset_thickness_ratio)),
                (assmbl_riks.surfaces['lip' + str(bolt_num)+'11'], ORIGINAL, (1 + parameters.gusset_thickness_ratio)),
                (assmbl_riks.surfaces['lip' + str(bolt_num)+'12'], ORIGINAL, (1 + parameters.gusset_thickness_ratio)),
                (assmbl_riks.surfaces['lip' + str(bolt_num)+'21'], ORIGINAL, (1 + parameters.gusset_thickness_ratio)),
                (assmbl_riks.surfaces['lip' + str(bolt_num)+'22'], ORIGINAL, (1 + parameters.gusset_thickness_ratio)),
                ), 
            stepName='Initial'
            )

# Overall buckling Imperfections

# Make assmblly instances independent o parts

assmbl_riks.makeIndependent(
    instances=(
        assmbl_riks.instances['sector-1'], 
        assmbl_riks.instances['sector-1-rad-2'], 
        assmbl_riks.instances['sector-1-rad-3'], 
        assmbl_riks.instances['gusset-1'], 
        assmbl_riks.instances['gusset-2']
        )
    )

# Translate nodes for global imperfections

for i in assmbl_riks.allInstances.items():
    for j in range(len(i[1].nodes)):
        xi = i[1].nodes[j].coordinates[0]
        yi = i[1].nodes[j].coordinates[1]
        zi = i[1].nodes[j].coordinates[2]
        assmbl_riks.editNode(
            nodes=i[1].nodes[j],
            offset1=(l_tot/(parameters.bow_imperfections))*sin(pi*zi/l_tot)*cos(parameters.imperfections_angle),
            offset2=(l_tot/(parameters.bow_imperfections))*sin(pi*zi/l_tot)*sin(parameters.imperfections_angle)
            )

# Apply concentrated force

riks_mdl.ConcentratedForce(
    cf3=-N_b_rd,
    createStepName='RIKS',
    distributionType=UNIFORM,
    field='',
    localCsys=None,
    name='compression',
    region=riks_mdl.rootAssembly.sets['RP-2-set']
    )

# Field and History output requests

riks_mdl.historyOutputRequests.changeKey(
    fromName='H-Output-1',
    toName='load'
    )
    
riks_mdl.historyOutputRequests['load'].setValues(
    rebar=EXCLUDE,
    region=riks_mdl.rootAssembly.sets['RP-1-set'], 
    sectionPoints=DEFAULT, variables=('RF3', )
    )
    
riks_mdl.HistoryOutputRequest(
    createStepName='RIKS',
    name='disp',
    rebar=EXCLUDE,
    region=riks_mdl.rootAssembly.sets['RP-2-set'],
    sectionPoints=DEFAULT,
    variables=('U3', )
    )

riks_mdl.fieldOutputRequests.changeKey(
    fromName='F-Output-1', 
    toName='fields'
    )
riks_mdl.fieldOutputRequests['fields'].setValues(
    variables=('S', 'MISES', 'E', 'PEEQ', 'U')
    )

# Edit keywords, create and submit jobs
# The static analysis is first run. Then the maximum displacement is found in the static model
# results and used as a normalization factor for the imperfection amplitude

# Edit the keywords for the static model to write 'U' magnitude displacements
stc_mdl.keywordBlock.synchVersions(storeNodesAndElements=False)
stc_mdl.keywordBlock.insert(xtr.GetBlockPosition(stc_mdl,'*End Step')-1, '*NODE FILE\nU')


# Create and submit the jobs -------------------------------------------------------------------

# Static model
stc_job=mdb.Job(
    atTime=None,
    contactPrint=OFF,
    description='',
    echoPrint=OFF, 
    explicitPrecision=SINGLE,
    getMemoryFromAnalysis=True,
    historyPrint=OFF, 
    memory=90,
    memoryUnits=PERCENTAGE,
    model=static_model_name,
    modelPrint=OFF, 
    multiprocessingMode=DEFAULT,
    name=IDstring+'-imp',
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

# Submit the static job and wait for the results
stc_job.submit() 
stc_job.waitForCompletion()

# find the maximum displacement from the static analysis (open and close the IMP odb)
stc_odb = xtr.open_odb(IDstring+'-imp.odb')
Umax = xtr.max_result(stc_odb, ['U', 'Magnitude'])
odbAccess.closeOdb(stc_odb)

# Calculate the imperfection amplitude based on max displacement
a_factor = bolts_dist/(parameters.distortional_imperfections*Umax)

# Edit the keywords for the compression riks model to include imperfections from the static analysis
riks_mdl.keywordBlock.synchVersions(storeNodesAndElements=False)
riks_mdl.keywordBlock.replace(xtr.GetBlockPosition(riks_mdl, '*step')-1, 
'\n** ----------------------------------------------------------------\n** \n**********GEOMETRICAL IMPERFECTIONS\n*IMPERFECTION,FILE='
+ IDstring+'-imp' +',STEP=1\n1,'+ str(a_factor)+'\n**')

riks_mdl.keywordBlock.insert(xtr.GetBlockPosition(riks_mdl,'*End Step')-1, '\n*NODE FILE, GLOBAL=YES, NSET=RIKS_NODE\nU')

# Riks model
riks_job=mdb.Job(
    atTime=None,
    contactPrint=OFF,
    description='',
    echoPrint=OFF, 
    explicitPrecision=SINGLE,
    getMemoryFromAnalysis=True,
    historyPrint=OFF, 
    memory=90,
    memoryUnits=PERCENTAGE,
    model=riks_model_name,
    modelPrint=OFF, 
    multiprocessingMode=DEFAULT,
    name=IDstring+'-riks',
    nodalOutputPrecision=SINGLE,
    numCpus=1,
    numGPUs=0,
    queue=None,
    resultsFormat=ODB,
    scratch='', 
    type=ANALYSIS,
    userSubroutine=os.getcwd()+'/GN_Riks_killer.f',
    waitHours=0,
    waitMinutes=0
    )

## Save the model -------------------------------------------------------------------------------------------------------
mdb.saveAs(pathName=os.getcwd()+'/'+IDstring+'.cae')

# Submit the riks job and wait for the results
riks_job.submit()

# Delete initial model
if mdb.models.has_key('Model-1'):
    del mdb.models['Model-1']

# Save the model -------------------------------------------------------------------------------------------------------
mdb.saveAs(pathName=os.getcwd()+'/'+IDstring+'.cae')

# Return to parent directory
os.chdir('..')
