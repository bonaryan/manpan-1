# Polygoner

Script creating sectors of polygonal profiles.

# pcoords.m

Matlab function calculating the x, y coordinates of one sector of a semi-closed polygonal cross section. The parent cross section consists of 3 identical sectors (120 degrees each).

[x, y] = pcoords(n, d, fy, rcoef, nbend, lext, tg, slend)

Input:

n\t\t(integer)\tNumber of corners (entire polygon, only values 3*m)

d\t\t(double)\tPolygon diameter

slend\t(double)\tSlenderness

fy\t\t(double)\tYield strength

rcoef\t(double)\tBending radius to thickness ratio (r/t = rcoef)

nbend\t(integer)\tNumber of points along the bend

lext\t(double)\textension length

tg\t\t(double)\tThickness of the gusset plate


Output:

x\t\t(array)\tx-coordinates
y\t\t(array)\ty coordinates

# polygoner.m

Matlab function that executes pcoords and outputs a matrix of arrays of x, y coordinates of polygonal sectors for a given range of input values.

profile_matrix = polygoner([n], [d], [slend], fy, rcoef, nbend, lext, tg)

Input:

n\t\t(array)\t\tNumber of corners (e.g. [6, 9, 12])
d\t\t(array)\t\tPolygon diameters. Initial:step:final (e.g. [300:50:500])
slend\t(array)\tSlenderness range. Initial:step:final (e.g. [80:1:500])
fy\t\t(double)\tYield strength
rcoef\t(double)\tBending radius to thickness ratio (r/t = rcoef)
nbend\t(integer)\tNumber of points along the bend
lext\t(double)\textension length
tg\t\t(double)\tThickness of the gusset plate

Output:

profile_matrix\t(3d cell array)\t\tCell array containing all the x-y arrays.

# Polygoner_ABQS_out.m

Scripts that exports an abaqus journal file with shell models from the given x-y arrays.

# Polygoner_CFSM_out.m

Scripts that executes CUFSM for given x-y arrays.

