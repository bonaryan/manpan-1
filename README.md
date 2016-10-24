# Polygoner

Script creating sectors of polygonal profiles.

# pcoords.m

Matlab function calculating the x, y coordinates of one sector of a semi-closed polygonal cross section. The parent cross section consists of 3 identical sectors (120 degrees each).

[x, y] = pcoords(n, d, fy, rcoef, nbend, lext, tg, slend)

Input:

n		(integer)	Number of corners (entire polygon, only values 3*m)
d		(double)	Polygon diameter
slend	(double)	Slenderness
fy		(double)	Yield strength
rcoef	(double)	Bending radius to thickness ratio (r/t = rcoef)
nbend	(integer)	Number of points along the bend
lext	(double)	extension length
tg		(double)	Thickness of the gusset plate

Output:

x		(array)	x-coordinates
y		(array)	y coordinates

# polygoner.m

Matlab function that executes pcoords and outputs a matrix of arrays of x, y coordinates of polygonal sectors for a given range of input values.

profile_matrix = polygoner([n], [d], [slend], fy, rcoef, nbend, lext, tg)

Input:

n		(array)		Number of corners (e.g. [6, 9, 12])
d		(array)		Polygon diameters. Initial:step:final (e.g. [300:50:500])
slend	(array)	Slenderness range. Initial:step:final (e.g. [80:1:500])
fy		(double)	Yield strength
rcoef	(double)	Bending radius to thickness ratio (r/t = rcoef)
nbend	(integer)	Number of points along the bend
lext	(double)	extension length
tg		(double)	Thickness of the gusset plate

Output:

profile_matrix	(3d cell array)		Cell array containing all the x-y arrays.

# Polygoner_ABQS_out.m

Scripts that exports an abaqus journal file with shell models from the given x-y arrays.

# Polygoner_CFSM_out.m

Scripts that executes CUFSM for given x-y arrays.

