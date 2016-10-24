% Script that prints a .jnl ASCII file for Abaqus. The .jnl creates in
% Abaqus all the models from a given set of xy points.


% attempt to write to multiple lines to file
jnltext = '# Save by user on 2016_10_24-12.23.43; build 6.14-1 2014_06_05-00.11.02 134264\r\n% from part import *'
% from material import *
% from section import *
% from assembly import *
% from step import *
% from interaction import *
% from load import *
% from mesh import *
% from optimization import *
% from job import *
% from sketch import *
% from visualization import *
% from connectorBehavior import *
% mdb.models.changeKey(fromName='Model-1', toName='p111')
% mdb.models['p111'].ConstrainedSketch(name='__profile__', sheetSize=1200.0)
% mdb.models['p111'].sketches['__profile__'].Line(point1=(-245.0, 45.0), point2=(
%     -195.0, 0.0))
% mdb.models['p111'].sketches['__profile__'].Line(point1=(-195.0, 0.0), point2=(
%     0.0, 110.0))
% mdb.models['p111'].sketches['__profile__'].Line(point1=(0.0, 110.0), point2=(
%     205.0, 0.0))
% mdb.models['p111'].sketches['__profile__'].Line(point1=(205.0, 0.0), point2=(
%     260.0, 40.0))
% mdb.models['p111'].Part(dimensionality=THREE_D, name='sector', type=
%     DEFORMABLE_BODY)
% mdb.models['p111'].parts['sector'].BaseShellExtrude(depth=1200.0, sketch=
%     mdb.models['p111'].sketches['__profile__'])
% del mdb.models['p111'].sketches['__profile__']
% mdb.models['p111'].Material(name='pure-elastic')
% mdb.models['p111'].materials['pure-elastic'].Elastic(table=((210000.0, 0.3), ))
% # Save by user on 2016_10_24-12.25.31; build 6.14-1 2014_06_05-00.11.02 134264
% '

fid = fopen('polygons.jnl', 'w');

fwrite(fid, jnltext);

fclose(fid);