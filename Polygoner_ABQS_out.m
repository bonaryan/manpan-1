% Script that prints a .jnl ASCII file for Abaqus. The .jnl creates in
% Abaqus all the models from a given set of xy points.


% attempt to write to multiple lines to file (not finished yet, need to include variable in the text and combine all models) 
jnltext = ['from part import *' char(10) 'from material import *' char(10) 'from section import *' char(10) 'from assembly import *' char(10) 'from step import *' char(10) 'from interaction import *' char(10) 'from load import *' char(10) 'from mesh import *' char(10) 'from optimization import *' char(10) 'from job import *' char(10) 'from sketch import *' char(10) 'from visualization import *' char(10) 'from connectorBehavior import *' char(10) 'mdb.models.changeKey(fromName=''Model-1'',toName=''p111'')' char(10) 'mdb.models[''p111''].ConstrainedSketch(name=''__profile__'', sheetSize=1200.0)' char(10) 'mdb.models[''p111''].sketches[''__profile__''].Line(point1=(-245.0, 45.0), point2=(-195.0, 0.0))' char(10) 'mdb.models[''p111''].sketches[''__profile__''].Line(point1=(-195.0, 0.0), point2=(0.0, 110.0))' char(10) 'mdb.models[''p111''].sketches[''__profile__''].Line(point1=(0.0, 110.0), point2=(205.0, 0.0))' char(10) 'mdb.models[''p111''].sketches[''__profile__''].Line(point1=(205.0, 0.0), point2=(260.0, 40.0))' char(10) 'mdb.models[''p111''].Part(dimensionality=THREE_D, name=''sector'', type=DEFORMABLE_BODY)' char(10) 'mdb.models[''p111''].parts[''sector''].BaseShellExtrude(depth=1200.0, sketch=mdb.models[''p111''].sketches[''__profile__''])' char(10) 'del mdb.models[''p111''].sketches[''__profile__'']' char(10) 'mdb.models[''p111''].Material(name=''pure-elastic'')' char(10) 'mdb.models[''p111''].materials[''pure-elastic''].Elastic(table=((210000.0, 0.3), ))']

fid = fopen('polygons.jnl', 'w');

fwrite(fid, jnltext);

fclose(fid);

