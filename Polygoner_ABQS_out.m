% Script that prints a .jnl ASCII file for Abaqus. The .jnl creates in
% Abaqus all the models from a given set of xy points.


% attempt to write to multiple lines to file (not finished yet, need to include variable in the text and combine all models)

jnltext_head = {'from part import *'
'from material import *'
'from section import *'
'from assembly import *'
'from step import *'
'from interaction import *'
'from load import *'
'from mesh import *'
'from optimization import *'
'from job import *'
'from sketch import *'
'from visualization import *'
'from connectorBehavior import *'}

jnltext_data = {['mdb.models.changeKey(fromName=''Model-1'',toName=''p', num2str(111), ''')']
    ['mdb.models[''p111''].ConstrainedSketch(name=''__profile__'', sheetSize=1200.0)']}
for i = 1 : length(profiles{1, 1, 1}(1, :))-1;
jnltext_data(end+1, 1) = {['mdb.models[''', num2str(111), '''].sketches[''__profile__''].Line(point1=(', num2str(profiles{1, 1, 1}(1, i)), ', ', num2str(profiles{1, 1, 1}(2, i)), '), point2=(', num2str(profiles{1, 1, 1}(1, i+1)), ', ', num2str(profiles{1, 1, 1}(2, i+1)), '))']};
end;

jnltext_data(end+1, 1) = {'mdb.models[''p111''].Part(dimensionality=THREE_D, name=''sector'', type=DEFORMABLE_BODY)'};
jnltext_data(end+1, 1) = {'mdb.models[''p111''].parts[''sector''].BaseShellExtrude(depth=1200.0, sketch=mdb.models[''p111''].sketches[''__profile__''])'};
jnltext_data(end+1, 1) = {'del mdb.models[''p111''].sketches[''__profile__'']'};
jnltext_data(end+1, 1) = {'mdb.models[''p111''].Material(name=''pure-elastic'')'};
jnltext_data(end+1, 1) = {'mdb.models[''p111''].materials[''pure-elastic''].Elastic(table=((210000.0, 0.3), ))'};

jnltext = [jnltext_head(:, 1); jnltext_data(:, 1)];

fid = fopen('polygons.jnl', 'w');

for i = 1:length(jnltext);
    fwrite(fid, jnltext{i});
    fwrite(fid, char(10));
end;

fclose(fid);

