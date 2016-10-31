function [profiles] = Polygoner_ABQS_out(n, d, slend, fy, rcoef, nbend, lext, tg)
% Script that prints a .jnl ASCII file for Abaqus. The .jnl creates in
% Abaqus all the models from a given set of xy points.


% attempt to write to multiple lines to file (not finished yet, need to include variable in the text and combine all models) 
n = 9;
d = 300;
slend = 90;
fy = 355;
rcoef = 6;
nbend = 5;
lext = 30;
tg = 10;

% Initialise the cell array to host the profiles
profiles = cell(length(n), length(d), length(slend));

% Loop through the values 
for i = 1:length(n); 
    for j = 1:length(d);
        for k = 1:length(slend);
            [x, y] = pcoords(n(i), d(j), slend(k), fy, rcoef, nbend, lext, tg);
            profiles{i, j, k} = [x; y];
        end
    end
end

%aoutput for Abaqus

        
    
p = max(1:1:length(x))


r = [1:1:p];

aa = 1;
ab = 2;
ac = 3;
ad = 4;
ae = 5;
af = 6;
ag = 7;
ah = 8;
ai = 9;
aj = 10;
ak = 11;
al = 12;
am = 13;
an = 14;
ao = 15;
ap = 16;
aq = 17;
ar = 18;
as = 19;
at = 20;
au = 21;
av = 22; 
aw = 23;
ax = 24;
ay = 25;
az = 26;
ba = 27;
bb = 28;
bc = 29;
bd = 30;
be = 31;
bf = 32;

if aa < p
    aa = 1;
else
    aa = p;
end

if ab < p
    ab = 2;
else
    ab = p;
end

if ac < p
    ac = 3;
else
    ac = p;
end

if ad < p
    ad = 4;
else
    ad = p;
end

if ae < p
    ae = 5;
else
    ae = p;
end

if af < p
    af = 6;
else
    af = p;
end

if ag < p
    ag = 7;
else
    ag = p;
end

if ah < p
    ah = 8;
else
    ah = p;
end

if ai < p
    ai = 9;
else
    ai = p;
end

if aj < p
    aj = 10;
else
    aj = p;
end

if ak < p
    ak = 11;
else
    ak = p;
end

if al < p
    al = 12;
else
    al = p;
end

if am < p
    am = 13;
else
    am = p;
end

if an < p
    an = 14;
else
    an = p;
end

if ao < p
    ao = 15;
else
    ao = p;
end

if ap < p
    ap = 16;
else
    ap = p;
end

if aq < p
    aq = 17;
else
    aq = p;
end

if ar < p
    ar = 18;
else
    ar = p;
end

if as < p
    as = 19;
else
    as = p;
end

if at < p
    at = 20;
else
    at = p;
end

if au < p
    au = 21;
else
    au = p;
end

if av < p
    av = 22;
else
    av = p;
end

if aw < p
    aw = 23;
else
    aw = p;
end

if ax < p
    ax = 24;
else
    ax = p;
end

if ay < p
    ay = 25;
else
    ay = p;
end

if az < p
    az = 26;
else
    az = p;
end

if ba < p
    ba = 27;
else
    ba = p;
end

if bb < p
    bb = 28;
else
    bb = p;
end

if bc < p
    bc = 29;
else
    bc = p;
end

if bd < p
    bd = 30;
else
    bd = p;
end

if be < p
    be = 31;
else
    be = p;
end

if bf < p
    bf = 32;
else
    bf = p;
end


pointsline{r(aa)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(aa)) ', ' num2str(y(aa)) '), point2=(' num2str(x(ab)) ', ' num2str(y(ab)) '))'];
pointsline{r(ab)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ab)) ', ' num2str(y(ab)) '), point2=(' num2str(x(ac)) ', ' num2str(y(ac)) '))'];
pointsline{r(ac)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ac)) ', ' num2str(y(ac)) '), point2=(' num2str(x(ad)) ', ' num2str(y(ad)) '))'];
pointsline{r(ad)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ad)) ', ' num2str(y(ad)) '), point2=(' num2str(x(ae)) ', ' num2str(y(ae)) '))'];
pointsline{r(ae)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ae)) ', ' num2str(y(ae)) '), point2=(' num2str(x(af)) ', ' num2str(y(af)) '))'];
pointsline{r(af)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(af)) ', ' num2str(y(af)) '), point2=(' num2str(x(ag)) ', ' num2str(y(ag)) '))'];
pointsline{r(ag)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ag)) ', ' num2str(y(ag)) '), point2=(' num2str(x(ah)) ', ' num2str(y(ah)) '))'];
pointsline{r(ah)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ah)) ', ' num2str(y(ah)) '), point2=(' num2str(x(ai)) ', ' num2str(y(ai)) '))'];
pointsline{r(ai)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ai)) ', ' num2str(y(ai)) '), point2=(' num2str(x(aj)) ', ' num2str(y(aj)) '))'];
pointsline{r(aj)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(aj)) ', ' num2str(y(aj)) '), point2=(' num2str(x(ak)) ', ' num2str(y(ak)) '))'];
pointsline{r(ak)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ak)) ', ' num2str(y(ak)) '), point2=(' num2str(x(al)) ', ' num2str(y(al)) '))'];
pointsline{r(al)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(al)) ', ' num2str(y(al)) '), point2=(' num2str(x(am)) ', ' num2str(y(am)) '))'];
pointsline{r(am)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(am)) ', ' num2str(y(am)) '), point2=(' num2str(x(an)) ', ' num2str(y(an)) '))'];
pointsline{r(an)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(an)) ', ' num2str(y(an)) '), point2=(' num2str(x(ao)) ', ' num2str(y(ao)) '))'];
pointsline{r(ao)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ao)) ', ' num2str(y(ao)) '), point2=(' num2str(x(ap)) ', ' num2str(y(ap)) '))'];
pointsline{r(ap)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ap)) ', ' num2str(y(ap)) '), point2=(' num2str(x(aq)) ', ' num2str(y(aq)) '))'];
pointsline{r(aq)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(aq)) ', ' num2str(y(aq)) '), point2=(' num2str(x(ar)) ', ' num2str(y(ar)) '))'];
pointsline{r(ar)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ar)) ', ' num2str(y(ar)) '), point2=(' num2str(x(as)) ', ' num2str(y(as)) '))'];
pointsline{r(as)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(as)) ', ' num2str(y(as)) '), point2=(' num2str(x(at)) ', ' num2str(y(at)) '))'];
pointsline{r(at)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(at)) ', ' num2str(y(at)) '), point2=(' num2str(x(au)) ', ' num2str(y(au)) '))'];
pointsline{r(au)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(au)) ', ' num2str(y(au)) '), point2=(' num2str(x(av)) ', ' num2str(y(av)) '))'];
pointsline{r(av)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(av)) ', ' num2str(y(av)) '), point2=(' num2str(x(aw)) ', ' num2str(y(aw)) '))'];
pointsline{r(aw)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(aw)) ', ' num2str(y(aw)) '), point2=(' num2str(x(ax)) ', ' num2str(y(ax)) '))'];
pointsline{r(ax)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ax)) ', ' num2str(y(ax)) '), point2=(' num2str(x(ay)) ', ' num2str(y(ay)) '))'];
pointsline{r(ay)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ay)) ', ' num2str(y(ay)) '), point2=(' num2str(x(az)) ', ' num2str(y(az)) '))'];
pointsline{r(az)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(az)) ', ' num2str(y(az)) '), point2=(' num2str(x(ba)) ', ' num2str(y(ba)) '))'];
pointsline{r(ba)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(ba)) ', ' num2str(y(ba)) '), point2=(' num2str(x(bb)) ', ' num2str(y(bb)) '))'];
pointsline{r(bb)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(bb)) ', ' num2str(y(bb)) '), point2=(' num2str(x(bc)) ', ' num2str(y(bc)) '))'];
pointsline{r(bc)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(bc)) ', ' num2str(y(bc)) '), point2=(' num2str(x(bd)) ', ' num2str(y(bd)) '))'];
pointsline{r(bd)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(bd)) ', ' num2str(y(bd)) '), point2=(' num2str(x(be)) ', ' num2str(y(be)) '))'];
pointsline{r(be)} = ['mdb.models[''p111''].sketches[''__profile__''].Line(point1=(' num2str(x(be)) ', ' num2str(y(be)) '), point2=(' num2str(x(bf)) ', ' num2str(y(bf)) '))'];


     
    jnltext = ['from part import *' char(10) 'from material import *' char(10) 'from section import *' char(10) 'from assembly import *' char(10) 'from step import *' char(10) 'from interaction import *' char(10) 'from load import *' char(10) 'from mesh import *' char(10) 'from optimization import *' char(10) 'from job import *' char(10) 'from sketch import *' char(10) 'from visualization import *' char(10) 'from connectorBehavior import *' char(10) 'mdb.models.changeKey(fromName=''Model-1'',toName=''p111'')' char(10) 'mdb.models[''p111''].ConstrainedSketch(name=''__profile__'', sheetSize=1200.0)' char(10) sprintf(pointsline{aa}) char(10) sprintf(pointsline{ab}) char(10) sprintf(pointsline{ac}) char(10) sprintf(pointsline{ad}) char(10) sprintf(pointsline{ae}) char(10) sprintf(pointsline{af}) char(10) sprintf(pointsline{ag}) char(10) sprintf(pointsline{ah}) char(10) sprintf(pointsline{ai}) char(10) sprintf(pointsline{aj}) char(10) sprintf(pointsline{ak}) char(10) sprintf(pointsline{al}) char(10) sprintf(pointsline{am}) char(10) sprintf(pointsline{an}) char(10) sprintf(pointsline{ao}) char(10) sprintf(pointsline{ap}) char(10) sprintf(pointsline{aq}) char(10) sprintf(pointsline{ar}) char(10) sprintf(pointsline{as}) char(10) sprintf(pointsline{at}) char(10) sprintf(pointsline{au}) char(10) sprintf(pointsline{av}) char(10) sprintf(pointsline{aw}) char(10) sprintf(pointsline{ax}) char(10) sprintf(pointsline{ay}) char(10) sprintf(pointsline{az}) char(10) sprintf(pointsline{ba}) char(10) sprintf(pointsline{bb}) char(10) sprintf(pointsline{bc}) char(10) sprintf(pointsline{bd}) char(10) sprintf(pointsline{be}) char(10) 'mdb.models[''p111''].Part(dimensionality=THREE_D, name=''sector'', type=DEFORMABLE_BODY)' char(10) 'mdb.models[''p111''].parts[''sector''].BaseShellExtrude(depth=1200.0, sketch=mdb.models[''p111''].sketches[''__profile__''])' char(10) 'del mdb.models[''p111''].sketches[''__profile__'']' char(10) 'mdb.models[''p111''].Material(name=''pure-elastic'')' char(10) 'mdb.models[''p111''].materials[''pure-elastic''].Elastic(table=((210000.0, 0.3), ))'];

fid = fopen('polygons.jnl', 'w');

fwrite(fid, jnltext);

fclose(fid);



   





