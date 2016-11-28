function [x_out, y_out, t, tg, l_lip] = pcoords(n, d, slend, fy, rcoef, nbend, l_ratio, t_ratio)
% Return x, y coords of points of a 1/3 of folded polygonal cross section.
% input args: number of corners, CS diameter, slenderness, yield strength, 
% bending arc radius r/t, no. of points along the bending arcs,lip length 
% to diameter ratio, gusset plate thickness to sector thickness ratio.
% output: [x, y], sector thickness, gussetplate thickness

%% Input (recomended values)
% % Number of corners (entire polygon, only 3*m)
% n = 9;
% 
% % Polygon diameter
% d = 500;
% 
% % Yield strength
% fy = 355;
% 
% % Bending radius to thickness ratio
% % (r/t = rcoef)
% rcoef = 6;
% 
% % Number of points along the bend
% nbend = 6;
% 
% % extension length to diameter ratio
% l_ratio = 0.1;
% 
% % Thickness of the gusset plate to sector thickness ratio
% t_ratio = 1.20;
%
% % Slenderness
% slend = 90;


% Calculated characteristics
R = d/2;
epsilon = sqrt(fy/235);
t = round(epsilon^2 * d / slend);
tg = round(t_ratio*t);
l_lip = l_ratio*d;

%% Polygon sector
% Angle corresponding to one edge of the polygon
theta = 2*pi/n;

% Angles of radii (measured from x-axis)
phi=5*pi/6:-theta:pi/6;

% xy coords of the polygon's corners
x = R*cos(phi);
y = R*sin(phi);

%% Bends
% Bending radius
rbend = rcoef*t;

% Distance between bending centre and corner
lc = rbend/cos(theta/2);

% Centers of bending arcs
xc  = (x(2:end-1) - lc*cos(phi(2:end-1)));
yc  = (y(2:end-1) - lc*sin(phi(2:end-1)));

% Bending arc angle
theta_b = pi - theta;

% Angles of the edges' midlines (measured from x-axis)
phi_mids = phi(1:end-1) - theta/2 ;

% xy coords of the arc's points
for i = 1:n/3 -1;
    for j = 1:nbend+1;
        xarc(i, j) = xc(i) + rbend*cos(phi_mids(i)-(j-1)*(theta/nbend));
        yarc(i, j) = yc(i) + rbend*sin(phi_mids(i)-(j-1)*(theta/nbend));
    end;
end;

%% Start-end extensions
% Bending radius
rs = rbend/5;

% First bend
v1 = phi_mids(1)-pi/2;
v2 = (phi(1)+phi_mids(1)-pi/2)/2;
l1 = (t+tg)/(2*cos(phi(1)-phi_mids(1)));
l2 = rs/sin(v2-phi_mids(1)+pi/2);
x1 = x(1)+l1*cos(v1);
y1 = y(1)+l1*sin(v1);

% First bend centre coords
xcs(1) = x1+l2*cos(v2);
ycs(1) = y1+l2*sin(v2);

% Last bend
v1 = phi_mids(end)+pi/2;
v2 = (v1+phi(end))/2;
l1 = (t+tg)/(2*cos(v1-phi(end)-pi/2));
l2 = rs/sin(v2-phi(end));
x1 = x(end)+l1*cos(v1);
y1 = y(end)+l1*sin(v1);

% Last bend centre coords
xcs(2) = x1+l2*cos(v2);
ycs(2) = y1+l2*sin(v2);

% First and last bend arc points coords
for j = 1:nbend+1;
    xsarc(1, j) = xcs(1) + rs*cos(4*pi/3+(j-1)*((phi_mids(1)-pi/3)/nbend));
    ysarc(1, j) = ycs(1) + rs*sin(4*pi/3+(j-1)*((phi_mids(1)-pi/3)/nbend));
    xsarc(2, j) = xcs(2) + rs*cos(phi_mids(end)+pi+(j-1)*((phi(end)+pi/2-phi_mids(end))/nbend));
    ysarc(2, j) = ycs(2) + rs*sin(phi_mids(end)+pi+(j-1)*((phi(end)+pi/2-phi_mids(end))/nbend));
end;

%% Points of the lips
% First lip
xstart = [xsarc(1, 1) + l_lip*cos(phi(1)), xsarc(1, 1) + l_lip*cos(phi(1))/2];
ystart = [ysarc(1, 1) + l_lip*sin(phi(1)), ysarc(1, 1) + l_lip*sin(phi(1))/2];


% Last point
xend = [xsarc(2, end) + l_lip*cos(phi(end))/2, xsarc(2, end) + l_lip*cos(phi(end))];
yend = [ysarc(2, end) + l_lip*sin(phi(end))/2, ysarc(2, end) + l_lip*sin(phi(end))];

%% Collect the x, y values in a sorted 2xn array
xarc = xarc';
yarc = yarc';

x_out = [xstart, xsarc(1, :), xarc(:)', xsarc(2, :), xend];
y_out = [ystart, ysarc(1, :), yarc(:)', ysarc(2, :), yend];

% Plot result
% plot(x_out, y_out);
