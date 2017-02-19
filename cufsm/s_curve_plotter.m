function [ax] = s_curve_plotter(curves)
% Signature curve plotter for CUFSM resilts

% Total number lengths
n_lengths = length(curves);

% Total number of eigenvalues
neigs = length(curves{1});

% Initialise a cell arrays to host x, y data
curvex = cell(n_lengths, neigs);
curvey = cell(n_lengths, neigs);

% Open a figure window
figure;

% Loop collecting the y values (load factor) for the different lengths and
% plotting them for all the eigenvalues in the same figure
for i_eig = 1:neigs;
    for j_hwave = 1:n_lengths;
        curvex{i_eig}(j_hwave) = curves{j_hwave}(1, 1);
        curvey{i_eig}(j_hwave) = curves{j_hwave}(i_eig, 2);
    end
    semilogx(curvex{i_eig}, curvey{i_eig});
hold on;
end;

% Plot limits
xmin=curvex{1}(1)*10/11;
xmax=curvex{1}(n_lengths)*11/10;
ymin=0;
ymax=3*median(curvey{1});
axis([xmin, xmax, ymin, ymax]);