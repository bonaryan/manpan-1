function [ax] = s_curve_plotter(curves, i, j, k)
% Signature curve plotter for CUFSM resilts

n_lengths = length(curves{1, 1, 1});
neigs = length(curves{1, 1, 1}{1});
curvex = cell(n_lengths, neigs);
figure;
for i_eig = 1:neigs;
    for j_hwave = 1:n_lengths;
        curvex{i_eig}(j_hwave) = curves{i, j, k}{j_hwave}(1, 1);
        curvey{i_eig}(j_hwave) = curves{i, j, k}{j_hwave}(i_eig, 2);
    end
    semilogx(curvex{i_eig}, curvey{i_eig});
hold on;
end
xmin=curvex{1}(1)*10/11;
ymin=0;
xmax=curvex{1}(n_lengths)*11/10;
ymax=3*median(curvey{1});
axis([xmin, xmax, ymin, ymax]);