function [profiles, meta] = polygoner(nrange, drange, slendrange, fy, rcoef, nbend, l_ratio, t_ratio, lambda)
% Return a cell array with the points of all the profiles within a range of
% values.
% input args: numbers of corners, CS diameters, slenderness', yield strength,
% bending arc radius r/t, no. of points along the bending arcs, end
% extensions length, gusset plate thickness.
% profiles output  : [x; y], [diameter; plate thicness; gusset plate thickness; fy]
% meta outpur      : d; t; tg; fy; A; Ixx; Izz; Ixz

% Example input
% nrange = [6, 9, 12];
% drange = [300:200:900];
% slendrange = linspace(70, 150, 10);
% lambda = [0.65, 1, 1.25];
%
% fy = 355;
% rcoef = 6;
% nbend = 4;
% l_ratio = 0.1;
% t_ratio = 1.2;

E = 210000;

% Initialise a cell array to host the profiles' xy values
profiles = cell(length(nrange), length(drange), length(slendrange));

% Initialise a cell array to host the profile metadata
meta = cell(length(nrange), length(drange), length(slendrange), length(lambda));

% Loop through the values within the given ranges
for i = 1:length(nrange);
    for j = 1:length(drange);
        for k = 1:length(slendrange);
            
            % Call pcoords to get data for a profile
            [x, y, t, tg] = pcoords(nrange(i), drange(j), slendrange(k), fy, rcoef, nbend, l_ratio, t_ratio);
            
            % Collect the xy values in a database
            profiles{i, j, k} = [x; y];
            
            % Metadata of the profiles
            % Crate node and elem arrays for the profile appropriate for
            % input to cutwp_prop2 function which returns cs properties
            
            % Current profile xy
            c_prof1 = profiles{i, j, k}';
            
            % Number of vertices on the current profile
            l_prof = length(c_prof1);
            
            % Construct the 2 extra parts by rotating the imported one
            R2 = [cos(-2*pi/3), -sin(-2*pi/3); sin(-2*pi/3), cos(-2*pi/3)];
            R3 = [cos(2*pi/3), -sin(2*pi/3); sin(2*pi/3), cos(2*pi/3)];
            for a = 1:l_prof;
                c_prof2(a, :) = (R2*c_prof1(a, :)')';
                c_prof3(a, :) = (R3*c_prof1(a, :)')';
            end;
            
            % A column of ones
            col1 = ones(l_prof', 1);
            
            % Construct the 'node' array
            node = [c_prof1(:, 1), c_prof1(:, 2);
                c_prof2(:, 1), c_prof2(:, 2);
                c_prof3(:, 1), c_prof3(:, 2)];
            
            % Construct the 'elem' array
            elem = [(1:l_prof-1)', (2:l_prof)', t*ones(l_prof-1', 1);
                l_prof, l_prof+1, 0.1;
                l_prof+(1:l_prof-1)', l_prof+(2:l_prof)', t*ones(l_prof-1', 1);
                2*l_prof, 2*l_prof+1, 0.1;
                2*l_prof+(1:l_prof-1)', 2*l_prof+(2:l_prof)', t*ones(l_prof-1', 1);
                3*l_prof, 1, 0.1];
            
            % Return cs properties using cutwp
            [A, ~, ~, Iyy, Izz] = cutwp_prop2(node, elem);
            
            % Current profile area and moment of inertia
            I = min(Iyy, Izz);
            
            % Loop through the different member slendernesses. The 'meta'
            % array has one more dimension (4D)
            for l = 1:length(lambda);
                
                % Current profile length
                len = lambda*pi*sqrt(E*I/(A*fy));
                
                % Store the metadata in a cell array
                meta{i, j, k, l} = [drange(j); t; tg; fy; A; Iyy; Izz; len(l)];
            end
        end
    end
end

% Save the profile database and metadata to the current directory as .mat
save('profiles.mat', 'profiles');
save('meta.mat', 'meta');