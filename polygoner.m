function [profiles, meta] = polygoner(nrange, drange, slendrange, fy, rcoef, nbend, l_ratio, t_ratio)
% Return a cell array with the points of all the profiles within a range of
% values.
% input args: numbers of corners, CS diameters, slenderness', yield strength, 
% bending arc radius r/t, no. of points along the bending arcs, end 
% extensions length, gusset plate thickness.
% output: [x; y], [diameter; plate thicness; gusset plate thickness; fy]

% Example input
% nrange = [6, 9, 12];
% drange = [300:200:900];
% slendrange = linspace(70, 150, 10);
% 
% fy = 355;
% rcoef = 6;
% nbend = 4;
% l_ratio = 0.1;
% t_ratio = 1.2;

% Initialise a cell array to host the profiles' xy values
profiles = cell(length(nrange), length(drange), length(slendrange));

% Initialise a cell array to host the profile metadata
meta = cell(length(nrange), length(drange), length(slendrange));

% Loop through the values within the given ranges
for i = 1:length(nrange); 
    for j = 1:length(drange);
        for k = 1:length(slendrange);
            [x, y, t, tg] = pcoords(nrange(i), drange(j), slendrange(k), fy, rcoef, nbend, l_ratio, t_ratio);
            profiles{i, j, k} = [x; y];
            meta{i, j, k} = [drange(j); t; tg; fy];
        end
    end
end

% Save the profile database and metadata to the current directory as .mat
save('profiles.mat', 'profiles');
save('meta.mat', 'meta');