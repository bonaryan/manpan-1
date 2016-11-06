function [profiles] = polygoner(nrange, drange, slendrange, fy, rcoef, nbend, lext, tg)
% Return a cell array with the points of all the profiles within a range of
% values.
% input args: numbers of corners, CS diameters, slenderness', yield strength, 
% bending arc radius r/t, no. of points along the bending arcs, end 
% extensions length, gusset plate thickness.
% output: [x, y]


% Example input
% nrange = [6, 9, 12];
% drange = [300:200:500];
% slendrange = [70:5:120];
% 
% fy = 355;
% rcoef = 6;
% nbend = 4;
% lext = 20;
% tg = 10;

% Initialise the cell array to host the profiles
profiles = cell(length(nrange), length(drange), length(slendrange));

% Loop through the values within the given ranges
for i = 1:length(nrange); 
    for j = 1:length(drange);
        for k = 1:length(slendrange);
            [x, y] = pcoords(nrange(i), drange(j), slendrange(k), fy, rcoef, nbend, lext, tg);
            profiles{i, j, k} = [x; y];
        end
    end
end

% Save the profile database to the current directory as "profiles.mat"
% save('profiles.mat', 'profiles');