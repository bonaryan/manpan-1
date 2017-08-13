% Assignment 2 for the FEM in structural engineering course
% Full factorial design, 3 factors

parameters = {
%    'Number of corners'
    'Diameter'
    'Profile slenderness'
    'Member slenderness'
    'Bolt spacing'
    }';
% Create generators, design space and confounding
dff = fullfact([3 2 3 3]);

% Create the generators with the real values of the factors
dff2 = zeros(54, 4);
dff2(dff(:, 1)==1, 1) = 500;
dff2(dff(:, 1)==2, 1) = 700;
dff2(dff(:, 1)==3, 1) = 900;
dff2(dff(:, 2)==1, 2) = 58;
dff2(dff(:, 2)==2, 2) = 70;
dff2(dff(:, 3)==1, 3) = 0.65;
dff2(dff(:, 3)==2, 3) = 1;
dff2(dff(:, 3)==3, 3) = 1.25;
dff2(dff(:, 4)==1, 4) = 3;
dff2(dff(:, 4)==2, 4) = 4;
dff2(dff(:, 4)==3, 4) = 5;


% Import max load results
filename = 'C:\Users\manpan\Documents\LTU\Polygoner\maxforcedispldata-NM.txt';
delimiter = '\t';
formatSpec = '%*s%*s%*s%*s%f%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
maxload = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray;

% Plot the interactions
figure;
interactionplot(maxload', dff2, 'varnames', parameters);
