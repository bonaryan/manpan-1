function [curves, shapes] = Polygoner_CFSM_out(profiles)
% Function that is called for a giver 3D cell array with profile coordinate
% data, executes CUFSM and returns the curves and shapes

% Size of the output matrix
matrix_size = size(profiles);

% Initialise the curves and shapes cell arrays to host the results
curves = cell(matrix_size(1), matrix_size(2), matrix_size(3));
shapes = cell(matrix_size(1), matrix_size(2), matrix_size(3));

% Define constants
E = 210000;
v = 0.3;
G = E/(2*(1+v));

% Define conditions
springs = 0;
constraints = 0;
BC = 'S-S';

% Define other parameters (to be explained)
GBTcon = struct('glob', 0, 'dist', 0, 'local', 0 , 'other', 0, 'ospace', [1], 'couple', [1], 'orth', [2]);
n = 100;
lengths = logspace(0, 3, n);
m_all = num2cell(ones(1, n));
neigs = 10;

% Loop executing CUFSM
for i = [1:matrix_size(1)];
    for j = [1:matrix_size(2)];
        for k= [1:matrix_size(3)];
            
            % Current profile
            c_prof = profiles{i, j, k}';
            
            % Number of vertices on the current profile
            l_prof = length(c_prof);
            
            % A column of ones
            col1 = ones(l_prof', 1);
            
            % Construct the 'node' array
            node = [(1:l_prof)', c_prof(:, 1), c_prof(:, 2), col1, col1, col1, col1, 100*col1];
            
            % Construct the 'elem' array
            elem = [(1:l_prof-1)', (1:l_prof-1)', (2:l_prof)', 0.1*ones(l_prof-1', 1), 100*ones(l_prof-1', 1)];
            
            % Construct the 'prop' array
            prop = [100, E, E, v, v, G ];
            
            % Run the FSM strip analysis
            [curves{i, j, k}, shapes{i, j, k}] =strip(prop, node, elem, lengths, springs, constraints, GBTcon, BC, m_all, neigs);
            
        end
    end
end