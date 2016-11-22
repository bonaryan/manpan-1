% Single profile solver
% Script that runs the strip analysis for a single profile

% Input profile number
i = 1;
j = 1;
k = 1;

% Define constants
E = 210000;
v = 0.3;
G = E/(2*(1+v));

% Define conditions
% % springs = 0;
BC = 'S-S';

% Define other parameters (to be explained)
GBTcon = struct('glob', 0, 'dist', 0, 'local', 0 , 'other', 0, 'ospace', [1], 'couple', [1], 'orth', [2]);
n = 100;
l = 4.5;
lengths = logspace(0, l, n);
lengths = [lengths, 50000];
m_all = num2cell(ones(1, (n+1)));
neigs = 10;
            
% Selected profile xy
c_prof1 = profiles{i, j, k}';

% Selected profile plate thickness
t = meta{i, j, k}(2);

% Number of vertices on the selected profile
l_prof = length(c_prof1);

% A column of ones
col1 = ones(l_prof', 1);

% % Construct the 2 extra parts by rotating the imported one
% % R2 = [cos(-2*pi/3), -sin(-2*pi/3); sin(-2*pi/3), cos(-2*pi/3)];
% % R3 = [cos(2*pi/3), -sin(2*pi/3); sin(2*pi/3), cos(2*pi/3)];
% % for a = 1:l_prof;
% %     c_prof2(a, :) = (R2*c_prof1(a, :)')';
% %     c_prof3(a, :) = (R3*c_prof1(a, :)')';
% % end

% Construct the 'node' array
% % node = [(1:l_prof)', c_prof1(:, 1), c_prof1(:, 2), col1, col1, col1, col1, 100*col1;
% %     (1*l_prof+1:2*l_prof)', c_prof2(:, 1), c_prof2(:, 2), col1, col1, col1, col1, 100*col1
% %     (2*l_prof+1:3*l_prof)', c_prof3(:, 1), c_prof3(:, 2), col1, col1, col1, col1, 100*col1];

node = [(1:l_prof)', c_prof1(:, 1), c_prof1(:, 2), col1, col1, col1, col1, 100*col1];


% Construct the 'elem' array
% % elem = [(1:l_prof-1)', (1:l_prof-1)', (2:l_prof)', t*ones(l_prof-1', 1), 100*ones(l_prof-1', 1);
% %     (1*l_prof:2*l_prof-2)', l_prof+(1:l_prof-1)', l_prof+(2:l_prof)', t*ones(l_prof-1', 1), 100*ones(l_prof-1', 1);
% %     (2*l_prof-1:3*l_prof-3)', 2*l_prof+(1:l_prof-1)', 2*l_prof+(2:l_prof)', t*ones(l_prof-1', 1), 100*ones(l_prof-1', 1)];

    elem = [(1:l_prof-1)', (1:l_prof-1)', (2:l_prof)', t*ones(l_prof-1', 1), 100*ones(l_prof-1', 1)];
    

% Construct the 'prop' array
prop = [100, E, E, v, v, G ];

% Constructing general constraints, springs between the sectors
% % constraints = [l_prof+2 1 1.000 l_prof-1 1 0.000 0 0
% %     l_prof+2 2 1.000 l_prof-1 2 0.000 0 0
% %     l_prof+2 3 1.000 l_prof-1 3 0.000 0 0
% %     2*l_prof+2 1 1.000 2*l_prof-1 1 0.000 0 0
% %     2*l_prof+2 2 1.000 2*l_prof-1 2 0.000 0 0
% %     2*l_prof+2 3 1.000 2*l_prof-1 3 0.000 0 0
% %     2 1 1.000 3*l_prof-1 1 0.000 0 0
% %     2 2 1.000 3*l_prof-1 2 0.000 0 0
% %     2 3 1.000 3*l_prof-1 3 0.000 0 0];

constraints = 0;
ku = 5;    
kw = 5;
    kv = 5;
%     Springs
    springs = [2 1  ku 0 
    2  2.000 kw 0 
    2   3.00  kv 0
    l_prof-1 1 ku 0
    l_prof-1 2 kw 0
    l_prof-1 3 kv 0];

% Save the input variables from this run to a file. This file can be loaded
% in CUFSM gui
save('loadfile.mat', 'prop', 'node', 'elem', 'lengths', 'springs', 'constraints', 'GBTcon', 'BC', 'm_all', 'neigs');

% Run the FSM strip analysis
[curves, shapes] =strip(prop, node, elem, lengths, springs, constraints, GBTcon, BC, m_all, neigs);

% Clear unwanted variables
clear c_prof1 c_prof2 c_prof3 col1 E G i j k l l_prof n R2 R3 t v a

% Clear the rest of variables
clear prop node elem lengths springs constraints GBTcon BC m_all neigs