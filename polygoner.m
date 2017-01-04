function [profiles, meta] = polygoner(nrange, drange, slendrange, fy, rcoef, nbend, l_ratio, t_ratio, lambda);
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
% slendrange = linspace(80, 140, 10);
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
            
            % find the element properties on the current profile
            nele = size(elem,1);
            for v = 1:nele;
                sn = elem(v,1); fn = elem(v,2); 
                % thickness of the element
                tk(v,1) = t;
                % compute the coordinate of the mid point of the element
                xm(v) = mean(node([sn fn],1));
                ym(v) = mean(node([sn fn],2));
                % compute the dimension of the element
                xd(v) = diff(node([sn fn],1));
                yd(v) = diff(node([sn fn],2));
                % compute the length of the element
                L(v,1) = norm([xd(v) yd(v)]);
                Ao(v,1) = L(v)*tk(v);      
            end
                        
            % Calculating cross-sectional class and effective area if needed:
            for v = 1:nele;
                Ep = [Ao L tk];
                ne = length(elem);
                epsilon=sqrt(235/fy); Ep2=zeros(ne,2); lambdap=zeros(ne,1); ro=zeros(ne,1); 
                if Ep(v,1) == eps
                    Ep2(v,:)=[0 123];
                else
                    %EC3-1-5 Part 4.4
                    lambdap(v)=(Ep(v,2)/Ep(v,3))/(28.4*epsilon*2);
                    ro(v)=abs((lambdap(v)-0.055*4)/lambdap(v)^2);
                    if ro(v)>1
                        ro(v)=1;
                    end
                    %EC3-1-1 Table 5.2
                    if Ep(v,2)/Ep(v,3) <= 42*epsilon
                        Ep2(v,1)=Ep(v,1);
                        Ep2(v,2)=3;
                    else
                        Ep2(v,1)=Ep(v,1)*ro(v);
                        Ep2(v,2)=4;
                    end                     
                end
            end
            % compute the effective cross section area
            Aeff = sum(Ep2(v,1));
            Class = max(Ep2(v,2)); 
            
            
            % Classification according to EC3 1-1
            max_side = max(sqrt(diff(node(:, 2)).^2+diff(node(:, 1)).^2))
            if max_side/t <= 42*epsilon
                Class = 3;
            else
                Class = 4;
            end
            
            % Loop through the different member slendernesses. The 'meta'
            % array has one more dimension (4D)
            for l = 1:length(lambda);
                
                % Current profile length
                len = lambda*pi*sqrt(E*I/(A*fy));
                
                % Store the metadata in a cell array
                meta{i, j, k, l} = [drange(j); t; tg; fy; A; Iyy; Izz; len(l); Aeff; Class];    
            end
        end
    end
end


% Plot result
% plot(x, y);

% Save the profile database and metadata to the current directory as .mat
save('profiles.mat', 'profiles');
save('meta.mat', 'meta');