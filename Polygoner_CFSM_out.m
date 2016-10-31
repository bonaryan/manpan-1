% Script that creates and saves a matlab database in the right format for
% CUFSM xrom a given set of x, y arrays


curves = cell(size(profiles, 1),size(profiles, 2),size(profiles, 3));
shapes = cell(size(profiles, 1),size(profiles, 2),size(profiles, 3));

matrix_size = size(profiles);

E = 210000;
v = 0.3;
springs = 0;
constraints = 0;
GBTcon = struct('glob', 0, 'dist', 0, 'local', 0 , 'other', 0, 'ospace', [1], 'couple', [1], 'orth', [2]);
BC = 'S-S';
n = 100;
lengths = logspace(0, 3, n);
m_all = ones(1, n);
m_all = mat2cell(m_all, 1, ones(n, 1));

neigs = 10;


for i = [1:matrix_size(1)];
    for j = [1:matrix_size(2)];
        for k= [1:matrix_size(3)];
            %             elem = [according to the 'profiles']
            %             node = [according to the 'profiles']
            node = [(1:length(profiles{i, j, k}))', profiles{i, j, k}(1, :)', profiles{i, j, k}(2, :)', ones(length(profiles{i, j, k})', 1), ones(length(profiles{i, j, k})', 1), ones(length(profiles{i, j, k})', 1), ones(length(profiles{i, j, k})', 1)];
            elem = [(1:length(profiles{i, j, k})-1)', (1:length(profiles{i, j, k})-1)', (2:length(profiles{i, j, k}))', ones(length(profiles{i, j, k})-1', 1), ones(length(profiles{i, j, k})-1', 1)];
            prop = [100, E, E, v, v, E/(2*(1+v))];
            %             [curves(i, j, k), shapes(i, j, k)] = strip()
            %             you need to construct 'elem' and 'node' arrays from the
            %              'profiles' and feed them to the strip function.
            [curves(i, j, k), shapes(i, j, k)] = strip(prop,node,elem,lengths,springs,constraints,GBTcon,BC,m_all,neigs)
            %             matrix_size{i, j, k} = [curve; shapes];
        end
    end
end;