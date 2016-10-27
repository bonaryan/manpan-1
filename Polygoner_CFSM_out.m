% Script that creates and saves a matlab database in the right format for
% CUFSM xrom a given set of x, y arrays

curves = {100,100,100}
shapes = {100,100,100}

matrix_size = size(profiles);

for i = [1:matrix_size(1)];
    for j = [1:matrix_size(2)]
        for k= [1:matrix_size(3)];
            elem = [according to the 'profiles']
            node = [according to the 'profiles']
            node = [(1:length(profiles{i, j, k}))', profiles{1, 1, 1}(1, :)', profiles{1, 1, 1}(2, :)', ones(length(profiles{1, 1, 1})', 1), ones(length(profiles{1, 1, 1})', 1), ones(length(profiles{1, 1, 1})', 1), ones(length(profiles{1, 1, 1})', 1)]
%             [curves, shapes]=strip('the function has to be run for the values taken from profiles cell array')
%               you need to construct 'elem' and 'node' arrays from the
%               'profiles' and feed them to the strip function.
        end
    end
end;