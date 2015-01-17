function [cuboids, transformations, num_objects] = load_features(...
    feature_filename_prefix, transformation_filename_prefix, num_features)

i = 0;
cuboids = {};
transformations = {};
num_objects = -1;

while true
    %%
    attributes_filename = [feature_filename_prefix, num2str(i), '.csv'];
    if exist(attributes_filename, 'file') ~= 2
        break;
    end
    display(['Loading "', attributes_filename '"...']);
    
    x = csvread(attributes_filename);
    if num_objects < 0
        num_objects = size(x, 1);
    else
        assert(num_objects == size(x, 1));
    end
    
    % NOTE:
    % Remove the last empty column.
    if sum(x(:, end)) == 0
        x = x(:, 1:end-1);
    end
    assert(size(x, 2) == num_features);
    
    %%
    transformation_filename = [transformation_filename_prefix, num2str(i), '.csv'];
    if exist(transformation_filename, 'file') ~= 2
        assert false;
    end
    display(['Loading "', transformation_filename '"...']);
    
    T = csvread(transformation_filename);
    assert(num_objects == size(T, 1));

    % NOTE:
    % Remove the last empty column.
    if sum(T(:, end)) == 0
        T = T(:, 1:end-1);
    end

    %%
    % NOTE:
    % Matlab index starts from 1.
    cuboids{i+1} = x;
    transformations{i+1} = T;
    
    %%
    i = i + 1;
end
