%function [joint_mean, joint_inv_cov, conditional_mean_A, conditional_mean_b, conditional_inv_cov]...
function [joint_mean, joint_inv_cov] = train_joint_normal_features(...
    cuboids, transformations, num_features, num_local_points)

num_cuboids = length(cuboids);
assert(length(transformations) == num_cuboids);

cuboids_exist = cell(num_cuboids);
for i = 1:num_cuboids
    cuboids_exist{i} = (sum(isnan(cuboids{i}), 2) == 0);
end    


%%
joint_mean = cell(num_cuboids, num_cuboids);
joint_inv_cov = cell(num_cuboids, num_cuboids);
% conditional_mean_A = cell(num_cuboids, num_cuboids);
% conditional_mean_b = cell(num_cuboids, num_cuboids);
% conditional_inv_cov = cell(num_cuboids, num_cuboids);

for i = 1:num_cuboids
    for j = 1:num_cuboids
        if i == j; continue; end;
        
        assert(size(cuboids{i}, 1) == size(cuboids{j}, 1));
        assert(size(cuboids{i}, 2) == num_features);
        assert(size(cuboids{j}, 2) == num_features);

        cuboids_exist_ij = cuboids_exist{i} & cuboids_exist{j};
        cuboids_i = cuboids{i}(cuboids_exist_ij, :);
        cuboids_j = cuboids{j}(cuboids_exist_ij, :);


        assert(length(transformations{i}) == size(cuboids{i}, 1));
        transformation_i = transformations{i}(cuboids_exist_ij, :);
        transformation_j = transformations{j}(cuboids_exist_ij, :);
        num_objects = length(transformation_i);
        assert(num_objects == length(transformation_j));
        if num_objects == 0; continue; end;
        

        X_i = cuboids_i;
        X_j = cuboids_j;

        for object_index = 1:num_objects
            assert(size(transformation_i, 2) == 12);
            first_translation_i = transformation_i(object_index, 10:12)';
            second_rotation_i = reshape(transformation_i(object_index, 1:9), 3, 3);
            
            assert(size(transformation_j, 2) == 12);
            first_translation_j = transformation_j(object_index, 10:12)';
            second_rotation_j = reshape(transformation_j(object_index, 1:9), 3, 3);

            for value_index = 1:num_local_points
                start_index = 3 * (value_index - 1) + 1;
                end_index = 3 * value_index;

                X_j(object_index, start_index:end_index) = ...
                    (first_translation_i + X_j(object_index, start_index:end_index)')';

                X_j(object_index, start_index:end_index) = ...
                    (second_rotation_i * X_j(object_index, start_index:end_index)')';
                
                X_i(object_index, start_index:end_index) = ...
                    (first_translation_j + X_i(object_index, start_index:end_index)')';

                X_i(object_index, start_index:end_index) = ...
                    (second_rotation_j * X_i(object_index, start_index:end_index)')';
            end
        end

        assert(sum(sum(X_i(:, 3 * num_local_points+1:end)...
            - cuboids_i(:, 3 * num_local_points+1:end))) == 0);
        
        assert(sum(sum(X_j(:, 3 * num_local_points+1:end)...
            - cuboids_j(:, 3 * num_local_points+1:end))) == 0);

        X = [X_i, X_j];
        mean_X = mean(X);

        centered_X = X - repmat(mean_X, num_objects, 1);
        cov = (centered_X' * centered_X) / num_objects;
        %inv_cov = pca_inv_cov(cov);
        %inv_cov = inv(cov + 1.0E-3 * eye(size(cov)));
        inv_cov = graphicalLasso(cov + 1.0E-3 * eye(size(cov)), 1.0E-3);
        
        joint_mean{i}{j} = mean_X';
        joint_inv_cov{i}{j} = inv_cov;

        %%
        diff = (X - repmat(mean_X, num_objects, 1))';
        error = diag(diff' * inv_cov * diff);
        disp(['(', num2str(i), ', ', num2str(j), '): max_error = ', num2str(max(error))]);

        csvwrite(['test_', num2str(i-1), '_', num2str(j-1), '.csv'], inv_cov);
        inv_cov = csvread(['test_', num2str(i-1), '_', num2str(j-1), '.csv']);
        
        %%
%         mean_1 = mean_X(1:num_features)';
%         mean_2 = mean_X(num_features+1:end)';
%         inv_cov_11 = inv_cov(1:num_features, 1:num_features);
%         inv_cov_22 = inv_cov(num_features+1:end, num_features+1:end);
% 
%         mean_s1 = mean_X(3*num_local_points+1:num_features)';
%         mean_s2 = mean_X(num_features+3*num_local_points+1:end)';
%         inv_cov_2s1 = inv_cov(num_features+1:end, 3*num_local_points+1:num_features);
%         inv_cov_1s2 = inv_cov(1:num_features, num_features+3*num_local_points+1:end);
% 
%         inv_inv_cov_11 = inv(inv_cov_11);
%         inv_inv_cov_22 = inv(inv_cov_22);
% 
%         cond_21_mean_A = -inv_inv_cov_22 * inv_cov_2s1;
%         cond_21_mean_b = mean_2 - cond_21_mean_A * mean_s1;
%         cond_21_inv_cov = inv_cov_22;
% 
%         conditional_mean_A{i}{j} = cond_21_mean_A;
%         conditional_mean_b{i}{j} = cond_21_mean_b;
%         conditional_inv_cov{i}{j} = cond_21_inv_cov;
% 
%         cond_12_mean_A = -inv_inv_cov_11 * inv_cov_1s2;
%         cond_12_mean_b = mean_1 - cond_12_mean_A * mean_s2;
%         cond_12_inv_cov = inv_cov_11;
% 
%         conditional_mean_A{j}{i} = cond_12_mean_A;
%         conditional_mean_b{j}{i} = cond_12_mean_b;
%         conditional_inv_cov{j}{i} = cond_12_inv_cov;
    end
end

