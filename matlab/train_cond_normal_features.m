function [conditional_mean_A, conditional_mean_b, conditional_inv_cov]...
    = train_cond_normal_features(cuboids, transformations, num_features, num_local_points)

num_cuboids = length(cuboids);
assert(length(transformations) == num_cuboids);

cuboids_exist = cell(num_cuboids);
for i = 1:num_cuboids
    cuboids_exist{i} = (sum(isnan(cuboids{i}), 2) == 0);
end    


%%
conditional_mean_A = cell(num_cuboids, num_cuboids);
conditional_mean_b = cell(num_cuboids, num_cuboids);
conditional_inv_cov = cell(num_cuboids, num_cuboids);

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
		num_objects = length(transformation_i);

		X_i = cuboids_i(:, 3 * num_local_points+1:end);
		X_j = cuboids_j;

		for object_index = 1:num_objects
			assert(size(transformation_i, 2) == 12);
			first_translation = transformation_i(object_index, 10:12)';
			second_rotation = reshape(transformation_i(object_index, 1:9), 3, 3);

			for value_index = 1:num_local_points
				start_index = 3 * (value_index - 1) + 1;
				end_index = 3 * value_index;

				X_j(object_index, start_index:end_index) = ...
					(first_translation + X_j(object_index, start_index:end_index)')';

				X_j(object_index, start_index:end_index) = ...
					(second_rotation * X_j(object_index, start_index:end_index)')';
			end
		end
		
		assert(sum(sum(X_j(:, 3 * num_local_points+1:end)...
			- cuboids_j(:, 3 * num_local_points+1:end))) == 0);

		mean_i = mean(X_i);
		mean_j = mean(X_j);
		
		centered_X_i = X_i - repmat(mean_i, num_objects, 1);
		centered_X_j = X_j - repmat(mean_j, num_objects, 1);
		
		cov_ii = (centered_X_i' * centered_X_i) / num_objects;
		cov_jj = (centered_X_j' * centered_X_j) / num_objects;
		cov_ij = (centered_X_i' * centered_X_j) / num_objects;
		cov_ji = (centered_X_j' * centered_X_i) / num_objects;
		%inv_cov_ii = pca_inv_cov(cov_ii);
		inv_cov_ii = inv(cov_ii + 1.0E-3 * eye(size(cov_ii)));
		
		conditional_mean_A{i}{j} = cov_ji * inv_cov_ii;
		conditional_mean_b{i}{j} = mean_j' - conditional_mean_A{i}{j} * mean_i';
		conditional_cov_ji = cov_jj - cov_ji * inv_cov_ii * cov_ij;
		%conditional_inv_cov_ji = pca_inv_cov(conditional_cov_ji);
		conditional_inv_cov_ji = inv(conditional_cov_ji + 1.0E-3 * eye(size(conditional_cov_ji)));
		
		conditional_inv_cov{i}{j} = conditional_inv_cov_ji;
		
		%%
		mean_ij = conditional_mean_A{i}{j} * X_i' + repmat(conditional_mean_b{i}{j}, 1, num_objects);
		diff = X_j' - mean_ij;
		error = diag(diff' * conditional_inv_cov{i}{j} * diff);
		disp(['(', num2str(i), ', ', num2str(j), '): max_error = ', num2str(max(error))]);
		
		%csvwrite(['error_', num2str(i-1), '_', num2str(j-1), '.csv'], error);
	end
end

