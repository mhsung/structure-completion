%function [] = generate_dataset_2(path, file_num, is_test_data)

rand('seed', 1);

feature_0 = [];
feature_1 = [];

transformation_0 = [];
transformation_1 = [];

for i = 1:100
    param_c = 0.01;
    param_x = 0.1 * rand + 0.01;
    param_y = 0.1 * rand + 0.01;

    x0 = [0, +(param_c + param_y), 0,...
        -param_x, +(param_c + param_y)-param_y, -param_c,...
        +param_x, +(param_c + param_y)-param_y, -param_c,...
        -param_x, +(param_c + param_y)+param_y, -param_c,...
        +param_x, +(param_c + param_y)+param_y, -param_c,...
        -param_x, +(param_c + param_y)-param_y, +param_c,...
        +param_x, +(param_c + param_y)-param_y, +param_c,...
        -param_x, +(param_c + param_y)+param_y, +param_c,...
        +param_x, +(param_c + param_y)+param_y, +param_c,...
        0, -param_c, -param_c, -param_c, -param_c,...
        +param_c, +param_c, +param_c, +param_c];

    x1 = [0, -(param_c + param_y), 0,...
        -param_x, -(param_c + param_y)-param_y, -param_c,...
        +param_x, -(param_c + param_y)-param_y, -param_c,...
        -param_x, -(param_c + param_y)+param_y, -param_c,...
        +param_x, -(param_c + param_y)+param_y, -param_c,...
        -param_x, -(param_c + param_y)-param_y, +param_c,...
        +param_x, -(param_c + param_y)-param_y, +param_c,...
        -param_x, -(param_c + param_y)+param_y, +param_c,...
        +param_x, -(param_c + param_y)+param_y, +param_c,...
        0, -param_c, -param_c, -param_c, -param_c,...
        +param_c, +param_c, +param_c, +param_c];
    
    angle = (2 * rand - 1)  * (0 / 180) * pi;
    %Column-wise.
    r0 = [cos(angle),sin(angle),0,...
        -sin(angle),cos(angle),0,...
        0,0,1];
    
    x0_center = x0(1:3)';
    x0_corners = reshape(x0(4:27), 3, 8) - repmat(x0_center, 1, 8);
    x0_corners = reshape(r0, 3, 3) * x0_corners;
    x0_corners = x0_corners + repmat(x0_center, 1, 8);
    x0(4:27) = reshape(x0_corners, 1, 24);
    
    %Column-wise.
    r1 = [1,0,0,...
        0,1,0,...
        0,0,1];

    % Add one garbage feature at the end.
    feature_0 = [feature_0; [x0, 0]];
    feature_1 = [feature_1; [x1, 0]];
    
    % Add one garbage feature at the end.
    transformation_0 = [transformation_0; [r0, -x0(1:3), 0]];
    transformation_1 = [transformation_1; [r1, -x1(1:3), 0]];
end

avg_feature_0 = sum(feature_0, 1)./size(feature_0, 1);
avg_feature_1 = sum(feature_1, 1)./size(feature_1, 1);

% Add one garbage feature at the end.
transformed_avg_feature_0 = avg_feature_0 + [repmat(-avg_feature_1(1:3), 1, 9), zeros(1, 10)];
transformed_avg_feature_1 = avg_feature_1 + [repmat(-avg_feature_0(1:3), 1, 9), zeros(1, 10)];

csvwrite('feature_0.csv', feature_0);
csvwrite('feature_1.csv', feature_1);

csvwrite('transformation_0.csv', transformation_0);
csvwrite('transformation_1.csv', transformation_1);

% csvwrite('transformed_avg_feature_0.csv', transformed_avg_feature_0');
% csvwrite('transformed_avg_feature_1.csv', transformed_avg_feature_1');

csvwrite('cuboids.csv', [avg_feature_0(:, 4:27); avg_feature_1(:, 4:27)]);


%%
center = avg_feature_0(1:3);

corner_0 = avg_feature_0(4:6);
corner_1 = avg_feature_0(7:9);
corner_2 = avg_feature_0(10:12);
corner_3 = avg_feature_0(13:15);

corner_4 = avg_feature_0(16:18);
corner_5 = avg_feature_0(19:21);
corner_6 = avg_feature_0(22:24);
corner_7 = avg_feature_0(25:27);

x = corner_0(1);

min_y = corner_0(2);
max_y = corner_2(2);
dy = (max_y - min_y) / 10;

min_z = corner_4(3);
max_z = corner_0(3);
dz = (max_z - min_z) / 10;

[Y, Z] = meshgrid(min_y:dy:max_y, min_z:dz:max_z);
Y = reshape(Y, size(Y, 1)*size(Y, 2), 1);
Z = reshape(Z, size(Z, 1)*size(Z, 2), 1);
samples = [repmat(x, length(Y), 1), Y, Z];

file = fopen(['points.pts'], 'w');
for i = 1:size(samples, 2)
    fprintf(file, '%d %f %f %f %f %f %f\n',...
        0, 0, 0, 0,...
        samples(1, i), samples(2, i), samples(3, i));
end
fclose(file);

save_sample_labels('point_labels', [size(samples, 1)]);

