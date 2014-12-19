function [] = train_dataset_1(path, file_num, is_test_data)

if exist([path, '/obj/'], 'dir') ~= 7; mkdir([path, '/obj/']); end;
if exist([path, '/pts/'], 'dir') ~= 7; mkdir([path, '/pts/']); end;
if exist([path, '/comp/'], 'dir') ~= 7; mkdir([path, '/comp/']); end;
if exist([path, '/arff/'], 'dir') ~= 7; mkdir([path, '/arff/']); end;

cuboid_1 = create_cuboid;
cuboid_2 = create_cuboid;

param_11 = max(rand, 0.1);
param_12 = max(rand, 0.1);

if is_test_data
    param_21 = max(rand, 0.1);
    param_22 = max(rand, 0.1);
else
    param_21 = param_11;
    param_22 = param_12;
end

%% Scale
S_1 = eye(3);
S_1(1, 1) = param_11;
S_1(2, 2) = param_12;
S_1(3, 3) = 0.1;
cuboid_1 = S_1 * cuboid_1;

S_2 = eye(3);
S_2(1, 1) = param_21;
S_2(2, 2) = param_22;
S_2(3, 3) = 0.1;
cuboid_2 = S_2 * cuboid_2;

%% Translation
t_1 = zeros(3, 1);
t_1(2) = +param_12 + 0.5;
cuboid_1 = repmat(t_1, 1, size(cuboid_1, 2)) + cuboid_1;

t_2 = zeros(3, 1);
t_2(2) = -param_22 - 0.5;
cuboid_2 = repmat(t_2, 1, size(cuboid_2, 2)) + cuboid_2;


%% Random Rotation (cuboid_2)
angle = (2 * rand - 1) * pi / 12;
R = [cos(angle), -sin(angle), 0;...
    sin(angle), cos(angle), 0;...
    0, 0, 1];

cuboid_2 = R * cuboid_2;


%% Random Rotation
angle = (2 * rand - 1) * pi / 12;
R = [cos(angle), -sin(angle), 0;...
    sin(angle), cos(angle), 0;...
    0, 0, 1];

cuboid_1 = R * cuboid_1;
cuboid_2 = R * cuboid_2;

%% Random Translation
t = [rand - 0.5; rand - 0.5; 0];
cuboid_1 = repmat(t, 1, size(cuboid_1, 2)) + cuboid_1;
cuboid_2 = repmat(t, 1, size(cuboid_2, 2)) + cuboid_2;


%%
num_faces_1 = size(cuboid_1, 2) / 3;
num_faces_2 = size(cuboid_2, 2) / 3;

samples_1 = create_cuboid_samples(cuboid_1, 0);
samples_2 = create_cuboid_samples(cuboid_2, num_faces_1);

save_cuboid_obj([path, '/obj/', num2str(file_num)], [cuboid_1, cuboid_2]);
save_samples_pts([path, '/pts/', num2str(file_num)], [samples_1, samples_2]);

save_cuboid_labels([path, '/comp/', num2str(file_num)], [size(cuboid_1, 2)/3, size(cuboid_2, 2)/3]);
save_sample_labels([path, '/arff/', num2str(file_num)], [size(samples_1, 2), size(samples_2, 2)]);