close all;
clear all;
% ---- %
num_features = 36;
num_local_points = 9;
% ---- %
feature_filename_prefix = 'feature_';
transformation_filename_prefix = 'transformation_';

[cuboids, transformations] = load_features(...
    feature_filename_prefix, transformation_filename_prefix, num_features);

[joint_mean, joint_inv_cov]...
    = train_joint_normal_features(cuboids, transformations, num_features, num_local_points);
save_joint_normal_features_csv;

[conditional_mean_A, conditional_mean_b, conditional_inv_cov] = ...
    train_cond_normal_features(cuboids, transformations, num_features, num_local_points);
save_cond_normal_features_csv;