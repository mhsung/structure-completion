#include "MeshCuboidParameters.h"

// -- Parameters -- //
//
DEFINE_bool(param_optimize_with_non_linear_constraints, true, "");
DEFINE_bool(param_optimize_training_cuboids, true, "");

DEFINE_int32(param_num_sample_point_neighbors, 8, "");
DEFINE_int32(param_min_num_cuboid_sample_points, 10, "");
DEFINE_int32(param_min_num_symmetric_point_pairs, 50, "");
DEFINE_int32(param_num_cuboid_surface_points, 1000, "");
DEFINE_int32(param_intra_cuboid_symmetry_axis, 0, "");
DEFINE_int32(param_eval_num_neighbor_range_samples, 1001, "");
DEFINE_int32(param_opt_max_iterations, 5, "");

DEFINE_double(param_min_sample_point_confidence, 0.7, "");
DEFINE_double(param_min_num_confidence_tol_sample_points, 0.5, "");
DEFINE_double(param_min_cuboid_bbox_size, 0.1, "");
DEFINE_double(param_min_cuboid_bbox_diag_length, 0.3, "");
DEFINE_double(param_sparse_neighbor_distance, 0.05, "");
DEFINE_double(param_cuboid_split_neighbor_distance, 0.2, "");
DEFINE_double(param_occlusion_test_neighbor_distance, 0.01, "");
DEFINE_double(param_eval_min_neighbor_distance, 0.02, "");
DEFINE_double(param_eval_max_neighbor_distance, 0.05, "");
DEFINE_double(param_min_cuboid_overall_visibility, 0.8, "");
DEFINE_double(param_max_potential, 1.0E8, "");
DEFINE_double(param_dummy_potential, 1.0E4, "");
DEFINE_double(param_null_cuboid_probability, 0.1, "");
DEFINE_double(param_fusion_visibility_smoothing_prior, 1.0, "");
DEFINE_double(param_fusion_grid_size, 0.01, "");
DEFINE_double(param_opt_single_energy_term_weight, 1.0E4, "");
DEFINE_double(param_opt_symmetry_energy_term_weight, 1.0E6, "");
DEFINE_double(param_part_assembly_window_size, 0.1, "");
DEFINE_double(param_part_assembly_voxel_size, 0.01, "");
DEFINE_double(param_part_assembly_voxel_variance, 0.01, "");

//DEFINE_double(param_sim_abs_attr_tol, 0.2, "");
//DEFINE_double(param_zero_tol, 1.0E-6, "");
//

// Test 2D view plane mask for occlusion.
DEFINE_bool(use_view_plane_mask, false, "");
DEFINE_double(param_view_plane_mask_proportion, 0.3, "");
DEFINE_double(param_view_plane_mask_min_x, 0.0, "");
DEFINE_double(param_view_plane_mask_min_y, 0.0, "");
DEFINE_double(param_view_plane_mask_max_x, 0.0, "");
DEFINE_double(param_view_plane_mask_max_y, 0.0, "");
// ---- //


// -- Experiments -- //
//
DEFINE_bool(run_ground_truth_cuboids, false, "");
DEFINE_bool(run_training, false, "");
DEFINE_bool(run_prediction, false, "");
DEFINE_bool(run_part_assembly, false, "");
DEFINE_bool(run_symmetry_detection, false, "");
DEFINE_bool(run_baseline, false, "");
DEFINE_bool(run_render_assembly, false, "");
DEFINE_bool(run_render_output, false, "");
DEFINE_bool(run_render_evaluation, false, "");
DEFINE_bool(run_extract_symmetry_info, false, "");
DEFINE_bool(no_evaluation, true, "");
DEFINE_string(mesh_filename, "", "");

DEFINE_string(data_root_path, "D:/Data/shape2pose/", "");
DEFINE_string(label_info_path, "data/0_body/dataset_name/", "");
DEFINE_string(mesh_path, "data/1_input/dataset_name/off/", "");
DEFINE_string(sample_path, "data/2_analysis/dataset_name/points/even1000/", "");
DEFINE_string(dense_sample_path, "data/2_analysis/dataset_name/points/random100000/", "");
DEFINE_string(mesh_label_path, "data/1_input/dataset_name/gt/", "");
DEFINE_string(sample_label_path, "data/4_experiments/exp1_dataset_name/1_prediction/", "");

DEFINE_string(retrieval_label_info_path, "data/0_body/dataset_name/", "");
DEFINE_string(retrieval_mesh_path, "data/1_input/dataset_name/off/", "");
DEFINE_string(retrieval_sample_path, "data/2_analysis/dataset_name/points/even1000/", "");
DEFINE_string(retrieval_dense_sample_path, "data/2_analysis/dataset_name/points/random100000/", "");
DEFINE_string(retrieval_mesh_label_path, "data/1_input/dataset_name/gt/", "");
DEFINE_string(retrieval_sample_label_path, "data/4_experiments/exp1_dataset_name/1_prediction/", "");

DEFINE_string(output_dir, "output", "");
DEFINE_string(training_dir, "training", "");
DEFINE_string(part_assembly_dir, "part_assembly", "");
DEFINE_string(symmetry_detection_dir, "symmetry_detection", "");
DEFINE_string(baseline_dir, "baseline", "");

DEFINE_string(label_info_filename, "regions.txt", "");
DEFINE_string(label_symmetry_info_filename, "regions_symmetry.txt", "");
DEFINE_string(symmetry_group_info_filename, "symmetry_groups.txt", "");
DEFINE_string(pose_filename, "pose.txt", "");
DEFINE_string(occlusion_pose_filename, "", "");
//DEFINE_string(occlusion_pose_filename, "occlusion_pose.txt", "");

DEFINE_string(single_feature_filename_prefix, "single_feature_", "");
DEFINE_string(pair_feature_filename_prefix, "pair_feature_", "");
DEFINE_string(single_stats_filename_prefix, "single_stats_", "");
DEFINE_string(pair_stats_filename_prefix, "pair_stats_", "");
DEFINE_string(feature_filename_prefix, "feature_", "");
DEFINE_string(transformation_filename_prefix, "transformation_", "");
DEFINE_string(joint_normal_relation_filename_prefix, "joint_normal_", "");
DEFINE_string(cond_normal_relation_filename_prefix, "conditional_normal_", "");
DEFINE_string(object_list_filename, "object_list.txt", "");

DEFINE_double(occlusion_test_radius, 0.01, "");
DEFINE_int32(random_view_seed, 20150416, "");

// To be removed.
//DEFINE_bool(use_symmetric_group_cuboids, false, "");
//
// ---- //
