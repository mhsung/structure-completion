#ifndef _MESH_CUBOID_PARAMETERS_H_
#define _MESH_CUBOID_PARAMETERS_H_

#include <gflags/gflags.h>

// -- Parameters -- //
//
DECLARE_bool(param_optimize_with_non_linear_constraints);
DECLARE_bool(param_optimize_training_cuboids);

DECLARE_int32(param_num_sample_point_neighbors);
DECLARE_int32(param_min_num_cuboid_sample_points);
DECLARE_int32(param_min_num_symmetric_point_pairs);
DECLARE_int32(param_num_cuboid_surface_points);
DECLARE_int32(param_intra_cuboid_symmetry_axis);
DECLARE_int32(param_eval_num_neighbor_range_samples);
DECLARE_int32(param_opt_max_iterations);

DECLARE_double(param_min_sample_point_confidence);
DECLARE_double(param_min_num_confidence_tol_sample_points);
DECLARE_double(param_min_cuboid_bbox_size);
DECLARE_double(param_min_cuboid_bbox_diag_length);
DECLARE_double(param_sparse_neighbor_distance);
DECLARE_double(param_cuboid_split_neighbor_distance);
DECLARE_double(param_occlusion_test_neighbor_distance);
DECLARE_double(param_eval_max_neighbor_distance);
DECLARE_double(param_min_cuboid_overall_visibility);
DECLARE_double(param_max_potential);
DECLARE_double(param_dummy_potential);
DECLARE_double(param_opt_single_energy_term_weight);
DECLARE_double(param_opt_symmetry_energy_term_weight);
DECLARE_double(param_part_assembly_window_size);
DECLARE_double(param_part_assembly_voxel_size);
DECLARE_double(param_part_assembly_voxel_variance);

//DECLARE_double(param_sim_abs_attr_tol);
//DECLARE_double(param_zero_tol);
//
// ---- //


// -- Experiments -- //
//
DECLARE_bool(run_ground_truth_cuboids);
DECLARE_bool(run_training);
DECLARE_bool(run_prediction);
DECLARE_bool(run_part_assembly);
DECLARE_bool(run_symmetry_detection);
DECLARE_bool(run_baseline);
DECLARE_bool(run_assembly_render);
DECLARE_string(mesh_filename);

DECLARE_string(data_root_path);
DECLARE_string(label_info_path);
DECLARE_string(mesh_path);
DECLARE_string(sample_path);
DECLARE_string(dense_sample_path);
DECLARE_string(mesh_label_path);
DECLARE_string(sample_label_path);
DECLARE_string(output_dir);
DECLARE_string(training_dir);
DECLARE_string(part_assembly_dir);
DECLARE_string(symmetry_detection_dir);
DECLARE_string(baseline_dir);

DECLARE_string(label_info_filename);
DECLARE_string(label_symmetry_info_filename);
DECLARE_string(symmetry_group_info_filename);
DECLARE_string(pose_filename);
DECLARE_string(occlusion_pose_filename);
//DECLARE_string(occlusion_pose_filename);

DECLARE_string(single_feature_filename_prefix);
DECLARE_string(pair_feature_filename_prefix);
DECLARE_string(single_stats_filename_prefix);
DECLARE_string(pair_stats_filename_prefix);
DECLARE_string(feature_filename_prefix);
DECLARE_string(transformation_filename_prefix);
DECLARE_string(joint_normal_relation_filename_prefix);
DECLARE_string(cond_normal_relation_filename_prefix);
DECLARE_string(object_list_filename);

DECLARE_double(occlusion_test_radius);
DECLARE_int32(random_view_seed);

// To be removed.
//DECLARE_bool(use_symmetric_group_cuboids, false, "");
//
// ---- //

#endif // _MESH_CUBOID_PARAMETERS_H_