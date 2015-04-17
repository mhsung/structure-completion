#include "MeshCuboidParameters.h"

DEFINE_bool(param_optimize_with_non_linear_constraints, true, "");

DEFINE_int32(param_num_sample_point_neighbors, 8, "");
DEFINE_int32(param_min_num_cuboid_sample_points, 10, "");
DEFINE_int32(param_num_cuboid_surface_points, 1000, "");
DEFINE_int32(param_intra_cuboid_symmetry_axis, 0, "");
DEFINE_int32(param_eval_num_neighbor_range_samples, 1001, "");
DEFINE_int32(param_opt_max_iterations, 5, "");

DEFINE_double(param_min_sample_point_confidence, 0.7, "");
DEFINE_double(param_min_num_confidence_tol_sample_points, 0.5, "");
DEFINE_double(param_sample_point_neighbor_distance, 0.02, "");
DEFINE_double(param_min_cuboid_bbox_size, 0.1, "");
DEFINE_double(param_min_cuboid_bbox_diag_length, 0.3, "");
DEFINE_double(param_observed_point_radius, 1.0E-2, "");
DEFINE_double(param_min_cuboid_overall_visibility, 0.8, "");
DEFINE_double(param_max_potential, 1.0E8, "");
DEFINE_double(param_dummy_potential, 1.0E4, "");
DEFINE_double(param_eval_max_neighbor_range, 0.005, "");
DEFINE_double(param_opt_single_energy_term_weight, 1.0E4, "");
DEFINE_double(param_opt_symmetry_energy_term_weight, 1.0E6, "");

//DEFINE_double(param_sim_abs_attr_tol, 0.2, "");
//DEFINE_double(param_zero_tol, 1.0E-6, "");