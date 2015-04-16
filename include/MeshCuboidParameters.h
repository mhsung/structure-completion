#ifndef _MESH_CUBOID_PARAMETERS_H_
#define _MESH_CUBOID_PARAMETERS_H_

#include <gflags/gflags.h>

DECLARE_bool(param_optimize_with_non_linear_constraints);

DECLARE_int32(param_num_sample_point_neighbors);
DECLARE_int32(param_min_num_cuboid_sample_points);
DECLARE_int32(param_num_cuboid_surface_points);
DECLARE_int32(param_intra_cuboid_symmetry_axis);

DECLARE_double(param_min_sample_point_confidence);
DECLARE_double(param_min_num_confidence_tol_sample_points);
DECLARE_double(param_sample_point_neighbor_distance);
DECLARE_double(param_min_cuboid_bbox_size);
DECLARE_double(param_min_cuboid_bbox_diag_length);
DECLARE_double(param_observed_point_radius);
DECLARE_double(param_min_cuboid_overall_visibility);
DECLARE_double(param_max_potential);
DECLARE_double(param_dummy_potential);

DECLARE_double(param_eval_max_neighbor_range);
DECLARE_double(param_eval_num_neighbor_range_samples);

//DECLARE_double(param_sim_abs_attr_tol);
//DECLARE_double(param_zero_tol);

#endif // _MESH_CUBOID_PARAMETERS_H_