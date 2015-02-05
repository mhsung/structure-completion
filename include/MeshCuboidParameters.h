#ifndef _MESH_CUBOID_PARAMETERS_H_
#define _MESH_CUBOID_PARAMETERS_H_

#include <gflags/gflags.h>

DECLARE_int32(param_num_sample_point_neighbors);
DECLARE_int32(param_min_num_cuboid_sample_points);
DECLARE_int32(param_num_cuboid_surface_points);

DECLARE_double(param_sample_point_confidence_tol);
DECLARE_double(param_min_num_confidence_tol_sample_points);
DECLARE_double(param_sample_point_neighbor_distance);
DECLARE_double(param_min_cuboid_bbox_size);
DECLARE_double(param_min_cuboid_bbox_diag_length);
DECLARE_double(param_observed_point_radius);
DECLARE_double(param_max_potential);
DECLARE_double(param_dummy_potential);
//DECLARE_double(param_sim_abs_attr_tol);
//DECLARE_double(param_zero_tol);

#endif // _MESH_CUBOID_PARAMETERS_H_