#ifndef _MESH_CUBOID_NON_LINEAR_SOLVER_H_
#define _MESH_CUBOID_NON_LINEAR_SOLVER_H_

#include "MeshCuboid.h"
#include "MeshCuboidRelation.h"
#include "MeshCuboidPredictor.h"

Eigen::VectorXd solve_quadratic_programming_with_constraints(
	const std::vector<MeshCuboid *>& _cuboids,
	const Eigen::MatrixXd& _quadratic_term,
	const Eigen::VectorXd& _linear_term,
	const double _constant_term,
	Eigen::VectorXd* _init_values_vec = NULL);

#endif	// _MESH_CUBOID_NON_LINEAR_SOLVER_H_