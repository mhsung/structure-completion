#ifndef _MESH_CUBOID_NON_LINEAR_SOLVER_H_
#define _MESH_CUBOID_NON_LINEAR_SOLVER_H_

#define GLOG_NO_ABBREVIATED_SEVERITIES

#include "MeshCuboid.h"
#include "MeshCuboidRelation.h"
#include "MeshCuboidPredictor.h"

#include <vector>
#include <QtOpenGL/qgl.h>

class CostFunctor {
public:
	static const int k_center_index = 0;
	static const int k_size_index = 3;
	static const int k_rotation_index = 6;
	static const int k_num_attributes = 9;

	CostFunctor(std::vector<MeshCuboid *>& _cuboids,
		const MeshCuboidPredictor& _predictor,
		const double _quadprog_ratio,
		QGLWidget *_viewer);

	static double *get_cuboid_attributes(
		const std::vector<MeshCuboid *>& _cuboids);

	static void set_cuboid_attributes(
		const double *x, std::vector<MeshCuboid *>& _cuboids);

	bool operator()(double const* const* x, double* residual) const;

private:
	std::vector<MeshCuboid *>& cuboids_;
	const MeshCuboidPredictor& predictor_;
	const double quadprog_ratio_;
	QGLWidget *viewer_;
};

void non_linear_optimize_cuboids(
	std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor &_predictor,
	const double _quadprog_ratio = 1.0,
	QGLWidget *_viewer = NULL);

#endif	// _MESH_CUBOID_SOLVER_H_