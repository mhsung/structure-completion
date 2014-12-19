#include "MeshCuboidNonLinearSolver.h"

#include "MeshCuboid.h"
#include "MeshCuboidRelation.h"
#include "MeshCuboidPredictor.h"
#include "MeshCuboidSolver.h"
#include "Utilities.h"

#include "ceres/ceres.h"

#include <Eigen/Geometry> 


CostFunctor::CostFunctor(
	std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor& _predictor,
	const double _quadprog_ratio,
	QGLWidget *_viewer)
	: cuboids_(_cuboids)
	, predictor_(_predictor)
	, quadprog_ratio_(_quadprog_ratio)
	, viewer_(_viewer)
{
}

double *CostFunctor::get_cuboid_attributes(
	const std::vector<MeshCuboid *>& _cuboids)
{
	unsigned int num_cuboids = static_cast<unsigned int>(_cuboids.size());

	double *x = new double[num_cuboids * k_num_attributes];

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		MeshCuboid *cuboid = _cuboids[cuboid_index];

		MyMesh::Point bbox_center = cuboid->get_bbox_center();
		MyMesh::Normal bbox_size = cuboid->get_bbox_size();
		std::array<MyMesh::Normal, 3> bbox_axes = cuboid->get_bbox_axes();

		Eigen::Matrix3d bbox_rotation_mat;
		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
			for (unsigned int i = 0; i < 3; ++i)
				bbox_rotation_mat.row(axis_index)(i) = bbox_axes[axis_index][i];

		Eigen::AngleAxisd angle_axis;
		angle_axis.fromRotationMatrix(bbox_rotation_mat);
		double angle = angle_axis.angle();
		Eigen::Vector3d angle_axis_vec = angle * angle_axis.axis();
			
		for (unsigned int i = 0; i < 3; ++i)
		{
			x[cuboid_index * k_num_attributes + k_center_index + i] = bbox_center[i];
			x[cuboid_index * k_num_attributes + k_size_index + i] = bbox_size[i];
			x[cuboid_index * k_num_attributes + k_rotation_index + i] = angle_axis_vec[i];
		}
	}

	return x;
}

void CostFunctor::set_cuboid_attributes(
	const double* const x, std::vector<MeshCuboid *>& _cuboids)
{
	unsigned int num_cuboids = static_cast<unsigned int>(_cuboids.size());
	//std::cout << "Length of array = " << (sizeof(x) / sizeof(*x)) << std::endl;

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		MeshCuboid *cuboid = _cuboids[cuboid_index];

		MyMesh::Point bbox_center;
		MyMesh::Normal bbox_size;
		Eigen::Vector3d angle_axis_vec;

		for (unsigned int i = 0; i < 3; ++i)
		{
			bbox_center[i] = x[cuboid_index * k_num_attributes + k_center_index + i];
			bbox_size[i] = x[cuboid_index * k_num_attributes + k_size_index + i];
			angle_axis_vec[i] = x[cuboid_index * k_num_attributes + k_rotation_index + i];
		}

		double angle = angle_axis_vec.norm();
		Eigen::AngleAxisd angle_axis(angle, angle_axis_vec.normalized());
		if (angle == 0) angle_axis.axis() = Eigen::Vector3d::UnitX();
		Eigen::Matrix3d bbox_rotation_mat = angle_axis.toRotationMatrix();

		std::array<MyMesh::Normal, 3> bbox_axes;
		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
			for (unsigned int i = 0; i < 3; ++i)
				bbox_axes[axis_index][i] = bbox_rotation_mat.row(axis_index)(i);

		cuboid->set_bbox_center(bbox_center);
		cuboid->set_bbox_size(bbox_size, false);
		cuboid->set_bbox_axes(bbox_axes, false);

		cuboid->update_corner_points();
		std::array<MyMesh::Point, MeshCuboid::k_num_corners> new_bbox_corners =
			cuboid->get_bbox_corners();

		// Update cuboid surface points.
		for (unsigned int point_index = 0; point_index < cuboid->num_cuboid_surface_points();
			++point_index)
		{
			MeshCuboidSurfacePoint *cuboid_surface_point =
				cuboid->get_cuboid_surface_point(point_index);

			MyMesh::Point new_point(0.0);
			for (unsigned int corner_index = 0;
				corner_index < MeshCuboid::k_num_corners; ++corner_index)
				new_point += cuboid_surface_point->corner_weights_[corner_index] *
					new_bbox_corners[corner_index];

			cuboid_surface_point->point_ = new_point;
		}
	}
}

bool CostFunctor::operator()(double const* const* x, double* residual) const
{
	// Update cuboids.
	set_cuboid_attributes(x[0], cuboids_);

	// Compute error.
	double single_total_energy = 0, pair_total_energy = 0;
	get_optimization_error(cuboids_, predictor_, single_total_energy, pair_total_energy);

	Real total_energy = pair_total_energy + quadprog_ratio_ * single_total_energy;
	residual[0] = total_energy;

	if (viewer_) viewer_->updateGL();

	return true;
}

void non_linear_optimize_cuboids(
	std::vector<MeshCuboid *>& _cuboids,
	const MeshCuboidPredictor& _predictor,
	const double _quadprog_ratio,
	QGLWidget *_viewer)
{
	unsigned int num_cuboids = static_cast<unsigned int>(_cuboids.size());
	const int num_parameters = CostFunctor::k_num_attributes * num_cuboids;

	// The variable to solve for with its initial value. It will be
	// mutated in place by the solver.
	double *x = CostFunctor::get_cuboid_attributes(_cuboids);
	std::vector<double*> parameter_blocks;
	parameter_blocks.push_back(x);

	// Build the problem.
	ceres::Problem problem;

	// Set up the only cost function (also known as residual). This uses
	// numeric differentiation to obtain the derivative (jacobian).
	CostFunctor cost_functor(_cuboids, _predictor, _quadprog_ratio, _viewer);

	ceres::DynamicNumericDiffCostFunction<CostFunctor> cost_function(&cost_functor);
	cost_function.AddParameterBlock(num_parameters);
	cost_function.SetNumResiduals(1);
	problem.AddResidualBlock(&cost_function, NULL, parameter_blocks);

	// Run the solver!
	ceres::Solver::Options options;
	options.max_num_iterations = 100;
	options.minimizer_type = ceres::LINE_SEARCH;
	options.line_search_direction_type = ceres::STEEPEST_DESCENT;
	options.linear_solver_type = ceres::DENSE_QR;
	options.minimizer_progress_to_stdout = true;
	ceres::Solver::Summary summary;
	Solve(options, &problem, &summary);
	//std::cout << summary.BriefReport() << "\n";
	std::cout << summary.FullReport() << "\n";
	//std::cout << "x : " << initial_x << " -> " << x << "\n";

	CostFunctor::set_cuboid_attributes(x, _cuboids);
	delete[] x;
}