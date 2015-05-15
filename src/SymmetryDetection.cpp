#include "SymmetryDetection.h"
#include "simplerandom.h"

Eigen::VectorXd get_symmetric_point(
	const Eigen::VectorXd& _normal,
	const Eigen::VectorXd& _center,
	const Eigen::VectorXd& _point)
{
	// Assume that '_normal' is normalized.
	Eigen::VectorXd plane_to_point = _normal * _normal.dot(_point - _center);
	Eigen::VectorXd symmetric_point = _point - (plane_to_point * 2);
	return symmetric_point;
}

Eigen::MatrixXd get_symmetric_point(
	const Eigen::VectorXd& _normal,
	const Eigen::VectorXd& _center,
	const Eigen::MatrixXd& _points)
{
	// Assume that '_normal' is normalized.
	Eigen::MatrixXd plane_to_points = _normal * _normal.transpose() * (_points.colwise() - _center);
	Eigen::MatrixXd symmetric_point = _points - (plane_to_points * 2);

	// Debug.
	//Eigen::MatrixXd mid_points = 0.5 * (symmetric_point + _points);
	//Eigen::VectorXd distances = (mid_points.colwise() - _center).transpose() * _normal;
	//std::cout << distances.sum() / distances.rows() << std::endl;

	return symmetric_point;
}

void SymmetryDetection::detect_reflectional_symmetry(
	const Eigen::MatrixXd &_points,
	const double &_inlier_dist, const double &_inlier_ratio,
	std::list<ReflectionPlane> &_reflection_planes)
{
	detect_reflectional_symmetry(_points, _points, _inlier_dist, _inlier_ratio, _reflection_planes);
}

void SymmetryDetection::detect_reflectional_symmetry(
	const Eigen::MatrixXd &_sparse_points,
	const Eigen::MatrixXd &_dense_points,
	const double &_inlier_dist, const double &_inlier_ratio,
	std::list<ReflectionPlane> &_reflection_planes)
{
	const unsigned int dimension = 3;
	assert(_sparse_points.rows() == dimension);
	assert(_dense_points.rows() == dimension);
	assert(_inlier_ratio > 0.0);
	assert(_inlier_ratio < 1.0);

	_reflection_planes.clear();
	const unsigned int num_sparse_points = _sparse_points.cols();
	const unsigned int num_dense_points = _dense_points.cols();

	ANNpointArray sparse_ann_points;
	ANNkd_tree* sparse_kd_tree = ICP::create_kd_tree(_sparse_points, sparse_ann_points);

	const unsigned int num_iteration = 1000;

	static SimpleRandomCong_t rng_cong;
	// seed = num_iteration.
	simplerandom_cong_seed(&rng_cong, num_iteration);

	for (unsigned int iteration = 0; iteration < num_iteration; ++iteration)
	{
		double sparse_index_1 = simplerandom_cong_next(&rng_cong) % num_sparse_points;
		double sparse_index_2 = simplerandom_cong_next(&rng_cong) % num_sparse_points;

		Eigen::VectorXd sparse_point_1 = _sparse_points.col(sparse_index_1);
		Eigen::VectorXd sparse_point_2 = _sparse_points.col(sparse_index_2);

		ReflectionPlane reflection_plane;

		reflection_plane.normal_ = sparse_point_2 - sparse_point_1;
		reflection_plane.normal_.normalize();
		reflection_plane.point_ = 0.5 * (sparse_point_1 + sparse_point_2);

		Eigen::MatrixXd symmetric_sparse_points = get_symmetric_point(
			reflection_plane.normal_, reflection_plane.point_, _sparse_points);

		Eigen::VectorXd sparse_distances;
		ICP::get_closest_points(sparse_kd_tree, symmetric_sparse_points, sparse_distances);
		unsigned int num_sparse_inliers = (sparse_distances.array() <= _inlier_dist).count();

		if (static_cast<double>(num_sparse_inliers) / num_sparse_points >= _inlier_ratio)
		{
			reflection_plane.inlier_indices_.clear();
			for (unsigned int i = 0; i < num_sparse_points; ++i)
				if (sparse_distances[i] < _inlier_dist)
					reflection_plane.inlier_indices_.push_back(i);
			
			_reflection_planes.push_back(reflection_plane);

			std::cout << "Reflection plane detected [" << _reflection_planes.size() << "]:" << std::endl;
			std::cout << " - # of inliers = (" << num_sparse_inliers << " / " << num_sparse_points << ")" << std::endl;
			std::cout << " - Normal = (" << reflection_plane.normal_.transpose() << ")" << std::endl;
			std::cout << " - Point = (" << reflection_plane.point_.transpose() << ")" << std::endl;
			std::cout << std::endl;
		}
	}

	annDeallocPts(sparse_ann_points);
	delete sparse_kd_tree;
}

/*
		if (static_cast<double>(num_sparse_inliers) / num_sparse_points >= _inlier_ratio)
		{
			Eigen::MatrixXd symmetric_dense_points = get_symmetric_point(
				reflection_plane.normal_, reflection_plane.point_, _dense_points);

			Eigen::VectorXd dense_distances = ((symmetric_dense_points - _dense_points).colwise().norm()).transpose();
			unsigned int num_dense_inliers = (dense_distances.array() <= _inlier_dist).count();
			Eigen::MatrixXd X(dimension, num_dense_inliers), Y(dimension, num_dense_inliers);

			reflection_plane.inlier_indices_.clear();
			reflection_plane.inlier_indices_.reserve(num_dense_inliers);

			unsigned int count_inliers = 0;
			for (unsigned int i = 0; i < num_dense_points; ++i)
			{
				if (dense_distances[i] < _inlier_dist)
				{
					X.col(count_inliers) = symmetric_dense_points.col(i);
					Y.col(count_inliers) = _dense_points.col(i);
					reflection_plane.inlier_indices_.push_back(i);
					++count_inliers;
				}
			}
			assert(count_inliers == num_dense_inliers);

			Eigen::Matrix3d rotation_mat;
			Eigen::Vector3d translation_vec;
			double threshold = 0.1 * _inlier_ratio;
			ICP::run_iterative_closest_points(X, Y, rotation_mat, translation_vec, &threshold);

			reflection_plane.normal_ = rotation_mat * reflection_plane.normal_ + translation_vec;
			reflection_plane.point_ = rotation_mat * reflection_plane.point_ + translation_vec;

			std::cout << "Optimized..." << std::endl;
			std::cout << " - Normal = (" << reflection_plane.normal_.transpose() << ")" << std::endl;
			std::cout << " - Point = (" << reflection_plane.point_.transpose() << ")" << std::endl;

			std::cout << std::endl;

			_reflection_planes.push_back(reflection_plane);
		}
		*/
