#ifndef _ICP_H_
#define _ICP_H_

// Implemented based on the following paper:
// Olga Sorkine, "Least-Squares Rigid Motion Using SVD".

#define MIN_NUM_ICP_POINT_PAIRS		10
#define MAX_NUM_ICP_ITERATIONS		1000
#define MIN_ICP_ANGLE_DIFFERENCE	1.0E-2
#define MIN_ICP_TRANSLATION			1.0E-8

#include <ANN/ANN.h>
#include <Eigen/Core>


namespace ICP {
	ANNkd_tree* create_kd_tree(
		const Eigen::MatrixXd &_points,
		ANNpointArray &_ann_points);

	// NOTE:
	// In all matrices, "column" is an instance.
	
	// _data_ann_kd_tree: KD-tree of data points. If it is not given,
	//    it is created by considering 'data_values' as data points.
	// _query_points: Query points.
	// _data_values: Data values.
	// _closest_data_value: Associated data values for each query (based on the closest distance between query and data).
	template<typename T>
	void get_closest_points(
		ANNkd_tree *_data_ann_kd_tree,
		const Eigen::MatrixXd &_query_points,
		const Eigen::MatrixBase<T> &_data_values,
		Eigen::MatrixBase<T> &_closest_data_values);

	void get_closest_points(
		ANNkd_tree *_data_ann_kd_tree,
		const Eigen::MatrixXd &_query_points,
		Eigen::VectorXd &_distances);

	// Return: error (Minus error means that the computation is failed.)
	double compute_rigid_transformation(const Eigen::MatrixXd &_X, const Eigen::MatrixXd &_Y,
		Eigen::Matrix3d &_rotation_mat, Eigen::Vector3d &_translation_vec,
		const double *_distance_threshold = NULL);

	// Return: error (Minus error means that the computation is failed.)
	// '_X' is modified.
	double run_iterative_closest_points(Eigen::MatrixXd &_X, const Eigen::MatrixXd &_Y,
		Eigen::Matrix3d &_rotation_mat, Eigen::Vector3d &_translation_vec,
		const double *_distance_threshold = NULL);
}

#endif	// _ICP_H_