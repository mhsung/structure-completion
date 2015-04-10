#ifndef _ICP_H_
#define _ICP_H_

// Implemented based on the following paper:
// Olga Sorkine, "Least-Squares Rigid Motion Using SVD".

#define MAX_BBOX_ORIENTATION_NUM_ITERATIONS		100
#define MIN_BBOX_ORIENTATION_ANGLE_DIFFERENCE	1
#define MIN_BBOX_ORIENTATION_TRANSLATION		1.0E-6

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

	double compute_rigid_transformation(const Eigen::MatrixXd &_X, const Eigen::MatrixXd &_Y,
		Eigen::Matrix3d &_rotation_mat, Eigen::Vector3d &_translation_vec);

	// NOTE:
	// '_X' and '_Y' are modified.
	double run_iterative_closest_points(Eigen::MatrixXd &_X, Eigen::MatrixXd &_Y,
		Eigen::Matrix3d &_rotation_mat, Eigen::Vector3d &_translation_vec);
}

#endif	// _ICP_H_