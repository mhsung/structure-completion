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
	ANNkd_tree* create_kd_tree(const Eigen::MatrixXd &_points, ANNpointArray &_ann_points);

	// _X_points: Query points.
	// _Y_values: Data values.
	// _closest_Y_value: Associated Y values for each X (based on the closest distance between X and Y).
	// _Y_ann_kd_tree: KD-tree of Y points. If it is not given, it is created by considering _Y as Y points.
	void get_closest_points(const Eigen::MatrixXd &_X_points, const Eigen::MatrixXd &_Y_values,
		Eigen::MatrixXd &_closest_Y_values, ANNkd_tree *_Y_ann_kd_tree = NULL);

	double compute_rigid_transformation(const Eigen::MatrixXd &_X, const Eigen::MatrixXd &_Y,
		Eigen::Matrix3d &_rotation_mat, Eigen::Vector3d &_translation_vec);

	// NOTE:
	// '_X' and '_Y' are modified.
	double run_iterative_closest_points(Eigen::MatrixXd &_X, Eigen::MatrixXd &_Y,
		Eigen::Matrix3d &_rotation_mat, Eigen::Vector3d &_translation_vec);
}

#endif	// _ICP_H_