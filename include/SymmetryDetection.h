#ifndef _SYMMETRY_DETECTION_H_
#define _SYMMETRY_DETECTION_H_

#include <list>
#include <vector>
#include <ANN/ANN.h>
#include <Eigen/Core>
#include "ICP.h"


namespace SymmetryDetection {
	typedef struct ReflectionPlane
	{
		Eigen::Vector3d normal_;
		Eigen::Vector3d point_;
		std::vector<unsigned int> inlier_indices_;
	};

	void detect_reflectional_symmetry(
		const Eigen::MatrixXd &_points,
		const double &_inlier_dist, const double &_inlier_ratio,
		std::list<ReflectionPlane> &_reflection_planes);

	void detect_reflectional_symmetry(
		const Eigen::MatrixXd &_sparse_points,
		const Eigen::MatrixXd &_dense_points,
		const double &_inlier_dist, const double &_inlier_ratio,
		std::list<ReflectionPlane> &_reflection_planes);
}

#endif	// _SYMMETRY_DETECTION_H_