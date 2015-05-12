#include "ICP.h"

#include <assert.h>
#include <Eigen/Geometry>
#include <Eigen/LU> 
#include <Eigen/SVD>


namespace ICP {
	ANNkd_tree* create_kd_tree(
		const Eigen::MatrixXd &_points,
		ANNpointArray &_ann_points)
	{
		assert(_points.rows() == 3);
		assert(_points.cols() > 0);

		int dimension = 3;
		int num_points = static_cast<int>(_points.cols());
		ANNkd_tree* ann_kd_tree;

		_ann_points = annAllocPts(num_points, dimension);	// allocate data points

		// read data points
		for (int point_index = 0; point_index < num_points; point_index++)
		{
			for (int i = 0; i < dimension; i++)
				_ann_points[point_index][i] = _points.col(point_index)[i];
		}

		ann_kd_tree = new ANNkd_tree(		// build search structure
			_ann_points,					// the data points
			num_points,						// number of points
			dimension);						// dimension of space

		return ann_kd_tree;
	}

	template<typename T>
	void get_closest_points(
		ANNkd_tree *_data_ann_kd_tree,
		const Eigen::MatrixXd &_query_points,
		const Eigen::MatrixBase<T> &_data_values,
		Eigen::MatrixBase<T> &_closest_data_values)
	{
		assert(_data_ann_kd_tree->nPoints() > 0);
		assert(_query_points.rows() == 3);
		assert(_query_points.cols() > 0);
		assert(_data_values.rows() > 0);
		assert(_data_values.cols() == _data_ann_kd_tree->nPoints());

		int num_queries = static_cast<int>(_query_points.cols());
		int num_data = static_cast<int>(_data_values.cols());
		int dimension = static_cast<int>(_data_values.rows());
		_closest_data_values = Eigen::MatrixBase<T>::Zero(dimension, num_queries);

		ANNpoint q = annAllocPt(3);
		ANNidxArray nn_idx = new ANNidx[1];
		ANNdistArray dd = new ANNdist[1];

		for (int point_index = 0; point_index < num_queries; point_index++)
		{
			for (unsigned int i = 0; i < 3; ++i)
				q[i] = _query_points.col(point_index)(i);

			_data_ann_kd_tree->annkSearch(q, 1, nn_idx, dd);
			int closest_Y_point_index = nn_idx[0];
			assert(closest_Y_point_index < num_data);
			assert(point_index < _closest_data_values.cols());

			_closest_data_values.col(point_index) = _data_values.col(closest_Y_point_index);
		}

		// Deallocate ANN.
		annDeallocPt(q);
		delete[] nn_idx;
		delete[] dd;
	}

	void get_closest_points(
		ANNkd_tree *_data_ann_kd_tree,
		const Eigen::MatrixXd &_query_points,
		Eigen::VectorXd &_distances)
	{
		assert(_data_ann_kd_tree->nPoints() > 0);
		assert(_query_points.rows() == 3);
		assert(_query_points.cols() > 0);

		int num_queries = static_cast<int>(_query_points.cols());
		_distances = Eigen::VectorXd::Zero(num_queries);

		ANNpoint q = annAllocPt(3);
		ANNidxArray nn_idx = new ANNidx[1];
		ANNdistArray dd = new ANNdist[1];

		for (int point_index = 0; point_index < num_queries; point_index++)
		{
			for (unsigned int i = 0; i < 3; ++i)
				q[i] = _query_points.col(point_index)(i);

			_data_ann_kd_tree->annkSearch(q, 1, nn_idx, dd);
			assert(point_index < _distances.rows());
			_distances[point_index] = std::sqrt(dd[0]);
		}

		// Deallocate ANN.
		annDeallocPt(q);
		delete[] nn_idx;
		delete[] dd;
	}

	double compute_rigid_transformation(const Eigen::MatrixXd &_X, const Eigen::MatrixXd &_Y,
		Eigen::Matrix3d &_rotation_mat, Eigen::Vector3d &_translation_vec,
		const double *_distance_threshold)
	{
		// Return: error (Minus error means that the computation is failed.)

		assert(_X.rows() == 3);
		assert(_Y.rows() == 3);

		if (_X.cols() == 0 || _Y.cols() == 0)
		{
			std::cerr << "Warning: No point exists." << std::endl;
			return -1;
		}

		Eigen::MatrixXd::Index X_num_points = _X.cols();
		Eigen::MatrixXd::Index Y_num_points = _Y.cols();

		ANNpointArray X_ann_points;
		ANNkd_tree* X_ann_kd_tree = create_kd_tree(_X, X_ann_points);

		ANNpointArray Y_ann_points;
		ANNkd_tree* Y_ann_kd_tree = create_kd_tree(_Y, Y_ann_points);
		
		Eigen::MatrixXd closest_X; 
		Eigen::MatrixXd closest_Y;

		get_closest_points(X_ann_kd_tree, _Y, _X, closest_X);
		get_closest_points(Y_ann_kd_tree, _X, _Y, closest_Y);

		annDeallocPts(X_ann_points); delete X_ann_kd_tree;
		annDeallocPts(Y_ann_points); delete Y_ann_kd_tree;

		assert(closest_X.cols() == Y_num_points);
		assert(closest_Y.cols() == X_num_points);

		Eigen::MatrixXd::Index n = X_num_points + Y_num_points;
		Eigen::MatrixXd all_X(3, n);
		Eigen::MatrixXd all_Y(3, n);

		all_X << _X, closest_X;
		all_Y << closest_Y, _Y;

		assert(all_X.cols() == n);
		assert(all_Y.cols() == n);

		Eigen::MatrixXd diff = all_X - all_Y;
		Eigen::RowVectorXd squared_dists = (diff.array() * diff.array()).colwise().sum();

		// Hausdorff distance.
		double prev_error = std::sqrt(squared_dists.maxCoeff());

		// Partial ICP.
		if (_distance_threshold)
		{
			//std::cout << "max_squared_dists = " << squared_dists.maxCoeff() << std::endl;
			//std::cout << "avg_squared_dists = " << squared_dists.sum() / n << std::endl;

			double squared_distance_threshold = (*_distance_threshold) * (*_distance_threshold);
			Eigen::MatrixXd::Index subset_n = (squared_dists.array() <= squared_distance_threshold).count();
			Eigen::MatrixXd subset_all_X(3, subset_n);
			Eigen::MatrixXd subset_all_Y(3, subset_n);
			Eigen::MatrixXd subset_diff(3, subset_n);
			Eigen::RowVectorXd subset_squared_dists(subset_n);

			int count = 0;
			for (int i = 0; i < n; ++i)
			{
				if (squared_dists[i] <= squared_distance_threshold)
				{
					assert(count < subset_n);
					subset_all_X.col(count) = all_X.col(i);
					subset_all_Y.col(count) = all_Y.col(i);
					subset_diff.col(count) = diff.col(i);
					subset_squared_dists.col(count) = squared_dists.col(i);
					++count;
				}
			}
			assert(count == subset_n);

			n = subset_n;
			all_X.swap(subset_all_X);
			all_Y.swap(subset_all_Y);
			diff.swap(subset_diff);
			squared_dists.swap(subset_squared_dists);
		}


		if (n < MIN_NUM_ICP_POINT_PAIRS)
		{
			std::cerr << "Warning: Too few point pairs." << std::endl;
			return -1;
		}

		Eigen::Vector3d all_X_mean = all_X.rowwise().mean();
		Eigen::Vector3d all_Y_mean = all_Y.rowwise().mean();

		Eigen::MatrixXd centered_all_X = all_X.colwise() - all_X_mean;
		Eigen::MatrixXd centered_all_Y = all_Y.colwise() - all_Y_mean;

		Eigen::MatrixXd S = centered_all_X * centered_all_Y.transpose();
		assert(S.rows() == 3);
		assert(S.cols() == 3);

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);

		Eigen::Matrix3d det_mat = Eigen::Matrix3d::Identity();
		det_mat(2, 2) = (svd.matrixV() * svd.matrixU().transpose()).determinant();

		_rotation_mat = svd.matrixV() * det_mat * svd.matrixU().transpose();
		_translation_vec = all_Y_mean - _rotation_mat * all_X_mean;

		diff = ((_rotation_mat * all_X).colwise() + _translation_vec) - all_Y;
		squared_dists = (diff.array() * diff.array()).colwise().sum();

		// Hausdorff distance.
		double error = std::sqrt(squared_dists.maxCoeff());
		if (error > prev_error)
			return -1;

		return error;
	}

	double run_iterative_closest_points(Eigen::MatrixXd &_X, const Eigen::MatrixXd &_Y,
		Eigen::Matrix3d &_rotation_mat, Eigen::Vector3d &_translation_vec,
		const double *_distance_threshold)
	{
		_rotation_mat = Eigen::Matrix3d::Identity();
		_translation_vec = Eigen::Vector3d::Zero();
		double prev_error = std::numeric_limits<double>::max();

		const unsigned int max_num_iterations = MAX_NUM_ICP_ITERATIONS;
		const double min_angle_difference = static_cast<double>(MIN_ICP_ANGLE_DIFFERENCE) / 180.0 * M_PI;
		const double min_translation = MIN_ICP_TRANSLATION;

		for (unsigned int iteration = 0; iteration < max_num_iterations; ++iteration)
		{
			Eigen::Matrix3d new_rotation_mat;
			Eigen::Vector3d new_translation_vec;
			double error = compute_rigid_transformation(_X, _Y,
				new_rotation_mat, new_translation_vec, _distance_threshold);

			if (error < 0)
			{
				// The transformation computation is failed.
				if (iteration == 0) return -1;
				//std::cout << "[FINAL] error = " << prev_error << std::endl;
				return prev_error;
			}
			else if (error > prev_error)
			{
				// NOTE:
				// Continue only when the error value is decreasing.
				//std::cout << "[FINAL] error = " << prev_error << std::endl;
				return prev_error;
			}

			prev_error = error;
			//std::cout << "error = " << prev_error << std::endl;

			Eigen::AngleAxisd rotation_angle;
			rotation_angle.fromRotationMatrix(new_rotation_mat);
			if (rotation_angle.angle() < min_angle_difference && new_translation_vec.norm() < min_translation)
				return prev_error;
			
			_rotation_mat = new_rotation_mat * _rotation_mat;
			_translation_vec = new_rotation_mat * _translation_vec + new_translation_vec;

			_X = ((new_rotation_mat * _X).colwise() + new_translation_vec);
		}

		return prev_error;
	}

}
