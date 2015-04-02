#include "MeshCuboidPredictor.h"

#include "ICP.h"
#include "MeshCuboidParameters.h"
#include "Utilities.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <ANN/ANN.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>


MeshCuboidPredictor::MeshCuboidPredictor(unsigned int _num_labels)
	: num_labels_(_num_labels)
{

}

Real MeshCuboidPredictor::get_single_potential(
	const MeshCuboid *_cuboid,
	const MeshCuboidAttributes *_attributes,
	const MeshCuboidTransformation *_transformation,
	const LabelIndex _label_index)const
{
	assert(_cuboid);
	assert(_attributes);
	assert(_transformation);

	Real potential = 0.0;

	// FIXME:
	// This should be done in pre-processing.
	unsigned int num_confidence_tol_sample_point = 0;
	for (std::vector<MeshSamplePoint *>::const_iterator sample_it = _cuboid->get_sample_points().begin();
		sample_it != _cuboid->get_sample_points().end(); ++sample_it)
	{
		assert(_label_index < (*sample_it)->label_index_confidence_.size());
		if ((*sample_it)->label_index_confidence_[_label_index] >= FLAGS_param_min_sample_point_confidence)
			++num_confidence_tol_sample_point;
	}

	if (num_confidence_tol_sample_point <
		FLAGS_param_min_num_confidence_tol_sample_points * _cuboid->get_sample_points().size())
	{
		potential = FLAGS_param_max_potential;
	}

	return potential;
}

Real MeshCuboidPredictor::get_pair_potential(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const MeshCuboidAttributes *_attributes_1, const MeshCuboidAttributes *_attributes_2,
	const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2)const
{
	assert(_cuboid_1); assert(_cuboid_2);
	assert(_attributes_1); assert(_attributes_2);
	assert(_transformation_1); assert(_transformation_2);

	// Not implemented.
	Real potential = 0.0;
	return potential;
}

void MeshCuboidPredictor::get_single_quadratic_form(
	MeshCuboid *_cuboid, const unsigned int _cuboid_index,
	Eigen::MatrixXd &_quadratic_term, Eigen::VectorXd &_linear_term, double& _constant_term) const
{
	assert(_cuboid);

	_quadratic_term.setZero();
	_linear_term.setZero();
	_constant_term = 0;

	unsigned int num_sample_points = _cuboid->num_sample_points();
	unsigned int num_cuboid_surface_points = _cuboid->num_cuboid_surface_points();
	const unsigned int mat_size = _quadratic_term.cols();
	const unsigned int num_corners = MeshCuboid::k_num_corners;
	const unsigned int num_attributes = MeshCuboidAttributes::k_num_attributes;

	// X: sample points, Y: cuboid surface points.
	unsigned int num_X_points = num_sample_points;
	unsigned int num_Y_points = num_cuboid_surface_points;

	if (num_X_points == 0 || num_Y_points == 0)
		return;

	double sum_visibility = 0.0;
	for (int cuboid_surface_point_index = 0;
		cuboid_surface_point_index < num_cuboid_surface_points;
		++cuboid_surface_point_index)
	{
		sum_visibility += _cuboid->get_cuboid_surface_point(cuboid_surface_point_index)->visibility_;
	}

	if (sum_visibility == 0)
		return;

	Eigen::MatrixXd all_X_points(3, num_X_points + num_Y_points);
	Eigen::MatrixXd all_Y_points(3, num_X_points + num_Y_points);
	Eigen::MatrixXd all_Y_coeffs(num_corners, num_X_points + num_Y_points);
	Eigen::VectorXd all_weights(num_X_points + num_Y_points);

	unsigned int pair_index = 0;


	// Sample point -> Cuboid surface point.
	for (int sample_point_index = 0; sample_point_index < num_sample_points; ++sample_point_index)
	{
		int cuboid_surface_point_index =
			_cuboid->get_sample_to_cuboid_surface_correspondences(sample_point_index);
		assert(cuboid_surface_point_index >= 0);

		//
		Eigen::Vector3d X_point;
		Eigen::Vector3d Y_point;
		Eigen::VectorXd Y_coeff(MeshCuboid::k_num_corners);

		for (unsigned int i = 0; i < 3; ++i)
		{
			X_point(i) = _cuboid->get_sample_point(sample_point_index)->point_[i];
			Y_point(i) = _cuboid->get_cuboid_surface_point(cuboid_surface_point_index)->point_[i];
		}

		for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
			Y_coeff(corner_index) = _cuboid->get_cuboid_surface_point(
			cuboid_surface_point_index)->corner_weights_[corner_index];

		all_X_points.col(pair_index) = X_point;
		all_Y_coeffs.col(pair_index) = Y_coeff;
		all_Y_points.col(pair_index) = Y_point;
		all_weights(pair_index) = 1.0 / num_sample_points;

		++pair_index;
		//
	}

	// Cuboid surface point -> Sample point.
	for (int cuboid_surface_point_index = 0;
		cuboid_surface_point_index < num_cuboid_surface_points;
		++cuboid_surface_point_index)
	{
		int sample_point_index =
			_cuboid->get_cuboid_surface_to_sample_correspondence(cuboid_surface_point_index);
		assert(sample_point_index >= 0);

		//
		Eigen::Vector3d X_point;
		Eigen::Vector3d Y_point;
		Eigen::VectorXd Y_coeff(MeshCuboid::k_num_corners);

		for (unsigned int i = 0; i < 3; ++i)
		{
			X_point(i) = _cuboid->get_sample_point(sample_point_index)->point_[i];
			Y_point(i) = _cuboid->get_cuboid_surface_point(cuboid_surface_point_index)->point_[i];
		}

		for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
			Y_coeff(corner_index) = _cuboid->get_cuboid_surface_point(
			cuboid_surface_point_index)->corner_weights_[corner_index];

		double visibility = _cuboid->get_cuboid_surface_point(cuboid_surface_point_index)->visibility_;

		all_X_points.col(pair_index) = X_point;
		all_Y_coeffs.col(pair_index) = Y_coeff;
		all_Y_points.col(pair_index) = Y_point;
		all_weights(pair_index) = visibility / sum_visibility;

		++pair_index;
		//
	}
	
	assert(pair_index == num_X_points + num_Y_points);

	assert(_quadratic_term.rows() == mat_size);
	assert(_linear_term.rows() == mat_size);
	assert(_cuboid_index * MeshCuboidAttributes::k_num_attributes <= mat_size);


	MeshCuboidAttributes attribute;
	attribute.compute_attributes(_cuboid);
	Eigen::VectorXd x = attribute.get_attributes();

	for (unsigned int point_index = 0; point_index < num_X_points + num_Y_points; ++point_index)
	{
		Eigen::Vector3d X_point = all_X_points.col(point_index);
		Eigen::VectorXd Y_coeff = all_Y_coeffs.col(point_index);
		Eigen::VectorXd Y_point = all_Y_points.col(point_index);
		assert(Y_coeff.size() == MeshCuboid::k_num_corners);
		double weight = all_weights(point_index);

		Eigen::MatrixXd A0(3, num_attributes);
		A0.setZero();

		for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
		{
			for (unsigned int i = 0; i < 3; ++i)
				A0.col(MeshCuboidAttributes::k_corner_index + 3 * corner_index + i)(i) =
				Y_coeff(corner_index);
		}

#ifdef DEBUG_TEST
		Real error = std::abs(((A0 * x) - Y_point).norm());
		CHECK_NUMERICAL_ERROR(__FUNCTION__, error);
#endif

		Eigen::Vector3d b = -X_point;

		Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3, mat_size);
		unsigned int start_index = (_cuboid_index * MeshCuboidAttributes::k_num_attributes);
		A.block<3, MeshCuboidAttributes::k_num_attributes>(0, start_index) = A0;

		_quadratic_term = _quadratic_term + weight * (A.transpose() * A);
		_linear_term = _linear_term + weight * (A.transpose() * b);
		_constant_term = _constant_term + weight * (b.transpose() * b)(0);
	}
}

Real MeshCuboidPredictor::get_pair_quadratic_form(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const unsigned int _cuboid_index_1, const unsigned int _cuboid_index_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2,
	Eigen::MatrixXd &_quadratic_term, Eigen::VectorXd &_linear_term, double& _constant_term)const
{
	assert(_cuboid_1); assert(_cuboid_2);
	assert(_label_index_1 < num_labels_);
	assert(_label_index_2 < num_labels_);

	const unsigned int num_attributes = MeshCuboidAttributes::k_num_attributes;
	const unsigned int mat_size = num_labels_ * num_attributes;

	_quadratic_term = Eigen::MatrixXd::Zero(mat_size, mat_size);
	_linear_term = Eigen::VectorXd::Zero(mat_size);
	_constant_term = 0.0;

	return 0;
}

MeshCuboidJointNormalRelationPredictor::MeshCuboidJointNormalRelationPredictor(
	const std::vector< std::vector<MeshCuboidJointNormalRelations *> > &_relations)
	: MeshCuboidPredictor(_relations.size())
	, relations_(_relations)
{
	for (unsigned int label_index = 0; label_index < num_labels_; ++label_index)
		assert(relations_[label_index].size() == num_labels_);
}

void MeshCuboidJointNormalRelationPredictor::get_missing_label_indices(
	const std::list<LabelIndex> &_given_label_indices,
	std::list<LabelIndex> &_missing_label_indices)const
{
	_missing_label_indices.clear();
	unsigned int num_labels = relations_.size();

	bool *is_missing_label = new bool[num_labels];
	memset(is_missing_label, true, num_labels * sizeof(bool));

	for (std::list<LabelIndex>::const_iterator it = _given_label_indices.begin();
		it != _given_label_indices.end(); ++it)
	{
		LabelIndex label_index = (*it);
		assert(label_index < num_labels);
		is_missing_label[label_index] = false;
	}

	for (LabelIndex label_index_1 = 0; label_index_1 < num_labels; ++label_index_1)
	{
		if (is_missing_label[label_index_1])
		{
			bool is_label_confliced = false;

			for (LabelIndex label_index_2 = 0; label_index_2 < num_labels; ++label_index_2)
			{
				if (!is_missing_label[label_index_2] && !relations_[label_index_1][label_index_2])
				{
					is_label_confliced = true;
					break;
				}
			}

			if (!is_label_confliced)
			{
				_missing_label_indices.push_back(label_index_1);
			}
		}
	}

	delete[] is_missing_label;
}

Real MeshCuboidJointNormalRelationPredictor::get_pair_potential(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const MeshCuboidAttributes *_attributes_1, const MeshCuboidAttributes *_attributes_2,
	const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2) const
{
	assert(_cuboid_1); assert(_cuboid_2);
	assert(_attributes_1); assert(_attributes_2);
	assert(_transformation_1); assert(_transformation_2);

	assert(_label_index_1 < num_labels_);
	assert(_label_index_2 < num_labels_);

	// NOTE:
	// Now considering only different label pairs.
	assert(_label_index_1 != _label_index_2);

	Real potential = FLAGS_param_max_potential;

	//const MeshCuboidJointNormalRelations *relation_12 = relations_[_label_index_1][_label_index_2];
	//const MeshCuboidJointNormalRelations *relation_21 = relations_[_label_index_2][_label_index_1];

	//if (relation_12 && relation_21)
	//{
	//	Real potential_12 = relation_12->compute_error(_cuboid_1, _cuboid_2, _transformation_1, _transformation_2);
	//	Real potential_21 = relation_21->compute_error(_cuboid_2, _cuboid_1, _transformation_2, _transformation_1);
	//	potential = potential_12 + potential_21;
	//}

	const MeshCuboidJointNormalRelations *relation_12 = relations_[_label_index_1][_label_index_2];
	if (relation_12)
	{
		Real potential_12 = relation_12->compute_error(_cuboid_1, _cuboid_2, _transformation_1, _transformation_2);
		potential = potential_12;

		// Using bilateral relations.
		// (A, B) and (B, A) pairs are the same.
#ifdef DEBUG_TEST
		const MeshCuboidJointNormalRelations *relation_21 = relations_[_label_index_2][_label_index_1];
		if (relation_21)
		{
			Real potential_21 = relation_12->compute_error(_cuboid_1, _cuboid_2, _transformation_1, _transformation_2);
			CHECK_NUMERICAL_ERROR(__FUNCTION__, potential_12, potential_21);
		}
	}
#endif

	return potential;
}

Real MeshCuboidJointNormalRelationPredictor::get_pair_quadratic_form(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const unsigned int _cuboid_index_1, const unsigned int _cuboid_index_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2,
	Eigen::MatrixXd &_quadratic_term, Eigen::VectorXd &_linear_term, Real& _constant_term) const
{
	assert(_cuboid_1); assert(_cuboid_2);
	assert(_label_index_1 < num_labels_);
	assert(_label_index_2 < num_labels_);

	// NOTE:
	// Now considering only different label pairs.
	assert(_label_index_1 != _label_index_2);

	const unsigned int num_attributes = MeshCuboidAttributes::k_num_attributes;
	const unsigned int num_features = MeshCuboidFeatures::k_num_features;
	const unsigned int mat_size = _quadratic_term.cols();

	assert(_quadratic_term.rows() == mat_size);
	assert(_linear_term.rows() == mat_size);
	assert(_cuboid_index_1 * num_attributes <= mat_size);
	assert(_cuboid_index_2 * num_attributes <= mat_size);

	_quadratic_term.setZero();
	_linear_term.setZero();
	_constant_term = 0;

	const MeshCuboidJointNormalRelations *relation_12 = relations_[_label_index_1][_label_index_2];
	if (!relation_12) return 0.0;


	MeshCuboidFeatures features_1, features_2;
	Eigen::MatrixXd attributes_to_features_map_1, attributes_to_features_map_2;

	features_1.compute_features(_cuboid_1, &attributes_to_features_map_1);
	features_2.compute_features(_cuboid_2, &attributes_to_features_map_2);

	assert(attributes_to_features_map_1.rows() == num_features);
	assert(attributes_to_features_map_2.rows() == num_features);

	assert(attributes_to_features_map_1.cols() == num_attributes);
	assert(attributes_to_features_map_2.cols() == num_attributes);

	MeshCuboidTransformation transformation_1;
	transformation_1.compute_transformation(_cuboid_1);
	Eigen::MatrixXd rotation_1;
	Eigen::MatrixXd translation_1;
	transformation_1.get_linear_map_transformation(rotation_1, translation_1);

	MeshCuboidTransformation transformation_2;
	transformation_2.compute_transformation(_cuboid_2);
	Eigen::MatrixXd rotation_2;
	Eigen::MatrixXd translation_2;
	transformation_2.get_linear_map_transformation(rotation_2, translation_2);


	// (Ax + b)'C(Ax + b) = x'(A'CA)x + 2*(b'CA)x.
	unsigned int start_index_1 = (_cuboid_index_1 * num_attributes);
	unsigned int start_index_2 = (_cuboid_index_2 * num_attributes);

	Eigen::MatrixXd A1_orig = Eigen::MatrixXd::Zero(2 * num_features, mat_size);
	Eigen::MatrixXd A2_orig = Eigen::MatrixXd::Zero(2 * num_features, mat_size);


	// The rotation is fixed, and attribute-to-feature map is a function of 'cuboid_1' attributes.
	A1_orig.block<num_features, num_attributes>(0, start_index_1)
		= A1_orig.block<num_features, num_attributes>(0, start_index_1)
		+ rotation_1 * attributes_to_features_map_1;

	// The translation is a function of 'cuboid_1' attributes.
	A1_orig.block<num_features, num_attributes>(0, start_index_1)
		= A1_orig.block<num_features, num_attributes>(0, start_index_1)
		+ translation_1;

	// The rotation is fixed, and attribute-to-feature map is a function of 'cuboid_2' attributes.
	A1_orig.block<num_features, num_attributes>(num_features, start_index_2)
		= A1_orig.block<num_features, num_attributes>(num_features, start_index_2)
		+ rotation_1 * attributes_to_features_map_2;

	// The translation is a function of 'cuboid_1' attributes.
	A1_orig.block<num_features, num_attributes>(num_features, start_index_1)
		= A1_orig.block<num_features, num_attributes>(num_features, start_index_1)
		+ translation_1;
	

	// The rotation is fixed, and attribute-to-feature map is a function of 'cuboid_2' attributes.
	A2_orig.block<num_features, num_attributes>(0, start_index_2)
		= A2_orig.block<num_features, num_attributes>(0, start_index_2)
		+ rotation_2 * attributes_to_features_map_2;

	// The translation is a function of 'cuboid_2' attributes.
	A2_orig.block<num_features, num_attributes>(0, start_index_2)
		= A2_orig.block<num_features, num_attributes>(0, start_index_2)
		+ translation_2;

	// The rotation is fixed, and attribute-to-feature map is a function of 'cuboid_1' attributes.
	A2_orig.block<num_features, num_attributes>(num_features, start_index_1)
		= A2_orig.block<num_features, num_attributes>(num_features, start_index_1)
		+ rotation_2 * attributes_to_features_map_1;

	// The translation is a function of 'cuboid_2' attributes.
	A2_orig.block<num_features, num_attributes>(num_features, start_index_2)
		= A2_orig.block<num_features, num_attributes>(num_features, start_index_2)
		+ translation_2;


	// NOTE:
	// Since the center point is always the origin in the local coordinates,
	// it is not used as the feature values.
	const unsigned int num_rows = MeshCuboidJointNormalRelations::k_mat_size;
	Eigen::MatrixXd A(num_rows, mat_size);
	A.topRows(2 * num_features - MeshCuboidFeatures::k_corner_index)
		= A1_orig.bottomRows(2 * num_features - MeshCuboidFeatures::k_corner_index);
	A.bottomRows(2 * num_features - MeshCuboidFeatures::k_corner_index)
		= A2_orig.bottomRows(2 * num_features - MeshCuboidFeatures::k_corner_index);

	Eigen::VectorXd b = -relation_12->get_mean();
	assert(b.rows() == MeshCuboidJointNormalRelations::k_mat_size);

	Eigen::MatrixXd C = relation_12->get_inv_cov();
	assert(C.rows() == MeshCuboidJointNormalRelations::k_mat_size);
	assert(C.cols() == MeshCuboidJointNormalRelations::k_mat_size);


	_quadratic_term = A.transpose() * C * A;
	_linear_term = (b.transpose() * C * A).transpose();
	_constant_term = (b.transpose() * C * b);

	
	MeshCuboidAttributes attributes_1, attributes_2;
	attributes_1.compute_attributes(_cuboid_1);
	attributes_2.compute_attributes(_cuboid_2);

	Eigen::VectorXd x1 = attributes_1.get_attributes();
	Eigen::VectorXd x2 = attributes_2.get_attributes();
	Eigen::VectorXd x = Eigen::VectorXd::Zero(mat_size);
	x.block<num_attributes, 1>(start_index_1, 0) = x1;
	x.block<num_attributes, 1>(start_index_2, 0) = x2;

	Real potential = 0;
	potential += (x.transpose() * _quadratic_term * x);
	potential += (2 * _linear_term.transpose() * x);
	potential += _constant_term;


#ifdef DEBUG_TEST
	Real same_potential = relation_12->compute_error(_cuboid_1, _cuboid_2, &transformation_1, &transformation_2);
	CHECK_NUMERICAL_ERROR(__FUNCTION__, potential, same_potential);
#endif

	//Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ",");
	//std::stringstream filename_sstr;
	//filename_sstr << std::string("quadratic_mat_")
	//	<< _cuboid_index_1 << std::string("_")
	//	<< _cuboid_index_2 << std::string(".csv");

	//std::ofstream csv_file(filename_sstr.str());
	//Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ",");
	//csv_file << _quadratic_term.format(csv_format) << std::endl;
	//csv_file.close();

	return potential;
}

MeshCuboidCondNormalRelationPredictor::MeshCuboidCondNormalRelationPredictor(
	const std::vector< std::vector<MeshCuboidCondNormalRelations *> > &_relations)
	: MeshCuboidPredictor(_relations.size())
	, relations_(_relations)
{
	for (unsigned int label_index = 0; label_index < num_labels_; ++label_index)
		assert(relations_[label_index].size() == num_labels_);
}

void MeshCuboidCondNormalRelationPredictor::get_missing_label_indices(
	const std::list<LabelIndex> &_given_label_indices,
	std::list<LabelIndex> &_missing_label_indices)const
{
	_missing_label_indices.clear();
	unsigned int num_labels = relations_.size();

	bool *is_missing_label = new bool[num_labels];
	memset(is_missing_label, true, num_labels * sizeof(bool));

	for (std::list<LabelIndex>::const_iterator it = _given_label_indices.begin();
		it != _given_label_indices.end(); ++it)
	{
		LabelIndex label_index = (*it);
		assert(label_index < num_labels);
		is_missing_label[label_index] = false;
	}

	for (LabelIndex label_index_1 = 0; label_index_1 < num_labels; ++label_index_1)
	{
		if (is_missing_label[label_index_1])
		{
			bool is_label_confliced = false;

			for (LabelIndex label_index_2 = 0; label_index_2 < num_labels; ++label_index_2)
			{
				if (!is_missing_label[label_index_2] && !relations_[label_index_1][label_index_2])
				{
					is_label_confliced = true;
					break;
				}
			}

			if (!is_label_confliced)
			{
				_missing_label_indices.push_back(label_index_1);
			}
		}
	}

	delete[] is_missing_label;
}

Real MeshCuboidCondNormalRelationPredictor::get_pair_potential(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const MeshCuboidAttributes *_attributes_1, const MeshCuboidAttributes *_attributes_2,
	const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2) const
{
	assert(_cuboid_1); assert(_cuboid_2);
	assert(_attributes_1); assert(_attributes_2);
	assert(_transformation_1); assert(_transformation_2);

	assert(_label_index_1 < num_labels_);
	assert(_label_index_2 < num_labels_);

	// NOTE:
	// Now considering only different label pairs.
	assert(_label_index_1 != _label_index_2);

	Real potential = FLAGS_param_max_potential;

	const MeshCuboidCondNormalRelations *relation_12 = relations_[_label_index_1][_label_index_2];
	const MeshCuboidCondNormalRelations *relation_21 = relations_[_label_index_2][_label_index_1];

	if (relation_12 && relation_21)
	{
		Real potential_12 = relation_12->compute_error(_cuboid_1, _cuboid_2, _transformation_1, _transformation_2);
		Real potential_21 = relation_21->compute_error(_cuboid_2, _cuboid_1, _transformation_2, _transformation_1);
		potential = potential_12 + potential_21;
	}

	return potential;
}

Real MeshCuboidCondNormalRelationPredictor::get_pair_quadratic_form(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const unsigned int _cuboid_index_1, const unsigned int _cuboid_index_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2,
	Eigen::MatrixXd &_quadratic_term, Eigen::VectorXd &_linear_term, Real& _constant_term) const
{
	assert(_cuboid_1); assert(_cuboid_2);
	assert(_label_index_1 < num_labels_);
	assert(_label_index_2 < num_labels_);

	// NOTE:
	// Now considering only different label pairs.
	assert(_label_index_1 != _label_index_2);

	const unsigned int num_attributes = MeshCuboidAttributes::k_num_attributes;
	const unsigned int num_features = MeshCuboidFeatures::k_num_features;
	const unsigned int num_global_coord_points = MeshCuboidFeatures::k_num_local_points;
	const unsigned int num_non_global_coord_features = num_features - 3 * num_global_coord_points;
	const unsigned int mat_size = _quadratic_term.cols();

	assert(_quadratic_term.rows() == mat_size);
	assert(_linear_term.rows() == mat_size);
	assert(_cuboid_index_1 * num_attributes <= mat_size);
	assert(_cuboid_index_2 * num_attributes <= mat_size);

	_quadratic_term.setZero();
	_linear_term.setZero();
	_constant_term = 0;


	MeshCuboidFeatures features_1, features_2;
	Eigen::MatrixXd attributes_to_features_map_1, attributes_to_features_map_2;

	features_1.compute_features(_cuboid_1, &attributes_to_features_map_1);
	features_2.compute_features(_cuboid_2, &attributes_to_features_map_2);

	assert(attributes_to_features_map_1.rows() == num_features);
	assert(attributes_to_features_map_2.rows() == num_features);

	assert(attributes_to_features_map_1.cols() == num_attributes);
	assert(attributes_to_features_map_2.cols() == num_attributes);

	MeshCuboidTransformation transformation_1;
	transformation_1.compute_transformation(_cuboid_1);
	Eigen::MatrixXd rotation_1;
	Eigen::MatrixXd translation_1;
	transformation_1.get_linear_map_transformation(rotation_1, translation_1);

	MeshCuboidTransformation transformation_2;
	transformation_2.compute_transformation(_cuboid_2);


	const MeshCuboidCondNormalRelations *relation_12 = relations_[_label_index_1][_label_index_2];
	if (!relation_12) return 0.0;


	// (Ax + b)'C(Ax + b) = x'(A'CA)x + 2*(b'CA)x.
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_features, mat_size);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(num_features);
	Eigen::MatrixXd C = Eigen::MatrixXd::Zero(num_features, num_features);

	unsigned int start_index_1 = (_cuboid_index_1 * num_attributes);
	unsigned int start_index_2 = (_cuboid_index_2 * num_attributes);

	Eigen::MatrixXd mean_A = Eigen::MatrixXd::Zero(num_features, num_features);
	mean_A.rightCols(num_non_global_coord_features) = relation_12->get_mean_A();

	A.block<num_features, num_attributes>(0, start_index_1)
		= A.block<num_features, num_attributes>(0, start_index_1)
		- mean_A * attributes_to_features_map_1;

	A.block<num_features, num_attributes>(0, start_index_2)
		= A.block<num_features, num_attributes>(0, start_index_2)
		+ rotation_1 * attributes_to_features_map_2;

	A.block<num_features, num_attributes>(0, start_index_1)
		= A.block<num_features, num_attributes>(0, start_index_1)
		+ translation_1;

	b = -relation_12->get_mean_b();

	C = relation_12->get_inv_cov();

	_quadratic_term = A.transpose() * C * A;
	_linear_term = (b.transpose() * C * A).transpose();
	_constant_term = (b.transpose() * C * b);


	MeshCuboidAttributes attributes_1, attributes_2;
	attributes_1.compute_attributes(_cuboid_1);
	attributes_2.compute_attributes(_cuboid_2);

	Eigen::VectorXd x1 = attributes_1.get_attributes();
	Eigen::VectorXd x2 = attributes_2.get_attributes();
	Eigen::VectorXd x = Eigen::VectorXd::Zero(mat_size);
	x.block<num_attributes, 1>(start_index_1, 0) = x1;
	x.block<num_attributes, 1>(start_index_2, 0) = x2;

	Real potential = 0;
	potential += (x.transpose() * _quadratic_term * x);
	potential += (2 * _linear_term.transpose() * x);
	potential += _constant_term;


#ifdef DEBUG_TEST
	Real same_potential = relation_12->compute_error(_cuboid_1, _cuboid_2, &transformation_1, &transformation_2);
	CHECK_NUMERICAL_ERROR(__FUNCTION__, potential, same_potential);
#endif
	
	//std::stringstream filename_sstr;
	//filename_sstr << std::string("quadratic_mat_")
	//	<< _cuboid_index_1 << std::string("_")
	//	<< _cuboid_index_2 << std::string(".csv");
	//std::ofstream csv_file(filename_sstr.str());
	//Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ",");
	//csv_file << quadratic_form.format(csv_format) << std::endl;
	//csv_file.close();
	
	return potential;
}

/*
MeshCuboidPCARelationPredictor::MeshCuboidPCARelationPredictor(
	const std::vector< std::vector<MeshCuboidPCARelations> > &_relations)
	: MeshCuboidPredictor(_relations.size())
	, relations_(_relations)
{
	for (unsigned int label_index = 0; label_index < num_labels_; ++label_index)
		assert(relations_[label_index].size() == num_labels_);
}

Real MeshCuboidPCARelationPredictor::get_pair_potential(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const MeshCuboidAttributes *_attributes_1, const MeshCuboidAttributes *_attributes_2,
	const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2) const
{
	assert(_cuboid_1); assert(_cuboid_2);
	assert(_attributes_1); assert(_attributes_2);
	assert(_transformation_1); assert(_transformation_2);

	assert(_label_index_1 < num_labels_);
	assert(_label_index_2 < num_labels_);

	// NOTE:
	// Now considering only different label pairs.
	assert(_label_index_1 != _label_index_2);

	Real potential = 0.0;

	MeshCuboidFeatures features_1, features_2;
	features_1.compute_features(_cuboid_1);
	features_2.compute_features(_cuboid_2);

	Eigen::VectorXd features_vec_1 = features_1.get_features();
	Eigen::VectorXd features_vec_2 = features_2.get_features();
	Eigen::VectorXd transformed_features_vec_1 = _transformation_2->get_transformed_features(_cuboid_1);
	Eigen::VectorXd transformed_features_vec_2 = _transformation_1->get_transformed_features(_cuboid_2);

	if (_label_index_1 < _label_index_2)
	{
		const MeshCuboidPCARelations &relation = relations_[_label_index_1][_label_index_2];
		potential += relation.compute_error(transformed_features_vec_1, transformed_features_vec_2);
	}
	else
	{
		const MeshCuboidPCARelations &relation = relations_[_label_index_2][_label_index_1];
		potential += relation.compute_error(transformed_features_vec_1, transformed_features_vec_2);
	}

	return potential;
}

Real MeshCuboidPCARelationPredictor::get_pair_quadratic_form(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const unsigned int _cuboid_index_1, const unsigned int _cuboid_index_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2,
	Eigen::MatrixXd &_quadratic_term, Eigen::VectorXd &_linear_term, Real& _constant_term) const
{
	assert(_cuboid_1); assert(_cuboid_2);
	assert(_label_index_1 < num_labels_);
	assert(_label_index_2 < num_labels_);

	// NOTE:
	// Now considering only different label pairs.
	assert(_label_index_1 != _label_index_2);

	const unsigned int num_attributes = MeshCuboidAttributes::k_num_attributes;
	const unsigned int num_features = MeshCuboidFeatures::k_num_features;
	const unsigned int num_global_coord_points = MeshCuboidFeatures::k_num_local_points;
	const unsigned int mat_size = _quadratic_term.cols();

	assert(_quadratic_term.rows() == mat_size);
	assert(_linear_term.rows() == mat_size);
	assert(_cuboid_index_1 * num_attributes <= mat_size);
	assert(_cuboid_index_2 * num_attributes <= mat_size);

	_quadratic_term.setZero();
	_linear_term.setZero();
	_constant_term = 0;


	MeshCuboidFeatures features_1, features_2;
	Eigen::MatrixXd attributes_to_features_map_1, attributes_to_features_map_2;

	features_1.compute_features(_cuboid_1, &attributes_to_features_map_1);
	features_2.compute_features(_cuboid_2, &attributes_to_features_map_2);

	assert(attributes_to_features_map_1.rows() == num_features);
	assert(attributes_to_features_map_2.rows() == num_features);

	assert(attributes_to_features_map_1.cols() == num_attributes);
	assert(attributes_to_features_map_2.cols() == num_attributes);

	MeshCuboidTransformation transformation_1;
	transformation_1.compute_transformation(_cuboid_1);
	Eigen::MatrixXd rotation_1;
	Eigen::MatrixXd translation_1;
	transformation_1.get_feature_transformation(rotation_1, translation_1);

	MeshCuboidTransformation transformation_2;
	transformation_2.compute_transformation(_cuboid_2);
	Eigen::MatrixXd rotation_2;
	Eigen::MatrixXd translation_2;
	transformation_2.get_feature_transformation(rotation_2, translation_2);

	MeshCuboidAttributes attributes_1, attributes_2;
	attributes_1.compute_attributes(_cuboid_1);
	attributes_2.compute_attributes(_cuboid_2);

	// (Ax + b)'C(Ax + b) = x'(A'CA)x + 2*(b'CA)x.
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2 * num_features, mat_size);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(2 * num_features);
	Eigen::MatrixXd C = Eigen::MatrixXd::Zero(2 * num_features, 2 * num_features);

	unsigned int start_index_1 = (_cuboid_index_1 * num_attributes);
	unsigned int start_index_2 = (_cuboid_index_2 * num_attributes);

	Eigen::VectorXd x1 = attributes_1.get_attributes();
	Eigen::VectorXd x2 = attributes_2.get_attributes();
	Eigen::VectorXd x0 = Eigen::VectorXd::Zero(mat_size);
	x0.block<num_attributes, 1>(start_index_1, 0) = x1;
	x0.block<num_attributes, 1>(start_index_2, 0) = x2;

	A.block<num_features, num_attributes>(0, start_index_1)
		= A.block<num_features, num_attributes>(0, start_index_1)
		+ rotation_2 * attributes_to_features_map_1;

	A.block<num_features, num_attributes>(0, start_index_2)
		= A.block<num_features, num_attributes>(0, start_index_2)
		+ translation_2;

	A.block<num_features, num_attributes>(num_features, start_index_2)
		= A.block<num_features, num_attributes>(num_features, start_index_2)
		+ rotation_1 * attributes_to_features_map_2;

	A.block<num_features, num_attributes>(num_features, start_index_1)
		= A.block<num_features, num_attributes>(num_features, start_index_1)
		+ translation_1;


	Eigen::MatrixXd P;
	if (_label_index_1 < _label_index_2)
	{
		const MeshCuboidPCARelations &relation = relations_[_label_index_1][_label_index_2];
		b = b - relation.get_mean();
		P = relation.get_pca_bases();
	}
	else
	{
		const MeshCuboidPCARelations &relation = relations_[_label_index_2][_label_index_1];
		b = b - relation.get_mean();
		P = relation.get_pca_bases();
	}

	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(2 * num_features, 2 * num_features);

	_quadratic_term = _quadratic_term + A.transpose() * (I - P).transpose() * (I - P) * A;
	_linear_term = _linear_term + ((b.transpose() * (I - P).transpose() * (I - P) * A).transpose());
	_constant_term = _constant_term + (b.transpose() * (I - P).transpose() * (I - P) * b);

	_quadratic_term = _quadratic_term + A.transpose() * P * A;
	_linear_term = _linear_term - (A.transpose() * P * A * x0);
	_constant_term = _constant_term + (x0.transpose() * A.transpose() * P * A * x0);


	Real potential = 0;
	potential += (x0.transpose() * _quadratic_term * x0);
	potential += (2 * _linear_term.transpose() * x0);
	potential += _constant_term;


	// DEBUG.
	Eigen::VectorXd transformed_features_vec_1 = transformation_2.get_transformed_features(_cuboid_1);
	Eigen::VectorXd transformed_features_vec_2 = transformation_1.get_transformed_features(_cuboid_2);
	Real same_potential = 0.0;
	if (_label_index_1 < _label_index_2)
	{
		const MeshCuboidPCARelations &relation = relations_[_label_index_1][_label_index_2];
		same_potential = relation.compute_error(transformed_features_vec_1, transformed_features_vec_2);
	}
	else
	{
		const MeshCuboidPCARelations &relation = relations_[_label_index_2][_label_index_1];
		same_potential = relation.compute_error(transformed_features_vec_1, transformed_features_vec_2);
	}

	CHECK_ERROR(__FUNCTION__, potential, same_potential);
	
	// DEBUG.
	//std::stringstream filename_sstr;
	//filename_sstr << std::string("quadratic_mat_")
	//	<< _cuboid_index_1 << std::string("_")
	//	<< _cuboid_index_2 << std::string(".csv");
	//std::ofstream csv_file(filename_sstr.str());
	//Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ",");
	//csv_file << quadratic_form.format(csv_format) << std::endl;
	//csv_file.close();

	return potential;
}

MeshCuboidCCARelationPredictor::MeshCuboidCCARelationPredictor(
	const std::vector< std::vector< std::vector<MeshCuboidCCARelations> > >& _relations)
	: MeshCuboidPredictor(_relations.size())
	, relations_(_relations)
{
	for (unsigned int label_index_1 = 0; label_index_1 < num_labels_; ++label_index_1)
	{
		assert(_relations[label_index_1].size() == num_labels_);
		for (unsigned int label_index_2 = 0; label_index_2 < num_labels_; ++label_index_2)
			assert(_relations[label_index_1][label_index_2].size() == num_labels_);
	}
}

Real MeshCuboidCCARelationPredictor::get_pair_potential(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const MeshCuboidAttributes *_attributes_1, const MeshCuboidAttributes *_attributes_2,
	const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2)const
{
	assert(_cuboid_1); assert(_cuboid_2);
	assert(_attributes_1); assert(_attributes_2);
	assert(_transformation_1); assert(_transformation_2);

	assert(_label_index_1 < num_labels_);
	assert(_label_index_2 < num_labels_);

	// NOTE:
	// Now considering only different label pairs.
	assert(_label_index_1 != _label_index_2);

	Real potential = 0.0;

	// Compute potential.
	//MeshCuboidCCARelations *relation_12;
	//// NOTE:
	//// Use increasing order of label indices for cuboid pairs.
	//if (label_index_1 <= label_index_2)
	//relation_12 = &(_relations[label_index_1][label_index_2]);
	//else
	//relation_12 = &(_relations[label_index_2][label_index_1]);
	//assert(relation_12);
	//potential = relation_12->compute_error(&attributes_1, &attributes_2);
	

	// NOTE:
	const std::vector<MeshCuboidCCARelations> *relation_12;
	// Use increasing order of label indices for cuboid pairs.
	if (_label_index_1 <= _label_index_2)
		relation_12 = &(relations_[_label_index_1][_label_index_2]);
	else
		relation_12 = &(relations_[_label_index_2][_label_index_1]);

	Eigen::VectorXd transformed_features_11 = _transformation_1->get_transformed_features(_cuboid_1);
	Eigen::VectorXd transformed_features_21 = _transformation_1->get_transformed_features(_cuboid_2);
	potential += (*relation_12)[_label_index_1].compute_error(
		transformed_features_11, transformed_features_21);

	Eigen::VectorXd transformed_features_12 = _transformation_2->get_transformed_features(_cuboid_1);
	Eigen::VectorXd transformed_features_22 = _transformation_2->get_transformed_features(_cuboid_2);
	potential += (*relation_12)[_label_index_2].compute_error(
		transformed_features_12, transformed_features_22);

	return potential;
}

MeshCuboidManualRelationPredictor::MeshCuboidManualRelationPredictor(
	const std::vector<MeshCuboidStats> &_single_stats,
	const std::vector< std::vector<MeshCuboidStats> > &_pair_stats)
	: MeshCuboidPredictor(_single_stats.size())
	, single_stats_(_single_stats)
	, pair_stats_(_pair_stats)
{
	assert(_pair_stats.size() == num_labels_);
	for (unsigned int label_index = 0; label_index < num_labels_; ++label_index)
		assert(_pair_stats[label_index].size() == num_labels_);
}

Real MeshCuboidManualRelationPredictor::get_single_potential(
	const MeshCuboid *_cuboid,
	const MeshCuboidAttributes *_attributes,
	const MeshCuboidTransformation *_transformation,
	const LabelIndex _label_index)const
{
	assert(_cuboid);
	assert(_attributes);
	assert(_transformation);

	assert(_label_index < num_labels_);

	Real potential = 0.0;

	// Compute potential.
	//
	MeshSingleCuboidManualFeatures single_features("single_potential");
	single_features.compute_values(_cuboid);
	potential -= single_features.compute_log_probability(single_stats_[_label_index]);

	// NOTE:
	// Pairwise attribute relations are also defined for a single part attributes.
	MeshPairCuboidManualFeatures pair_features("pair_potential");
	pair_features.compute_values(_cuboid, _cuboid);
	potential -= pair_features.compute_log_probability(pair_stats_[_label_index][_label_index]);
	//

	return potential;
}

Real MeshCuboidManualRelationPredictor::get_pair_potential(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const MeshCuboidAttributes *_attributes_1, const MeshCuboidAttributes *_attributes_2,
	const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
	const LabelIndex _label_index_1, const LabelIndex _label_index_2)const
{
	assert(_cuboid_1); assert(_cuboid_2);
	assert(_attributes_1); assert(_attributes_2);
	assert(_transformation_1); assert(_transformation_2);

	assert(_label_index_1 < num_labels_);
	assert(_label_index_2 < num_labels_);

	Real potential = 0.0;

	// Compute potential.
	//
	MeshPairCuboidManualFeatures pair_features("pair_potential");
	pair_features.compute_values(_cuboid_1, _cuboid_2);

	// NOTE:
	// Check symmetry:
	// pair_features.compute_values(label_2, cuboid_2, label_1, cuboid_1);

	// NOTE:
	// Use increasing order of label indices for cuboid pairs.
	if (_label_index_1 <= _label_index_2)
	{
		potential -= pair_features.compute_log_probability(
			pair_stats_[_label_index_1][_label_index_2]);
	}
	else
	{
		potential -= pair_features.compute_log_probability(
			pair_stats_[_label_index_2][_label_index_1]);
	}

	return potential;
}
*/
