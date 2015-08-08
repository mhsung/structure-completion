#include "MeshCuboidRelation.h"

#include "Utilities.h"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <Eigen/Eigenvalues> 
#include <Eigen/LU> 


// Note:
// Assume that z = 0 plane is the ground plane, and +z-axis is the normal direction.
const MyMesh::Normal MeshCuboidAttributes::k_up_direction = MyMesh::Normal(0.0, 0.0, 1.0);

MeshCuboidAttributes::MeshCuboidAttributes()
	: object_name_("")
{
	initialize();
}

MeshCuboidAttributes::MeshCuboidAttributes(const std::string _object_name)
	: object_name_(_object_name)
{
	initialize();
}

MeshCuboidAttributes::MeshCuboidAttributes(const MeshCuboidAttributes& _other)
	: object_name_(_other.object_name_)
	, attributes_(_other.attributes_)
{
}

MeshCuboidAttributes::~MeshCuboidAttributes()
{
}

void MeshCuboidAttributes::initialize()
{
	// Note:
	// Initialized with NaN.
	if (!std::numeric_limits<Real>::has_quiet_NaN)
	{
		std::cerr << "Error: " << typeid(Real).name() << " does not support quiet NaN." << std::endl;
		do {
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	const int num_attributes = static_cast<int>(k_num_attributes);

	attributes_ = Eigen::VectorXd::Zero(num_attributes);
	for (int attribute_index = 0; attribute_index < num_attributes; ++attribute_index)
	{
		attributes_[attribute_index] = std::numeric_limits<Real>::quiet_NaN();
	}
}

void MeshCuboidAttributes::compute_attributes(const MeshCuboid *_cuboid)
{
	assert(attributes_.size() == static_cast<int>(k_num_attributes));

	for (unsigned int i = 0; i < 3; ++i)
	{
		//attributes_[k_center_index + i] = _cuboid->get_bbox_center()[i];

		for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
			attributes_[k_corner_index + 3 * corner_index + i] = _cuboid->get_bbox_corner(corner_index)[i];
	}
}

void MeshCuboidAttributes::get_attribute_collection_matrix(const std::list<MeshCuboidAttributes *>& _stats,
	Eigen::MatrixXd& _values)
{
	unsigned int num_objects = _stats.size();

	const int num_attributes = static_cast<int>(k_num_attributes);
	_values = Eigen::MatrixXd(num_objects, num_attributes);

	unsigned int object_index = 0;
	for (std::list<MeshCuboidAttributes *>::const_iterator stat_it = _stats.begin(); stat_it != _stats.end();
		++stat_it, ++object_index)
	{
		assert(*stat_it);

		for (int attribute_index = 0; attribute_index < num_attributes; ++attribute_index)
			_values(object_index, attribute_index) = (*stat_it)->attributes_[attribute_index];
	}
}

bool MeshCuboidAttributes::save_attribute_collection(const std::list<MeshCuboidAttributes *>& _stats,
	const char* _filename)
{
	const int num_attributes = static_cast<int>(k_num_attributes);

	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}
	std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	std::cout << std::scientific;

	//unsigned int num_objects = _stats.size();

	// Write attribute values.
	unsigned int object_index = 0;
	for (std::list<MeshCuboidAttributes *>::const_iterator stat_it = _stats.begin(); stat_it != _stats.end();
		++stat_it, ++object_index)
	{
		assert(*stat_it);

		for (int attribute_index = 0; attribute_index < num_attributes; ++attribute_index)
		{
			Real value = (*stat_it)->attributes_[attribute_index];
			if (std::isnan(value))
			{
				// Note:
				// Undefined attributes are recorded as NaN.
				file << "NaN,";
			}
			else
			{
				file << value << ",";
			}
		}
		file << std::endl;
	}

	file.close();
	return true;
}


MeshCuboidFeatures::MeshCuboidFeatures()
	: object_name_("")
{
	initialize();
}

MeshCuboidFeatures::MeshCuboidFeatures(const std::string _object_name)
	: object_name_(_object_name)
{
	initialize();
}

MeshCuboidFeatures::MeshCuboidFeatures(const MeshCuboidFeatures& _other)
	: object_name_(_other.object_name_)
	, features_(_other.features_)
{
}

MeshCuboidFeatures::~MeshCuboidFeatures()
{
}

void MeshCuboidFeatures::initialize()
{
	// Note:
	// Initialized with NaN.
	if (!std::numeric_limits<Real>::has_quiet_NaN)
	{
		std::cerr << "Error: " << typeid(Real).name() << " does not support quiet NaN." << std::endl;
		do {
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	features_ = Eigen::VectorXd::Zero(MeshCuboidFeatures::k_num_features);
	for (int attribute_index = 0; attribute_index < MeshCuboidFeatures::k_num_features; ++attribute_index)
	{
		features_[attribute_index] = std::numeric_limits<Real>::quiet_NaN();
	}
}

void MeshCuboidFeatures::compute_features(const MeshCuboid *_cuboid,
	Eigen::MatrixXd *_attributes_to_features_map)
{
	const unsigned int num_attributes = MeshCuboidAttributes::k_num_attributes;
	assert(features_.size() == static_cast<int>(MeshCuboidFeatures::k_num_features));


	// Linear map.
	Eigen::MatrixXd attributes_to_features_map = 
		Eigen::MatrixXd::Zero(MeshCuboidFeatures::k_num_features, num_attributes);

	Eigen::RowVector3d up_direction_vec;
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		up_direction_vec(axis_index) = MeshCuboidAttributes::k_up_direction[axis_index];

	unsigned int next_feature_index = 0;


	// Center point.
	for (unsigned int i = 0; i < 3; ++i)
	{
		features_[next_feature_index] = _cuboid->get_bbox_center()[i];

		//attributes_to_features_map.row(next_feature_index)(
		//	MeshCuboidAttributes::k_center_index + i) = 1.0;
		for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
		{
			attributes_to_features_map.row(next_feature_index)(
				MeshCuboidAttributes::k_corner_index + 3 * corner_index + i) =
				1.0 / MeshCuboid::k_num_corners;
		}

		++next_feature_index;
	}


	// Corner points.
	for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
	{
		MyMesh::Point corner = _cuboid->get_bbox_corner(corner_index);

		for (unsigned int i = 0; i < 3; ++i)
		{
			features_[next_feature_index] = corner[i];
			attributes_to_features_map.row(next_feature_index)(
				MeshCuboidAttributes::k_corner_index + 3 * corner_index + i) = 1.0;
			++next_feature_index;
		}
	}


	// Center height.
	features_[next_feature_index] = dot(MeshCuboidAttributes::k_up_direction, _cuboid->get_bbox_center());
	for (unsigned int i = 0; i < 3; ++i)
	{
		//attributes_to_features_map.row(next_feature_index)(
		//	MeshCuboidAttributes::k_center_index + i) = MeshCuboidAttributes::k_up_direction[i];
		for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
		{
			attributes_to_features_map.row(next_feature_index)(
				MeshCuboidAttributes::k_corner_index + 3 * corner_index + i) =
				(1.0 / MeshCuboid::k_num_corners) * MeshCuboidAttributes::k_up_direction[i];
		}
	}
	++next_feature_index;


	// Corner height.
	for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
	{
		MyMesh::Point corner = _cuboid->get_bbox_corner(corner_index);
		features_[next_feature_index] = dot(MeshCuboidAttributes::k_up_direction, corner);

		for (unsigned int i = 0; i < 3; ++i)
		{
			attributes_to_features_map.row(next_feature_index)(
				MeshCuboidAttributes::k_corner_index + 3 * corner_index + i) =
				MeshCuboidAttributes::k_up_direction[i];
		}
		++next_feature_index;
	}

	assert(next_feature_index == static_cast<int>(MeshCuboidFeatures::k_num_features));


#ifdef DEBUG_TEST
	MeshCuboidAttributes attributes;
	attributes.compute_attributes(_cuboid);
	Eigen::VectorXd same_features = attributes_to_features_map * attributes.get_attributes();

	assert(same_features.rows() == MeshCuboidFeatures::k_num_features);
	Real error = (same_features - features_).array().abs().sum();

	CHECK_NUMERICAL_ERROR(__FUNCTION__, error);
#endif

	// Optional.
	if (_attributes_to_features_map)
		(*_attributes_to_features_map) = attributes_to_features_map;
}

void MeshCuboidFeatures::get_feature_collection_matrix(const std::list<MeshCuboidFeatures *>& _stats,
	Eigen::MatrixXd& _values)
{
	unsigned int num_objects = _stats.size();

	const int num_features = static_cast<int>(k_num_features);
	_values = Eigen::MatrixXd(num_objects, num_features);

	unsigned int object_index = 0;
	for (std::list<MeshCuboidFeatures *>::const_iterator stat_it = _stats.begin(); stat_it != _stats.end();
		++stat_it, ++object_index)
	{
		assert(*stat_it);

		for (int feature_index = 0; feature_index < num_features; ++feature_index)
			_values(object_index, feature_index) = (*stat_it)->features_[feature_index];
	}
}

bool MeshCuboidFeatures::load_feature_collection(const char* _filename,
	std::list<MeshCuboidFeatures *>& _stats)
{
	for (std::list<MeshCuboidFeatures *>::iterator it = _stats.begin(); it != _stats.end(); ++it)
		delete (*it);
	_stats.clear();

	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't load file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	std::string buffer;
	std::stringstream strstr;
	std::string token;
	bool succeded = true;

	while (!file.eof())
	{
		std::getline(file, buffer);
		if (buffer == "") break;

		strstr.str(std::string());
		strstr.clear();
		strstr.str(buffer);

		MeshCuboidFeatures *new_features = new MeshCuboidFeatures();
		assert(new_features);

		for (int feature_index = 0; feature_index < MeshCuboidFeatures::k_num_features; ++feature_index)
		{
			std::getline(strstr, token, ',');
			if (strstr.eof())
			{
				std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
				succeded = false;
				break;
			}

			if (token == "NaN")
			{
				// Note:
				// Undefined attributes are recorded as NaN.
				new_features->features_[feature_index] = std::numeric_limits<Real>::quiet_NaN();
			}
			else
			{
				new_features->features_[feature_index] = std::stof(token);
			}
		}

		if (!succeded) break;
		_stats.push_back(new_features);
	}

	if (!succeded)
	{
		for (std::list<MeshCuboidFeatures *>::iterator it = _stats.begin(); it != _stats.end(); ++it)
			delete (*it);
		_stats.clear();
		return false;
	}

	return true;
}

bool MeshCuboidFeatures::save_feature_collection(const char* _filename,
	const std::list<MeshCuboidFeatures *>& _stats)
{
	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}
	std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	std::cout << std::scientific;

	//unsigned int num_objects = _stats.size();

	// Write attribute values.
	unsigned int object_index = 0;
	for (std::list<MeshCuboidFeatures *>::const_iterator stat_it = _stats.begin(); stat_it != _stats.end();
		++stat_it, ++object_index)
	{
		assert(*stat_it);

		for (int feature_index = 0; feature_index < MeshCuboidFeatures::k_num_features; ++feature_index)
		{
			Real value = (*stat_it)->features_[feature_index];
			if (std::isnan(value))
			{
				// Note:
				// Undefined attributes are recorded as NaN.
				file << "NaN,";
			}
			else
			{
				file << value << ",";
			}
		}
		file << std::endl;
	}

	file.close();
	return true;
}

MeshCuboidTransformation::MeshCuboidTransformation()
: object_name_("")
{
	initialize();
}

MeshCuboidTransformation::MeshCuboidTransformation(const std::string _object_name)
: object_name_(_object_name)
{
	initialize();
}

MeshCuboidTransformation::~MeshCuboidTransformation()
{
}

void MeshCuboidTransformation::initialize()
{
	first_translation_ = Eigen::Vector3d::Zero();
	second_rotation_ = Eigen::Matrix3d::Identity();
}

void MeshCuboidTransformation::compute_transformation(const MeshCuboid *_cuboid)
{
	assert(_cuboid);

	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		for (unsigned int i = 0; i < 3; ++i)
		{
			second_rotation_.row(axis_index)(i) = _cuboid->get_bbox_axis(axis_index)[i];
		}
		first_translation_(axis_index) = -_cuboid->get_bbox_center()[axis_index];
	};
}

Eigen::VectorXd
MeshCuboidTransformation::get_transformed_features(
const MeshCuboidFeatures& _other_features)const
{
	Eigen::VectorXd transformed_features = _other_features.get_features();

	for (unsigned int i = 0; i < MeshCuboidFeatures::k_num_local_points; ++i)
	{
		Eigen::VectorXd sub_values = transformed_features.block(3 * i, 0, 3, 1);
		sub_values = first_translation_ + sub_values;
		sub_values = second_rotation_ * sub_values;
		transformed_features.block(3 * i, 0, 3, 1) = sub_values;
	}

//#ifdef DEBUG_TEST
//	Eigen::MatrixXd rotation;
//	Eigen::MatrixXd translation;
//	get_linear_map_transformation(rotation, translation);
//	Eigen::VectorXd same_transformed_features = rotation * (translation * _other_features.get_features());
//
//	assert(same_transformed_features.rows() == MeshCuboidFeatures::k_num_features);
//	Real error = (same_transformed_features - transformed_features).array().abs().sum();
//
//	CHECK_NUMERICAL_ERROR(__FUNCTION__, error);
//#endif

	return transformed_features;
}

Eigen::VectorXd
MeshCuboidTransformation::get_transformed_features(const MeshCuboid *_other_cuboid)const
{
	assert(_other_cuboid);

	MeshCuboidFeatures other_features;
	other_features.compute_features(_other_cuboid);

	return get_transformed_features(other_features);
}

Eigen::VectorXd
MeshCuboidTransformation::get_inverse_transformed_features(
	const MeshCuboidFeatures& _other_features)const
{
	Eigen::VectorXd inverse_transformed_features = _other_features.get_features();

	for (unsigned int i = 0; i < MeshCuboidFeatures::k_num_local_points; ++i)
	{
		Eigen::VectorXd sub_values = inverse_transformed_features.block(3 * i, 0, 3, 1);
		sub_values = second_rotation_.transpose() * sub_values;
		sub_values = -first_translation_ + sub_values;
		inverse_transformed_features.block(3 * i, 0, 3, 1) = sub_values;
	}

//#ifdef DEBUG_TEST
//	Eigen::MatrixXd rotation;
//	Eigen::MatrixXd translation;
//	get_linear_map_inverse_transformation(rotation, translation);
//	Eigen::VectorXd same_transformed_features = rotation * (translation * _other_features.get_features());
//
//	assert(same_transformed_features.rows() == MeshCuboidFeatures::k_num_features);
//	Real error = (same_transformed_features - inverse_transformed_features).array().abs().sum();
//
//	CHECK_NUMERICAL_ERROR(__FUNCTION__, error);
//#endif

	return inverse_transformed_features;
}

Eigen::VectorXd
MeshCuboidTransformation::get_inverse_transformed_features(const MeshCuboid *_other_cuboid)const
{
	assert(_other_cuboid);

	MeshCuboidFeatures other_features;
	other_features.compute_features(_other_cuboid);

	return get_inverse_transformed_features(other_features);
}

void MeshCuboidTransformation::get_transformation(
	Eigen::Matrix3d &_rotation, Eigen::Vector3d &_translation) const
{
	_rotation = second_rotation_;
	_translation = first_translation_;
	_translation = _rotation * _translation;
}

void MeshCuboidTransformation::get_inverse_transformation(
	Eigen::Matrix3d &_rotation, Eigen::Vector3d &_translation) const
{
	_rotation = second_rotation_.transpose();
	_translation = -first_translation_;
}

void MeshCuboidTransformation::get_linear_map_transformation(
	Eigen::MatrixXd &_rotation, Eigen::MatrixXd &_translation) const
{
	const unsigned int num_features = MeshCuboidFeatures::k_num_features;
	_rotation = Eigen::MatrixXd::Identity(num_features, num_features);
	_translation = Eigen::MatrixXd::Zero(num_features, MeshCuboidAttributes::k_num_attributes);

	for (unsigned int k = 0; k < MeshCuboidFeatures::k_num_local_points; ++k)
	{
		_rotation.block<3, 3>(3 * k, 3 * k) = second_rotation_;

		for (unsigned int i = 0; i < 3; ++i)
		{
			//_translation(3 * k + i, MeshCuboidAttributes::k_center_index + i) = -1;
			for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
			{
				_translation(3 * k + i, MeshCuboidAttributes::k_corner_index + 3 * corner_index + i) =
					-1.0 / MeshCuboid::k_num_corners;
			}
		}
	}

	_translation = _rotation * _translation;
}

void MeshCuboidTransformation::get_linear_map_inverse_transformation(
	Eigen::MatrixXd &_rotation, Eigen::MatrixXd &_translation) const
{
	const unsigned int num_features = MeshCuboidFeatures::k_num_features;
	_rotation = Eigen::MatrixXd::Identity(num_features, num_features);
	_translation = Eigen::MatrixXd::Zero(num_features, MeshCuboidAttributes::k_num_attributes);

	for (unsigned int k = 0; k < MeshCuboidFeatures::k_num_local_points; ++k)
	{
		_rotation.block<3, 3>(3 * k, 3 * k) = second_rotation_.transpose();

		for (unsigned int i = 0; i < 3; ++i)
		{
			//_translation(3 * k + i, MeshCuboidAttributes::k_center_index + i) = 1;
			for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
			{
				_translation(3 * k + i, MeshCuboidAttributes::k_corner_index + 3 * corner_index + i) =
					1.0 / MeshCuboid::k_num_corners;
			}
		}
	}
}

bool MeshCuboidTransformation::load_transformation_collection(const char* _filename,
	std::list<MeshCuboidTransformation *>& _stats)
{
	for (std::list<MeshCuboidTransformation *>::iterator it = _stats.begin(); it != _stats.end(); ++it)
		delete (*it);
	_stats.clear();

	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't load file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	std::string buffer;
	std::stringstream strstr;
	std::string token;
	bool succeded = true;

	while (!file.eof())
	{
		std::getline(file, buffer);
		if (buffer == "") break;

		strstr.str(std::string());
		strstr.clear();
		strstr.str(buffer);

		MeshCuboidTransformation *new_transformation = new MeshCuboidTransformation();
		assert(new_transformation);

		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		{
			for (unsigned int i = 0; i < 3; ++i)
			{
				std::getline(strstr, token, ',');
				if (strstr.eof())
				{
					std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
					succeded = false;
					break;
				}

				if (!succeded) break;
				new_transformation->second_rotation_.col(axis_index)(i) = std::stof(token);
			}
		}

		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		{
			std::getline(strstr, token, ',');
			if (strstr.eof())
			{
				std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
				succeded = false;
				break;
			}

			new_transformation->first_translation_(axis_index) = std::stof(token);
		}

		if (!succeded) break;
		_stats.push_back(new_transformation);
	}

	if (!succeded)
	{
		for (std::list<MeshCuboidTransformation *>::iterator it = _stats.begin(); it != _stats.end(); ++it)
			delete (*it);
		_stats.clear();
		return false;
	}

	return true;
}

bool MeshCuboidTransformation::save_transformation_collection(const char* _filename,
	const std::list<MeshCuboidTransformation *>& _stats)
{
	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}
	std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	std::cout << std::scientific;

	//unsigned int num_objects = _stats.size();

	// Write attribute values.
	unsigned int object_index = 0;
	for (std::list<MeshCuboidTransformation *>::const_iterator stat_it = _stats.begin(); stat_it != _stats.end();
		++stat_it, ++object_index)
	{
		assert(*stat_it);

		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
			for (unsigned int i = 0; i < 3; ++i)
				file << (*stat_it)->second_rotation_.col(axis_index)(i) << ",";

		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
			file << (*stat_it)->first_translation_(axis_index) << ",";

		file << std::endl;
	}

	file.close();
	return true;
}

MeshCuboidJointNormalRelations::MeshCuboidJointNormalRelations()
{
	mean_ = Eigen::VectorXd::Zero(k_mat_size);
	inv_cov_ = Eigen::MatrixXd::Zero(k_mat_size, k_mat_size);
}

MeshCuboidJointNormalRelations::~MeshCuboidJointNormalRelations()
{

}

bool MeshCuboidJointNormalRelations::load_joint_normal_csv(const char* _filename)
{
	Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ", ", "", "", "", "", "");

	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	std::string buffer;
	std::stringstream strstr;
	std::string token;

	std::getline(file, buffer);
	strstr.str(std::string());
	strstr.clear();
	strstr.str(buffer);

	for (int j = 0; j < k_mat_size; ++j)
	{
		if (strstr.eof())
		{
			std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}

		std::getline(strstr, token, ',');

		// Transposed.
		mean_(j) = std::stof(token);
	}

	for (int i = 0; i < k_mat_size; ++i)
	{
		if (file.eof())
		{
			std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}

		std::getline(file, buffer);
		strstr.str(std::string());
		strstr.clear();
		strstr.str(buffer);

		for (int j = 0; j < k_mat_size; ++j)
		{
			if (strstr.eof())
			{
				std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
				return false;
			}

			std::getline(strstr, token, ',');
			inv_cov_(i, j) = std::stof(token);
		}
	}

	file.close();

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(inv_cov_);
	Real min_eigenvalue = es.eigenvalues().minCoeff();
	if (min_eigenvalue < -1.0E-6)
	{
		std::cerr << "Error: The inverse covariance matrix is not positive-semidefinite: ("
			<< _filename << " = " << min_eigenvalue << ")" << std::endl;
		do {
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	Real symmetry_diff = (inv_cov_ - inv_cov_.transpose()).array().abs().sum();
	if (symmetry_diff > 1.0E-6)
	{
		std::cerr << "Error: The inverse covariance matrix is not symmetric: ("
			<< _filename << " = " << symmetry_diff << ")" << std::endl;
		do {
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	return true;
}

bool MeshCuboidJointNormalRelations::save_joint_normal_csv(const char* _filename)const
{
	Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ",");

	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}
	std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	std::cout << std::scientific;

	file << mean_.transpose().format(csv_format) << std::endl;
	file << inv_cov_.format(csv_format) << std::endl;

	file.close();
	return true;
}

bool MeshCuboidJointNormalRelations::load_joint_normal_dat(const char* _filename)
{
	std::ifstream file(_filename, std::ios::in | std::ios::binary);
	if (!file.good())
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	int16_t rows, cols;

	rows = 0, cols = 0;
	file.read((char*)(&rows), sizeof(int16_t));
	file.read((char*)(&cols), sizeof(int16_t));
	assert(static_cast<int>(rows) == k_mat_size);
	assert(static_cast<int>(cols) == 1);
	mean_.resize(rows);
	file.read((char *)mean_.data(), rows*cols*sizeof(Eigen::MatrixXd::Scalar));
	
	rows = 0, cols = 0;
	file.read((char*)(&rows), sizeof(int16_t));
	file.read((char*)(&cols), sizeof(int16_t));
	assert(static_cast<int>(rows) == k_mat_size);
	assert(static_cast<int>(cols) == k_mat_size);
	inv_cov_.resize(rows, cols);
	file.read((char *)inv_cov_.data(), rows*cols*sizeof(Eigen::MatrixXd::Scalar));

	file.close();

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(inv_cov_);
	Real min_eigenvalue = es.eigenvalues().minCoeff();
	if (min_eigenvalue < -1.0E-6)
	{
		std::cerr << "Error: The inverse covariance matrix is not positive-semidefinite: ("
			<< _filename << " = " << min_eigenvalue << ")" << std::endl;
		do {
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	Real symmetry_diff = (inv_cov_ - inv_cov_.transpose()).array().abs().sum();
	if (symmetry_diff > 1.0E-6)
	{
		std::cerr << "Error: The inverse covariance matrix is not symmetric: ("
			<< _filename << " = " << symmetry_diff << ")" << std::endl;
		do {
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	return true;
}

void MeshCuboidJointNormalRelations::get_pairwise_cuboid_features(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
	Eigen::VectorXd &_pairwise_features_vec)
{
	assert(_cuboid_1);
	assert(_cuboid_2);

	MeshCuboidFeatures features_1;
	features_1.compute_features(_cuboid_1);

	MeshCuboidFeatures features_2;
	features_2.compute_features(_cuboid_2);

	get_pairwise_cuboid_features(features_1, features_2,
		_transformation_1, _transformation_2, _pairwise_features_vec);
}

void MeshCuboidJointNormalRelations::get_pairwise_cuboid_features(
	const MeshCuboidFeatures &_features_1, const MeshCuboidFeatures &_features_2,
	const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
	Eigen::VectorXd &_pairwise_features_vec)
{
	assert(_transformation_1);
	assert(_transformation_2);

	Eigen::VectorXd transformed_features_vec_11 = _transformation_1->get_transformed_features(_features_1);
	assert(std::abs(transformed_features_vec_11[0]) < NUMERIAL_ERROR_THRESHOLD);
	assert(std::abs(transformed_features_vec_11[1]) < NUMERIAL_ERROR_THRESHOLD);
	assert(std::abs(transformed_features_vec_11[2]) < NUMERIAL_ERROR_THRESHOLD);
	assert(transformed_features_vec_11.rows() == MeshCuboidFeatures::k_num_features);

	Eigen::VectorXd transformed_features_vec_12 = _transformation_1->get_transformed_features(_features_2);
	assert(transformed_features_vec_12.rows() == MeshCuboidFeatures::k_num_features);


	Eigen::VectorXd transformed_features_vec_22 = _transformation_2->get_transformed_features(_features_2);
	assert(std::abs(transformed_features_vec_22[0]) < NUMERIAL_ERROR_THRESHOLD);
	assert(std::abs(transformed_features_vec_22[1]) < NUMERIAL_ERROR_THRESHOLD);
	assert(std::abs(transformed_features_vec_22[2]) < NUMERIAL_ERROR_THRESHOLD);
	assert(transformed_features_vec_22.rows() == MeshCuboidFeatures::k_num_features);

	Eigen::VectorXd transformed_features_vec_21 = _transformation_2->get_transformed_features(_features_1);
	assert(transformed_features_vec_21.rows() == MeshCuboidFeatures::k_num_features);


	// NOTE:
	// Since the center point is always the origin in the local coordinates,
	// it is not used as the feature values.
	assert(k_mat_size == 2 * (2 * MeshCuboidFeatures::k_num_features - MeshCuboidFeatures::k_corner_index));
	_pairwise_features_vec.resize(k_mat_size);
	_pairwise_features_vec <<
		transformed_features_vec_11.bottomRows(MeshCuboidFeatures::k_num_features - MeshCuboidFeatures::k_corner_index),
		transformed_features_vec_12,
		transformed_features_vec_22.bottomRows(MeshCuboidFeatures::k_num_features - MeshCuboidFeatures::k_corner_index),
		transformed_features_vec_21;
}

double MeshCuboidJointNormalRelations::compute_error(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2) const
{
	Eigen::VectorXd pairwise_cuboid_feature;
	get_pairwise_cuboid_features(_cuboid_1, _cuboid_2, _transformation_1, _transformation_2,
		pairwise_cuboid_feature);

	assert(mean_.rows() == pairwise_cuboid_feature.rows());
	assert(inv_cov_.rows() == pairwise_cuboid_feature.rows());
	assert(inv_cov_.cols() == pairwise_cuboid_feature.rows());

	Eigen::VectorXd diff = pairwise_cuboid_feature - mean_;

	// Mahalanobis norm.
	double error = diff.transpose() * inv_cov_ * diff;
	assert(error >= 0);

	//std::cerr << "Negative error value (error = " << error << ")" << std::endl;
	//Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ",");

	//Real inv_cov_det = inv_cov_.determinant();
	//std::cout << "inv_cov_det = " << inv_cov_det << std::endl;

	//Real symmetry_diff = (inv_cov_ - inv_cov_.transpose()).array().abs().sum();
	//std::cerr << "symmetry_diff = " << symmetry_diff << std::endl;

	//std::cout << "Norm(diff) = " << diff.transpose() * diff << std::endl;

	//std::cout << "diff = " << std::endl;
	//std::cout << diff.format(csv_format) << std::endl;

	//do {
	//	std::cout << '\n' << "Press the Enter key to continue.";
	//} while (std::cin.get() != '\n');

	return error;
}

double MeshCuboidJointNormalRelations::compute_conditional_error(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const MeshCuboidTransformation *_transformation_1) const
{
	Eigen::VectorXd conditional_pairwise_cuboid_feature;
	assert(_transformation_1);

	MeshCuboidFeatures features_1;
	features_1.compute_features(_cuboid_1);

	MeshCuboidFeatures features_2;
	features_2.compute_features(_cuboid_2);

	Eigen::VectorXd transformed_features_vec_11 = _transformation_1->get_transformed_features(features_1);
	assert(std::abs(transformed_features_vec_11[0]) < NUMERIAL_ERROR_THRESHOLD);
	assert(std::abs(transformed_features_vec_11[1]) < NUMERIAL_ERROR_THRESHOLD);
	assert(std::abs(transformed_features_vec_11[2]) < NUMERIAL_ERROR_THRESHOLD);
	assert(transformed_features_vec_11.rows() == MeshCuboidFeatures::k_num_features);

	Eigen::VectorXd transformed_features_vec_12 = _transformation_1->get_transformed_features(features_2);
	assert(transformed_features_vec_12.rows() == MeshCuboidFeatures::k_num_features);

	// NOTE:
	// Since the center point is always the origin in the local coordinates,
	// it is not used as the feature values.
	const int num_rows = 2 * MeshCuboidFeatures::k_num_features - MeshCuboidFeatures::k_corner_index;
	conditional_pairwise_cuboid_feature.resize(num_rows);
	conditional_pairwise_cuboid_feature <<
		transformed_features_vec_11.bottomRows(MeshCuboidFeatures::k_num_features - MeshCuboidFeatures::k_corner_index),
		transformed_features_vec_12;

	Eigen::VectorXd conditional_mean = mean_.segment(0, num_rows);
	Eigen::MatrixXd conditional_inv_cov = inv_cov_.block(0, 0, num_rows, num_rows);
	Eigen::VectorXd diff = conditional_pairwise_cuboid_feature - conditional_mean;

	// Mahalanobis norm.
	double error = diff.transpose() * conditional_inv_cov * diff;
	assert(error >= 0);

	return error;
}

MeshCuboidCondNormalRelations::MeshCuboidCondNormalRelations()
{
	assert(MeshCuboidFeatures::k_num_global_feature_values > 0);

	mean_A_ = Eigen::MatrixXd::Zero(MeshCuboidFeatures::k_num_features, MeshCuboidFeatures::k_num_global_feature_values);
	mean_b_ = Eigen::VectorXd::Zero(MeshCuboidFeatures::k_num_features);
	inv_cov_ = Eigen::MatrixXd::Zero(MeshCuboidFeatures::k_num_features, MeshCuboidFeatures::k_num_features);
}

MeshCuboidCondNormalRelations::~MeshCuboidCondNormalRelations()
{

}

bool MeshCuboidCondNormalRelations::load_cond_normal_csv(const char* _filename)
{
	Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ", ", "", "", "", "", "");

	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	std::string buffer;
	std::stringstream strstr;
	std::string token;

	for (int i = 0; i < MeshCuboidFeatures::k_num_global_feature_values; ++i)
	{
		if (file.eof())
		{
			std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}

		std::getline(file, buffer);
		strstr.str(std::string());
		strstr.clear();
		strstr.str(buffer);

		for (int j = 0; j < MeshCuboidFeatures::k_num_features; ++j)
		{
			if (strstr.eof())
			{
				std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
				return false;
			}

			std::getline(strstr, token, ',');

			// Transposed.
			mean_A_(j, i) = std::stof(token);
		}
	}

	std::getline(file, buffer);
	strstr.str(std::string());
	strstr.clear();
	strstr.str(buffer);

	for (int j = 0; j < MeshCuboidFeatures::k_num_features; ++j)
	{
		if (strstr.eof())
		{
			std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}

		std::getline(strstr, token, ',');

		// Transposed.
		mean_b_(j) = std::stof(token);
	}

	for (int i = 0; i < MeshCuboidFeatures::k_num_features; ++i)
	{
		if (file.eof())
		{
			std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}

		std::getline(file, buffer);
		strstr.str(std::string());
		strstr.clear();
		strstr.str(buffer);

		for (int j = 0; j < MeshCuboidFeatures::k_num_features; ++j)
		{
			if (strstr.eof())
			{
				std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
				return false;
			}

			std::getline(strstr, token, ',');
			inv_cov_(i, j) = std::stof(token);
		}
	}

	file.close();

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(inv_cov_);
	Real min_eigenvalue = es.eigenvalues().minCoeff();
	if (min_eigenvalue < -1.0E-6)
	{
		std::cerr << "Error: The inverse covariance matrix is not positive-semidefinite: ("
			<< _filename << " = " << min_eigenvalue << ")" << std::endl;
		do {
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	Real symmetry_diff = (inv_cov_ - inv_cov_.transpose()).array().abs().sum();
	if (symmetry_diff > 1.0E-6)
	{
		std::cerr << "Error: The inverse covariance matrix is not symmetric: ("
			<< _filename << " = " << symmetry_diff << ")" << std::endl;
		do {
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	return true;
}

bool MeshCuboidCondNormalRelations::save_cond_normal_csv(const char* _filename)const
{
	Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ",");

	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}
	std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	std::cout << std::scientific;

	// NOTE:
	// mean_A_. mean_b_: Transposed.
	file << mean_A_.transpose().format(csv_format) << std::endl;
	file << mean_b_.transpose().format(csv_format) << std::endl;
	file << inv_cov_.format(csv_format) << std::endl;

	file.close();
	return true;
}

bool MeshCuboidCondNormalRelations::load_cond_normal_dat(const char* _filename)
{
	std::ifstream file(_filename, std::ios::in | std::ios::binary);
	if (!file.good())
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	int16_t rows, cols;

	rows = 0, cols = 0;
	file.read((char*)(&rows), sizeof(int16_t));
	file.read((char*)(&cols), sizeof(int16_t));
	assert(static_cast<int>(rows) == MeshCuboidFeatures::k_num_features);
	assert(static_cast<int>(cols) == MeshCuboidFeatures::k_num_global_feature_values);
	mean_A_.resize(rows, cols);
	file.read((char *)mean_A_.data(), rows*cols*sizeof(Eigen::MatrixXd::Scalar));


	rows = 0, cols = 0;
	file.read((char*)(&rows), sizeof(int16_t));
	file.read((char*)(&cols), sizeof(int16_t));
	assert(static_cast<int>(rows) == MeshCuboidFeatures::k_num_features);
	assert(static_cast<int>(cols) == 1);
	mean_b_.resize(rows);
	file.read((char *)mean_b_.data(), rows*cols*sizeof(Eigen::MatrixXd::Scalar));


	rows = 0, cols = 0;
	file.read((char*)(&rows), sizeof(int16_t));
	file.read((char*)(&cols), sizeof(int16_t));
	assert(static_cast<int>(rows) == MeshCuboidFeatures::k_num_features);
	assert(static_cast<int>(cols) == MeshCuboidFeatures::k_num_features);
	inv_cov_.resize(rows, cols);
	file.read((char *)inv_cov_.data(), rows*cols*sizeof(Eigen::MatrixXd::Scalar));

	file.close();

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(inv_cov_);
	Real min_eigenvalue = es.eigenvalues().minCoeff();
	if (min_eigenvalue < -1.0E-6)
	{
		std::cerr << "Error: The inverse covariance matrix is not positive-semidefinite: ("
			<< _filename << " = " << min_eigenvalue << ")" << std::endl;
		do {
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	Real symmetry_diff = (inv_cov_ - inv_cov_.transpose()).array().abs().sum();
	if (symmetry_diff > 1.0E-6)
	{
		std::cerr << "Error: The inverse covariance matrix is not symmetric: ("
			<< _filename << " = " << symmetry_diff << ")" << std::endl;
		do {
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	return true;
}

void MeshCuboidCondNormalRelations::get_pairwise_cuboid_features(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
	Eigen::VectorXd &_global_features_vec_1, Eigen::VectorXd &_transformed_features_vec_12)
{
	assert(_cuboid_1);
	assert(_cuboid_2);

	MeshCuboidFeatures features_1;
	features_1.compute_features(_cuboid_1);

	MeshCuboidFeatures features_2;
	features_2.compute_features(_cuboid_2);

	get_pairwise_cuboid_features(features_1, features_2,
		_transformation_1, _transformation_2,
		_global_features_vec_1, _transformed_features_vec_12);
}

void MeshCuboidCondNormalRelations::get_pairwise_cuboid_features(
	const MeshCuboidFeatures &_features_1, const MeshCuboidFeatures &_features_2,
	const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2,
	Eigen::VectorXd &_global_features_vec_1, Eigen::VectorXd &_transformed_features_vec_12)
{
	assert(_transformation_1);
	assert(_transformation_2);

	Eigen::VectorXd features_vec_1 = _features_1.get_features();
	assert(features_vec_1.rows() == MeshCuboidFeatures::k_num_features);

	_transformed_features_vec_12 = _transformation_1->get_transformed_features(_features_2);
	assert(_transformed_features_vec_12.rows() == MeshCuboidFeatures::k_num_features);

	_global_features_vec_1 = features_vec_1.bottomRows(MeshCuboidFeatures::k_num_global_feature_values);
}

double MeshCuboidCondNormalRelations::compute_error(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2,
	const MeshCuboidTransformation *_transformation_1, const MeshCuboidTransformation *_transformation_2)const
{
	Eigen::VectorXd global_features_vec_1, transformed_features_vec_12;
	get_pairwise_cuboid_features(_cuboid_1, _cuboid_2, _transformation_1, _transformation_2,
		global_features_vec_1, transformed_features_vec_12);

	assert(mean_A_.rows() == transformed_features_vec_12.rows());
	assert(mean_A_.cols() == global_features_vec_1.rows());
	assert(mean_b_.rows() == transformed_features_vec_12.rows());
	assert(inv_cov_.rows() == transformed_features_vec_12.rows());
	assert(inv_cov_.cols() == transformed_features_vec_12.rows());

	const Eigen::VectorXd mean_2 = mean_A_ * global_features_vec_1 + mean_b_;
	Eigen::VectorXd diff = transformed_features_vec_12 - mean_2;

	// Mahalanobis norm.
	double error = diff.transpose() * inv_cov_ * diff;
	assert(error >= 0);

	//std::cerr << "Negative error value (error = " << error << ")" << std::endl;
	//Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ",");

	//Real inv_cov_det = inv_cov_.determinant();
	//std::cout << "inv_cov_det = " << inv_cov_det << std::endl;

	//std::cout << "Norm(diff) = " << diff.transpose() * diff << std::endl;

	//std::cout << "diff = " << std::endl;
	//std::cout << diff.format(csv_format) << std::endl;

	//do {
	//	std::cout << '\n' << "Press the Enter key to continue.";
	//} while (std::cin.get() != '\n');

	return error;
}

/*
MeshCuboidPCARelations::MeshCuboidPCARelations()
{
	mean_ = Eigen::VectorXd::Zero(k_mat_size);
	pca_bases_ = Eigen::MatrixXd::Zero(k_mat_size, k_mat_size);
}

MeshCuboidPCARelations::~MeshCuboidPCARelations()
{

}

bool MeshCuboidPCARelations::load_pca_csv(const char* _filename)
{
	Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ", ", "", "", "", "", "");

	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	std::string buffer;
	std::stringstream strstr;
	std::string token;

	std::getline(file, buffer);
	strstr.str(std::string());
	strstr.clear();
	strstr.str(buffer);

	for (int j = 0; j < k_mat_size; ++j)
	{
		if (strstr.eof())
		{
			std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}

		std::getline(strstr, token, ',');
		mean_(j) = std::stof(token);
	}

	for (int i = 0; i < k_mat_size; ++i)
	{
		if (file.eof())
		{
			std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}

		std::getline(file, buffer);
		strstr.str(std::string());
		strstr.clear();
		strstr.str(buffer);

		for (int j = 0; j < k_mat_size; ++j)
		{
			if (strstr.eof())
			{
				std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
				return false;
			}

			std::getline(strstr, token, ',');
			pca_bases_(i, j) = std::stof(token);
		}
	}

	file.close();

	return true;
}

bool MeshCuboidPCARelations::save_pca_csv(const char* _filename)const
{
	Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ",");

	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}
	std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	std::cout << std::scientific;

	file << mean_.transpose().format(csv_format) << std::endl;
	file << pca_bases_.format(csv_format) << std::endl;

	file.close();
	return true;
}

double MeshCuboidPCARelations::compute_error(
	const Eigen::VectorXd &_features_vec_1, const Eigen::VectorXd &_features_vec_2) const
{
	assert(_features_vec_1.rows() == MeshCuboidFeatures::k_num_features);
	assert(_features_vec_2.rows() == MeshCuboidFeatures::k_num_features);

	Eigen::VectorXd concat_features(2 * MeshCuboidFeatures::k_num_features);
	concat_features << _features_vec_1, _features_vec_2;
	assert(mean_.rows() == concat_features.rows());
	assert(pca_bases_.rows() == concat_features.rows());
	assert(pca_bases_.cols() == concat_features.rows());

	Eigen::VectorXd diff = concat_features - mean_;
	diff = diff - (pca_bases_ * diff);

	double error = diff.transpose() * diff;
	return error;
}

MeshCuboidCCARelations::MeshCuboidCCARelations()
{
	mean_1_ = Eigen::VectorXd::Zero(MeshCuboidFeatures::k_num_features);
	mean_2_ = Eigen::VectorXd::Zero(MeshCuboidFeatures::k_num_features);
	correlations_ = Eigen::VectorXd::Zero(MeshCuboidFeatures::k_num_features);
	bases_12_ = Eigen::MatrixXd::Zero(MeshCuboidFeatures::k_num_features, MeshCuboidFeatures::k_num_features);
	bases_21_ = Eigen::MatrixXd::Zero(MeshCuboidFeatures::k_num_features, MeshCuboidFeatures::k_num_features);
}

MeshCuboidCCARelations::~MeshCuboidCCARelations()
{

}

bool MeshCuboidCCARelations::load_cca_bases(const char* _filename)
{
	Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ",");

	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	std::string buffer;
	std::stringstream strstr;
	std::string token;

	std::getline(file, buffer);
	strstr.str(std::string());
	strstr.clear();
	strstr.str(buffer);

	for (int j = 0; j < MeshCuboidFeatures::k_num_features; ++j)
	{
		if (strstr.eof())
		{
			std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}

		std::getline(strstr, token, ',');
		mean_1_(j) = std::stof(token);
	}

	std::getline(file, buffer);
	strstr.str(std::string());
	strstr.clear();
	strstr.str(buffer);

	for (int j = 0; j < MeshCuboidFeatures::k_num_features; ++j)
	{
		if (strstr.eof())
		{
			std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}

		std::getline(strstr, token, ',');
		mean_2_(j) = std::stof(token);
	}

	std::getline(file, buffer);
	strstr.str(std::string());
	strstr.clear();
	strstr.str(buffer);

	for (int j = 0; j < MeshCuboidFeatures::k_num_features; ++j)
	{
		if (strstr.eof())
		{
			std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}

		std::getline(strstr, token, ',');
		correlations_(j) = std::stof(token);
	}

	for (int i = 0; i < MeshCuboidFeatures::k_num_features; ++i)
	{
		if (file.eof())
		{
			std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}

		std::getline(file, buffer);
		strstr.str(std::string());
		strstr.clear();
		strstr.str(buffer);

		for (int j = 0; j < MeshCuboidFeatures::k_num_features; ++j)
		{
			if (strstr.eof())
			{
				std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
				return false;
			}

			std::getline(strstr, token, ',');
			bases_12_(i, j) = std::stof(token);
		}
	}

	for (int i = 0; i < MeshCuboidFeatures::k_num_features; ++i)
	{
		if (file.eof())
		{
			std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}

		std::getline(file, buffer);
		strstr.str(std::string());
		strstr.clear();
		strstr.str(buffer);

		for (int j = 0; j < MeshCuboidFeatures::k_num_features; ++j)
		{
			if (strstr.eof())
			{
				std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
				return false;
			}

			std::getline(strstr, token, ',');
			bases_21_(i, j) = std::stof(token);
		}
	}

	file.close();
	return true;
}

bool MeshCuboidCCARelations::save_cca_bases(const char* _filename)const
{
	Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ",");

	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}
	std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	std::cout << std::scientific;

	file << mean_1_.transpose().format(csv_format) << std::endl;
	file << mean_2_.transpose().format(csv_format) << std::endl;
	file << correlations_.transpose().format(csv_format) << std::endl;
	file << bases_12_.format(csv_format) << std::endl;
	file << bases_21_.format(csv_format) << std::endl;


	file.close();
	return true;
}

double MeshCuboidCCARelations::compute_error(
	const Eigen::VectorXd &_features_vec_1, const Eigen::VectorXd &_features_vec_2) const
{
	assert(bases_12_.cols() == _features_vec_1.rows());
	assert(bases_21_.cols() == _features_vec_2.rows());

	Eigen::VectorXd transformed_attributes_1 = bases_12_ * (_features_vec_1 - mean_1_);
	Eigen::VectorXd transformed_attributes_2 = bases_21_ * (_features_vec_2 - mean_2_);
	Eigen::VectorXd diff = transformed_attributes_1 - transformed_attributes_2;

	// Mahalanobis norm.
	double error = diff.transpose() * correlations_.asDiagonal() * diff;
	return error;
}
*/