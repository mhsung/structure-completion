#include "MeshCuboidStats.h"
#include "MeshCuboidRelation.h"
#include "Parameters.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>


double standard_normal_distribution(double _z);


MeshCuboidStats::MeshCuboidStats()
{

}

MeshCuboidStats::~MeshCuboidStats()
{

}

void MeshCuboidStats::initialize()
{
	keys_.clear();
	stats_.clear();
}

unsigned int MeshCuboidStats::num_stats() const
{
	assert(keys_.size() == stats_.size());
	return keys_.size();
}

bool MeshCuboidStats::load_stats(const char* _filename, bool _verbose)
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Loading " << _filename << "..." << std::endl;

	initialize();

	std::string buffer;
	std::stringstream strstr;
	std::string token;
	std::string category;

	keys_.clear();
	stats_.clear();

	std::getline(file, buffer);

	strstr.str(std::string());
	strstr.clear();
	strstr.str(buffer);

	// NOTE:
	// The first column is not a key.
	std::getline(strstr, category, ',');

	while (std::getline(strstr, token, ','))
		keys_.push_back(token);

	if (keys_.empty())
	{
		std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
		return false;
	}
	stats_.resize(keys_.size());


	while (!file.eof())
	{
		std::getline(file, buffer);
		if (buffer == "") break;

		strstr.str(std::string());
		strstr.clear();
		strstr.str(buffer);

		std::getline(strstr, category, ',');

		unsigned int count = 0;
		for (; std::getline(strstr, token, ',') && count < stats_.size(); ++count)
		{
			if (category == "Count")
			{
				stats_[count].count_ = atoi(token.c_str());
			}
			else if (category == "Mean")
			{
				stats_[count].mean_ = atof(token.c_str());
			}
			else if (category == "Stdev")
			{
				stats_[count].stdev_ = atof(token.c_str());
			}
			else
			{
				assert(false);
			}
		}

		if (count < stats_.size())
		{
			std::cerr << "Wrong file format: \"" << _filename << "\"" << std::endl;
			return false;
		}
	}


	if (_verbose)
		std::cout << "Done." << std::endl;

	file.close();
	return true;
}

bool MeshCuboidStats::save_stats(const char* _filename, bool _verbose)
{
	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Saving " << _filename << "..." << std::endl;

	std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	std::cout << std::scientific;

	assert(keys_.size() == stats_.size());


	// Write attribute keys.
	file << ",";
	for (std::vector<std::string>::const_iterator it = keys_.begin();
		it != keys_.end(); ++it)
	{
		file << (*it) << ",";
	}
	file << std::endl;


	file << "Count,";
	for (std::vector<StatParams>::const_iterator it = stats_.begin();
		it != stats_.end(); ++it)
	{
		file << (*it).count_ << ",";
	}
	file << std::endl;

	file << "Mean,";
	for (std::vector<StatParams>::const_iterator it = stats_.begin();
		it != stats_.end(); ++it)
	{
		file << (*it).mean_ << ",";
	}
	file << std::endl;

	file << "Stdev,";
	for (std::vector<StatParams>::const_iterator it = stats_.begin();
		it != stats_.end(); ++it)
	{
		file << (*it).stdev_ << ",";
	}
	file << std::endl;

	if (_verbose)
		std::cout << "Done." << std::endl;

	file.close();
	return true;
}

void MeshCuboidStats::compute_stats(const std::vector<std::string>& _keys,
	const std::vector<std::vector<Real>>& _values)
{
	keys_.clear();
	stats_.clear();

	unsigned int num_keys = _keys.size();
	unsigned int num_objects = _values.size();

	keys_ = _keys;
	stats_.resize(num_keys);

	for (unsigned int key_index = 0; key_index < num_keys; ++key_index)
	{
		stats_[key_index].count_ = 0;
		stats_[key_index].mean_ = 0;
		stats_[key_index].stdev_ = 0;
	}

	// Count.
	for (unsigned int key_index = 0; key_index < num_keys; ++key_index)
	{
		unsigned int count = 0;
		for (unsigned int object_index = 0; object_index < num_objects; ++object_index)
		{
			assert(_values[object_index].size() == num_keys);
			Real value = _values[object_index][key_index];
			if (!std::isnan(value))
			{
				++count;
			}
		}
		stats_[key_index].count_ = count;
	}

	// Mean.
	for (unsigned int key_index = 0; key_index < num_keys; ++key_index)
	{
		unsigned int count = stats_[key_index].count_;

		if (count == 0)
		{
			// NOTE:
			// Should be NaN.
			stats_[key_index].mean_ = 0;
		}
		else
		{
			Real sum_values = 0.0;
			for (unsigned int object_index = 0; object_index < num_objects; ++object_index)
			{
				assert(_values[object_index].size() == num_keys);
				Real value = _values[object_index][key_index];
				if (!std::isnan(value))
				{
					sum_values += value;
				}
			}

			Real mean = sum_values / static_cast<Real>(count);
			stats_[key_index].mean_ = mean;
		}
	}

	// Stdev.
	for (unsigned int key_index = 0; key_index < num_keys; ++key_index)
	{
		unsigned int count = stats_[key_index].count_;

		if (count == 0)
		{
			// NOTE:
			// Should be NaN.
			stats_[key_index].stdev_ = 0;
		}
		else
		{
			Real mean = stats_[key_index].mean_;

			Real sq_sum_values = 0.0;
			for (unsigned int object_index = 0; object_index < num_objects; ++object_index)
			{
				assert(_values[object_index].size() == num_keys);
				Real value = _values[object_index][key_index];
				if (!std::isnan(value))
				{
					sq_sum_values += (value * value);
				}
			}

			Real stdev = std::sqrt(sq_sum_values / static_cast<Real>(count)-mean * mean);
			stats_[key_index].stdev_ = stdev;
		}
	}
}


MeshCuboidManualFeatures::MeshCuboidManualFeatures(const std::string _object_name)
	: object_name_(_object_name)
{
}

MeshCuboidManualFeatures::~MeshCuboidManualFeatures()
{

}

void MeshCuboidManualFeatures::clear()
{
	keys_.clear();
	values_.clear();
}

bool MeshCuboidManualFeatures::get_keys_and_values(const std::list<MeshCuboidManualFeatures *>& _stats,
	std::vector<std::string>& _keys,
	std::vector<std::vector<Real>>& _values)
{
	_keys.clear();
	_values.clear();

	unsigned int num_objects = _stats.size();
	if (num_objects == 0) return true;

	_values.resize(num_objects);


	std::list<MeshCuboidManualFeatures *>::const_iterator stat_it = _stats.begin();
	assert(*stat_it);
	_keys = (*stat_it)->keys_;
	++stat_it;

	for (; stat_it != _stats.end(); ++stat_it)
	{
		assert(*stat_it);
		if (!std::equal(_keys.begin(), _keys.end(), (*stat_it)->keys_.begin()))
		{
			std::cerr << "Error: Stats do not share the same set of keys." << std::endl;
			do {
				std::cout << '\n' << "Press the Enter key to continue.";
			} while (std::cin.get() != '\n');
			return false;
		}
	}

	unsigned int num_keys = _keys.size();


	stat_it = _stats.begin();
	unsigned int object_index = 0;
	for (; stat_it != _stats.end(); ++stat_it, ++object_index)
	{
		assert(*stat_it);
		_values[object_index].resize(num_keys, 0);
		assert((*stat_it)->values_.size() == num_keys);

		for (unsigned int key_index = 0; key_index < num_keys; ++key_index)
		{
			_values[object_index][key_index] = (*stat_it)->values_[key_index];
		}
	}

	return true;
}

void MeshCuboidManualFeatures::get_object_names(const std::list<MeshCuboidManualFeatures *>& _stats,
	std::vector<std::string>& _object_names)
{
	unsigned int num_objects = _stats.size();

	_object_names.clear();
	_object_names.resize(num_objects);

	unsigned int object_index = 0;
	for (std::list<MeshCuboidManualFeatures *>::const_iterator stat_it = _stats.begin(); stat_it != _stats.end();
		++stat_it, ++object_index)
	{
		assert(*stat_it);
		_object_names[object_index] = (*stat_it)->object_name_;
	}
}

bool MeshCuboidManualFeatures::save_keys_and_values(
	const std::list<MeshCuboidManualFeatures *>& _stats, const char* _filename)
{
	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}
	std::setprecision(std::numeric_limits<long double>::digits10 + 1);
	std::cout << std::scientific;

	std::vector<std::string> keys;
	std::vector<std::vector<Real>> values;
	std::vector<std::string> object_names;

	bool ret = get_keys_and_values(_stats, keys, values);
	assert(ret);
	get_object_names(_stats, object_names);

	unsigned int num_keys = keys.size();
	unsigned int num_objects = values.size();
	assert(num_objects == object_names.size());


	// Write attribute keys.
	file << ",";
	for (std::vector<std::string>::const_iterator key_it = keys.begin();
		key_it != keys.end(); ++key_it)
	{
		file << (*key_it) << ",";
	}
	file << std::endl;


	// Write attribute values.
	for (unsigned int object_index = 0; object_index < num_objects; ++object_index)
	{
		file << object_names[object_index] << ",";
		assert(values[object_index].size() == num_keys);
		for (unsigned int key_index = 0; key_index < num_keys; ++key_index)
		{
			Real value = values[object_index][key_index];
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

bool MeshCuboidManualFeatures::save_stats(
	const std::list<MeshCuboidManualFeatures *>& _stats, const char* _filename)
{
	bool ret;
	std::vector<std::string> keys;
	std::vector<std::vector<Real>> values;

	ret = get_keys_and_values(_stats, keys, values);
	assert(ret);

	MeshCuboidStats stats;
	stats.compute_stats(keys, values);
	ret = stats.save_stats(_filename);
	if (!ret) return false;

	return true;
}

Real MeshCuboidManualFeatures::compute_log_probability(const MeshCuboidStats& _stats)const
{
	const Real half_bandwidth = 0.01;
	const Real min_probability = 1.0E-12;


	if (!std::equal(keys_.begin(), keys_.end(), _stats.keys_.begin()))
	{
		std::cerr << "Error: Stats do not share the same set of keys." << std::endl;
		do {
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
		return false;
	}

	assert(keys_.size() == values_.size());
	assert(values_.size() == _stats.stats_.size());
	unsigned int num_keys = keys_.size();
	

	Real accumulated_log_probability = 0.0;

	for (unsigned int key_index = 0; key_index < num_keys; ++key_index)
	{
		const Real value = values_[key_index];

		MeshCuboidStats::StatParams stat_params = _stats.stats_[key_index];
		const unsigned int count = stat_params.count_;
		const Real mean = stat_params.mean_;
		const Real stdev = stat_params.stdev_;

		// Gaussian distribution.
		// Compute only when the number of samples is greater than zero,
		// and standard deviation is greater than zero.
		if (count > 0 && stdev > 0)
		{
			// Gaussian distribution.
			Real z = (value - mean) / stdev;
			Real probability = standard_normal_distribution(z + half_bandwidth)
				- standard_normal_distribution(z - half_bandwidth);
			assert(probability >= 0.0);
			assert(probability <= 1.0);

			// Avoid that log(probability) becomes -infinity.
			probability = std::max(probability, min_probability);
			Real log_probability = std::log(probability);
			assert(log_probability <= 0.0);
			accumulated_log_probability += log_probability;
		}
	}

	assert(accumulated_log_probability <= 0.0);
	return accumulated_log_probability;
}

MeshSingleCuboidManualFeatures::MeshSingleCuboidManualFeatures(const std::string _object_name)
	: MeshCuboidManualFeatures(_object_name)
{
}

MeshSingleCuboidManualFeatures::~MeshSingleCuboidManualFeatures()
{
}

//std::string MeshSingleCuboidFeatures::get_key_bbox_axes(
//	const Label _label, const unsigned int _axis_index, const unsigned int _i) const
//{
//	std::stringstream key;
//	// Do not compute probability.
//	key << k_distribution_none
//		<< "_label_" << _label << "_axis_" << _axis_index << "_" << _i;
//	return key.str();
//}
//
//std::string MeshSingleCuboidFeatures::get_key_bbox_center(
//	const Label _label, const unsigned int _axis_index) const
//{
//	std::stringstream key;
//	// Do not compute probability.
//	key << k_distribution_none
//		<< "_label_" << _label << "_center_" << _axis_index;
//	return key.str();
//}

std::string MeshSingleCuboidManualFeatures::get_key_bbox_size(
	const unsigned int _axis_index) const
{
	std::stringstream key;
	//key << "_label_" << _label;
	key	<< "_size_" << _axis_index;
	return key.str();
}

//std::string MeshSingleCuboidFeatures::get_key_bbox_height(
//	const unsigned int _i) const
//{
//	std::stringstream key;
//	//key << "_label_" << _label;
//	key << "_height_" << _i;
//	return key.str();
//}

std::string MeshSingleCuboidManualFeatures::get_key_bbox_corner_height(
	const unsigned int _corner_index) const
{
	std::stringstream key;
	//key << "_label_" << _label;
	key << "_corner_height_" << _corner_index;
	return key.str();
}

std::string MeshSingleCuboidManualFeatures::get_key_bbox_axes_up_angle(
	const unsigned int _axis_index) const
{
	std::stringstream key;
	//key << "_label_" << _label;
	key << "_axis_up_angle_" << _axis_index;
	return key.str();
}

//Real MeshSingleCuboidFeatures::get_value_bbox_axes(
//	const MeshCuboid *_cuboid, const unsigned int _axis_index, const unsigned int _i) const
//{
//	MyMesh::Normal axis = _cuboid->get_bbox_axis(_axis_index).normalized();
//	Real value = axis[_i];
//	return value;
//}
//
//Real MeshSingleCuboidFeatures::get_value_bbox_center(
//	const MeshCuboid *_cuboid, const unsigned int _axis_index) const
//{
//	Real value = _cuboid->get_bbox_center()[_axis_index];
//	return value;
//}

Real MeshSingleCuboidManualFeatures::get_value_bbox_size(
	const MeshCuboid *_cuboid, const unsigned int _axis_index) const
{
	Real value = _cuboid->get_bbox_size()[_axis_index];

#ifdef COMPUTE_LINEAR_MAPS
	// Linear map.
	const int num_attributes = static_cast<int>(MeshCuboidAttributes::NUM_ATTRIBUTES);
	Eigen::VectorXd linear_map_A(num_attributes);
	linear_map_A.setZero();
	linear_map_A[static_cast<int>(MeshCuboidAttributes::SIZE_X) + _axis_index] = 1.0;
#endif

	return value;
}

//Real MeshSingleCuboidFeatures::get_value_bbox_height(
//	const MeshCuboid *_cuboid, const unsigned int _i) const
//{
//	std::array<MyMesh::Point, MeshCuboid::k_num_corners> corners = _cuboid->get_bbox_corners();
//	const unsigned int num_corners = MeshCuboid::k_num_corners;
//	Real value;
//
//	if (_i == 0)		value = std::numeric_limits<Real>::max();	// Min.
//	else if (_i == 1)	value = -std::numeric_limits<Real>::max();	// Max.
//	assert(_i < 2);
//
//	for (unsigned int corner_index = 0; corner_index < num_corners; ++corner_index)
//	{
//		MyMesh::Point corner = corners[corner_index];
//		Real height = dot(corner, k_up_direction);
//
//		if (_i == 0)		value = std::min(value, height);
//		else if (_i == 1)	value = std::max(value, height);
//	}
//
//	return value;
//}

Real MeshSingleCuboidManualFeatures::get_value_bbox_corner_height(
	const MeshCuboid *_cuboid, const unsigned int _corner_index) const
{
	MyMesh::Normal axis_direction;
	MyMesh::Point corner = _cuboid->get_bbox_corner(_corner_index, &axis_direction);
	Real value = dot(corner, k_up_direction);

#ifdef COMPUTE_LINEAR_MAPS
	// Linear map.
	Eigen::Vector3d up_direction_vec;
	Eigen::Vector3d axis_direction_vec;
	Eigen::Matrix3d bbox_axes_mat;

	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		up_direction_vec[axis_index] = k_up_direction[axis_index];
		axis_direction_vec[axis_index] = axis_direction[axis_index];
		for (unsigned int i = 0; i < 3; i++)
			bbox_axes_mat.col(axis_index)[i] = _cuboid->get_bbox_axis(axis_index)[i];
	}

	const int num_attributes = static_cast<int>(MeshCuboidAttributes::NUM_ATTRIBUTES);
	Eigen::VectorXd linear_map_A(num_attributes);
	linear_map_A.setZero();
	linear_map_A.block(0, 0, 3, 1) = -up_direction_vec;
	linear_map_A.block(3, 0, 3, 1) =
		(-0.5 * axis_direction_vec.asDiagonal()) * bbox_axes_mat.transpose() * up_direction_vec;
#endif

	return value;
}

Real MeshSingleCuboidManualFeatures::get_value_bbox_axis_up_angle(
	const MeshCuboid *_cuboid, const unsigned int _axis_index) const
{
	MyMesh::Normal axis = _cuboid->get_bbox_axis(_axis_index).normalized();
	Real value = dot(axis, k_up_direction);
	return value;
}

void MeshSingleCuboidManualFeatures::compute_values(const MeshCuboid *_cuboid)
{
	clear();

	// NOTE:
	// For efficiency, allocate memory for keys and values.
	const unsigned int num_keys = 3 + 3 + 2;
	keys_.reserve(num_keys);
	values_.reserve(num_keys);

	std::string key;
	Real value;

	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		//// Bounding box axes.
		//MyMesh::Normal axis = _cuboid->get_bbox_axis(axis_index);
		//for (unsigned int i = 0; i < 3; ++i)
		//{
		//	key = get_key_bbox_axes(axis_index, i);
		//	value = get_value_bbox_axes(_cuboid, axis_index, i);
		//	assert(maps_.find(key) == maps_.end());
		//	maps_[key] = value;
		//}

		//// Bounding box center.
		//key = get_key_bbox_center(axis_index);
		//value = get_value_bbox_center(_cuboid, axis_index);
		//assert(maps_.find(key) == maps_.end());
		//maps_[key] = value;

		// Bounding box size.
		key = get_key_bbox_size(axis_index);
		value = get_value_bbox_size(_cuboid, axis_index);
		keys_.push_back(key);
		values_ .push_back(value);

		// Bounding box axis up-direction angle.
		key = get_key_bbox_axes_up_angle(axis_index);
		value = get_value_bbox_axis_up_angle(_cuboid, axis_index);
		keys_.push_back(key);
		values_.push_back(value);
	}

	//// Bounding box min/max height.
	//for (unsigned int i = 0; i < 2; ++i)
	//{
	//	key = get_key_bbox_height(i);
	//	value = get_value_bbox_height(_cuboid, i);
	//	keys_.push_back(key);
	//	values_.push_back(value);
	//}

	const unsigned int num_corners = MeshCuboid::k_num_corners;
	for (unsigned int corner_index = 0; corner_index < num_corners; ++corner_index)
	{
		key = get_key_bbox_corner_height(corner_index);
		value = get_value_bbox_corner_height(_cuboid, corner_index);
		keys_.push_back(key);
		values_.push_back(value);
	}

	assert(keys_.size() == values_.size());
}

/*
void MeshSingleCuboidFeatures::compute_values(MeshCuboidStructure *_cuboid_structure)
{
	assert(_cuboid_structure);
	assert(_cuboid_structure->label_cuboids_.size() == _cuboid_structure->num_labels());

	for (LabelIndex label_index = 0; label_index < _cuboid_structure->num_labels(); ++label_index)
	{
		// NOTE:
		// The current implementation assumes that there is only one part for each label.
		assert(_cuboid_structure->label_cuboids_[label_index].size() <= 1);
		if (_cuboid_structure->label_cuboids_[label_index].empty())
			continue;

		//Label label = _cuboid_structure->get_label(label_index);
		MeshCuboid *cuboid = _cuboid_structure->label_cuboids_[label_index].front();
		assert(cuboid);
		assert(cuboid->get_label_index() == label_index);

		//compute_values(label, cuboid);
		compute_values(cuboid);
	}
}
*/

MeshPairCuboidManualFeatures::MeshPairCuboidManualFeatures(const std::string _object_name)
	: MeshCuboidManualFeatures(_object_name)
{
}

MeshPairCuboidManualFeatures::~MeshPairCuboidManualFeatures()
{
}

std::string MeshPairCuboidManualFeatures::get_key_bbox_axis_angle(
	const unsigned int _axis_index_1, const unsigned int _axis_index_2) const
{
	std::stringstream key;
	//key << "_label_" << _label_1 << "_" << _label_2;
	key << "_axis_angle_(axis_axis)_" << _axis_index_1 << "_" << _axis_index_2;
	return key.str();
}

//std::string MeshPairCuboidFeatures::get_key_bbox_size_ratio(
//	const unsigned int _axis_index_1, const unsigned int _axis_index_2) const
//{
//	// FIXME:
//	// This should be a ratio distribution, not a Gaussian distribution.
//	std::stringstream key;
//	//key << "_label_" << _label_1 << "_" << _label_2;
//	key << "_size_ratio_(axis_axis)_" << _axis_index_1 << "_" << _axis_index_2;
//	return key.str();
//}

std::string MeshPairCuboidManualFeatures::get_key_bbox_face_corner_distance(
	const CuboidPairOrder _pair_order,
	const unsigned int _face_index_1, const unsigned int _corner_index_2) const
{
	std::stringstream key;
	//key << "_label_" << _label_1 << "_" << _label_2;
	switch (_pair_order)
	{
	case CUBOID_1_2: key << "_12"; break;
	case CUBOID_2_1: key << "_21"; break;
	}
	key	<< "_face_corner_dist_(face_corner)_" << _face_index_1 << "_" << _corner_index_2;
	return key.str();
}

std::string MeshPairCuboidManualFeatures::get_key_bbox_axis_center_distance(
	const CuboidPairOrder _pair_order,
	const unsigned int _axis_index_1) const
{
	std::stringstream key;
	//key << "_label_" << _label_1 << "_" << _label_2;
	switch (_pair_order)
	{
		case CUBOID_1_2: key << "_12"; break;
		case CUBOID_2_1: key << "_21"; break;
	}
	key	<< "_axis_center_dist_(axis)_" << _axis_index_1;
	return key.str();
}

Real MeshPairCuboidManualFeatures::get_value_bbox_axis_angle(
	const MeshCuboid* _cuboid_1, const MeshCuboid* _cuboid_2,
	const unsigned int _axis_index_1, const unsigned int _axis_index_2) const
{
	MyMesh::Normal axis_1 = _cuboid_1->get_bbox_axis(_axis_index_1).normalized();
	MyMesh::Normal axis_2 = _cuboid_2->get_bbox_axis(_axis_index_2).normalized();

	Real cos_value = dot(axis_1, axis_2);
	assert(std::abs(cos_value) <= 1 + param_zero_tol);
	cos_value = std::min(std::max(cos_value, -1.0), 1.0);
	Real value = acos(cos_value) / M_PI;

	return value;
}

//Real MeshPairCuboidFeatures::get_value_bbox_size_ratio(
//	const MeshCuboid* _cuboid_1, const MeshCuboid* _cuboid_2,
//	const unsigned int _axis_index_1, const unsigned int _axis_index_2) const
//{
//	double size_1 = _cuboid_1->get_bbox_size()[_axis_index_1];
//	double size_2 = _cuboid_2->get_bbox_size()[_axis_index_2];
//	Real value = size_1 / size_2;
//
//	return value;
//}

Real MeshPairCuboidManualFeatures::get_value_bbox_face_corner_distance(
	const MeshCuboid* _cuboid_1, const MeshCuboid* _cuboid_2,
	const unsigned int _face_index_1, const unsigned int _corner_index_2) const
{
	MyMesh::Normal axis_direction_1, axis_direction_2;

	std::pair<MyMesh::Normal, MyMesh::Point> face_1 =
		_cuboid_1->get_bbox_faces(_face_index_1, &axis_direction_1);
	MyMesh::Normal face_normal_1 = face_1.first.normalized();
	MyMesh::Point face_point_1 = face_1.second;

	MyMesh::Point corner_2 =
		_cuboid_2->get_bbox_corner(_corner_index_2, &axis_direction_2);

	Real value = dot(corner_2 - face_point_1, face_normal_1);

#ifdef COMPUTE_LINEAR_MAPS
	// Linear map.
	Eigen::Vector3d axis_direction_vec_1, axis_direction_vec_2;
	Eigen::Matrix3d bbox_axes_mat_1, bbox_axes_mat_2;

	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		axis_direction_vec_1[axis_index] = axis_direction_1[axis_index];
		axis_direction_vec_2[axis_index] = axis_direction_2[axis_index];
		for (unsigned int i = 0; i < 3; i++)
		{
			bbox_axes_mat_1.col(axis_index)[i] = _cuboid_1->get_bbox_axis(axis_index)[i];
			bbox_axes_mat_2.col(axis_index)[i] = _cuboid_2->get_bbox_axis(axis_index)[i];
		}
	}

	Eigen::Vector3d face_normal_vec_1 = bbox_axes_mat_1 * axis_direction_vec_1;

	const int num_attributes = static_cast<int>(MeshCuboidAttributes::NUM_ATTRIBUTES);
	Eigen::VectorXd linear_map_A(2 * num_attributes);
	linear_map_A.setZero();
	Real linear_map_b = 0;
	
	linear_map_A.block(0, 0, 3, 1) = -face_normal_vec_1;
	linear_map_A.block(3, 0, 3, 1) =
		(-0.5 * axis_direction_vec_1.asDiagonal()) * bbox_axes_mat_1.transpose() * face_normal_vec_1;
	linear_map_A.block(6, 0, 3, 1) = face_normal_vec_1;
	linear_map_A.block(9, 0, 3, 1) =
		(0.5 * axis_direction_vec_2.asDiagonal()) * bbox_axes_mat_2.transpose() * face_normal_vec_1;
#endif

	return value;
}

Real MeshPairCuboidManualFeatures::get_value_bbox_axis_center_distance(
	const MeshCuboid* _cuboid_1, const MeshCuboid* _cuboid_2,
	const unsigned int _axis_index_1) const
{
	MyMesh::Point center_1 = _cuboid_1->get_bbox_center();
	MyMesh::Point center_2 = _cuboid_2->get_bbox_center();
	MyMesh::Normal axis_1 = _cuboid_1->get_bbox_axis(_axis_index_1).normalized();
	Real value = dot(center_2 - center_1, axis_1);

#ifdef COMPUTE_LINEAR_MAPS
	// Linear map.
	Eigen::Vector3d axis_vec_1;
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		axis_vec_1[axis_index] = axis_1[axis_index];

	const int num_attributes = static_cast<int>(MeshCuboidAttributes::NUM_ATTRIBUTES);
	Eigen::VectorXd linear_map_A(2 * num_attributes);
	linear_map_A.setZero();
	Real linear_map_b = 0;

	linear_map_A.block(0, 0, 3, 1) = -axis_vec_1;
	linear_map_A.block(6, 0, 3, 1) = axis_vec_1;
#endif

	return value;
}

void MeshPairCuboidManualFeatures::compute_values_bbox_axis_angle(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2)
{
	// NOTE:
	// Different-parts-only relationships.
	//assert(_label_1 != _label_2);

	// NOTE:
	// Two-way relationships.
	// label_2_1 is the same with label_1_2.
	// Have only one case of them.
	//assert(_label_1 <= _label_2);

	for (unsigned int axis_index_1 = 0; axis_index_1 < 3; ++axis_index_1)
	{
		for (unsigned int axis_index_2 = 0; axis_index_2 < 3; ++axis_index_2)
		{
			std::string key = get_key_bbox_axis_angle(axis_index_1, axis_index_2);
			Real value = get_value_bbox_axis_angle(
				_cuboid_1, _cuboid_2, axis_index_1, axis_index_2);
			keys_.push_back(key);
			values_.push_back(value);
		}
	}
}

//void MeshPairCuboidFeatures::compute_values_bbox_size_ratio(
//	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2)
//{
//	// NOTE:
//	// Two-way relationships.
//	// label_2_1 is the same with label_1_2.
//	// Have only one case of them.
//	//assert(_label_1 <= _label_2);
//
//	for (unsigned int axis_index_1 = 0; axis_index_1 < 3; ++axis_index_1)
//	{
//		for (unsigned int axis_index_2 = 0; axis_index_2 < 3; ++axis_index_2)
//		{
//			std::string key = get_key_bbox_size_ratio(axis_index_1, axis_index_2);
//			Real value = get_value_bbox_size_ratio(
//				_cuboid_1, _cuboid_2, axis_index_1, axis_index_2);
//			keys_.push_back(key);
//			values_.push_back(value);
//		}
//	}
//}

void MeshPairCuboidManualFeatures::compute_values_bbox_face_corner_distance(
	const CuboidPairOrder _pair_order,
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2)
{
	// NOTE:
	// Different-parts-only relationships.
	//assert(_label_1 != _label_2);

	// NOTE:
	// One-way relationships.

	const unsigned int num_faces_1 = _cuboid_1->get_bbox_faces().size();
	const unsigned int num_corners_2 = MeshCuboid::k_num_corners;

	for (unsigned int face_index_1 = 0; face_index_1 < num_faces_1; ++face_index_1)
	{
		for (unsigned int corner_index_2 = 0; corner_index_2 < num_corners_2; ++corner_index_2)
		{
			std::string key = get_key_bbox_face_corner_distance(
				_pair_order, face_index_1, corner_index_2);
			Real value = get_value_bbox_face_corner_distance(
				_cuboid_1, _cuboid_2, face_index_1, corner_index_2);
			keys_.push_back(key);
			values_.push_back(value);
		}
	}
}

void MeshPairCuboidManualFeatures::compute_values_bbox_axis_center_distance(
	const CuboidPairOrder _pair_order,
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2)
{
	// NOTE:
	// Different-parts-only relationships.
	//assert(_label_1 != _label_2);

	// NOTE:
	// One-way relationships.

	for (unsigned int axis_index_1 = 0; axis_index_1 < 3; ++axis_index_1)
	{
		std::string key = get_key_bbox_axis_center_distance(
			_pair_order, axis_index_1);
		Real value = get_value_bbox_axis_center_distance(
			_cuboid_1, _cuboid_2, axis_index_1);
		keys_.push_back(key);
		values_.push_back(value);
	}
}

void MeshPairCuboidManualFeatures::compute_values(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2)
{
	clear();

	// NOTE:
	// For efficiency, allocate memory for keys and values.
	const unsigned int num_keys = 3*3 + 3*3 + 2*48 + 2*3;
	keys_.reserve(num_keys);
	values_.reserve(num_keys);


	//compute_values_bbox_size_ratio(_cuboid_1, _cuboid_2);

	// NOTE:
	// Different-parts-only relationships.
	//if (_label_1 != _label_2)
	{
		compute_values_bbox_axis_angle(_cuboid_1, _cuboid_2);

		// NOTE:
		// One-way relationships.
		// Make both label_1_2 and label_2_1 cases.
		compute_values_bbox_face_corner_distance(CUBOID_1_2, _cuboid_1, _cuboid_2);
		compute_values_bbox_face_corner_distance(CUBOID_2_1, _cuboid_2, _cuboid_1);

		compute_values_bbox_axis_center_distance(CUBOID_1_2, _cuboid_1, _cuboid_2);
		compute_values_bbox_axis_center_distance(CUBOID_2_1, _cuboid_2, _cuboid_1);
	}

	assert(keys_.size() == values_.size());
}

/*
void MeshPairCuboidFeatures::compute_values(MeshCuboidStructure *_cuboid_structure)
{
	assert(_cuboid_structure);
	assert(_cuboid_structure->label_cuboids_.size() == _cuboid_structure->num_labels());

	for (LabelIndex label_index_1 = 0; label_index_1 < _cuboid_structure->num_labels();
		++label_index_1)
	{
		// NOTE:
		// The current implementation assumes that there is only one part for each label.
		assert(_cuboid_structure->label_cuboids_[label_index_1].size() <= 1);
		if (_cuboid_structure->label_cuboids_[label_index_1].empty())
			continue;

		//Label label_1 = _cuboid_structure->get_label(label_index_1);
		MeshCuboid *cuboid_1 = _cuboid_structure->label_cuboids_[label_index_1].front();
		assert(cuboid_1);
		assert(cuboid_1->get_label_index() == label_index_1);

		for (LabelIndex label_index_2 = label_index_1; label_index_2 < _cuboid_structure->num_labels();
			++label_index_2)
		{
			// NOTE:
			// The current implementation assumes that there is only one part for each label.
			assert(_cuboid_structure->label_cuboids_[label_index_2].size() <= 1);
			if (_cuboid_structure->label_cuboids_[label_index_2].empty())
				continue;

			//Label label_2 = _cuboid_structure->get_label(label_index_2);
			MeshCuboid *cuboid_2 = _cuboid_structure->label_cuboids_[label_index_2].front();
			assert(cuboid_2);
			assert(cuboid_2->get_label_index() == label_index_2);

			//compute_values(label_1, cuboid_1, label_2, cuboid_2);
			compute_values(cuboid_1, cuboid_2);
		}
	}
}
*/

double standard_normal_distribution(double _z)
{
	const double table[31][10] = {
			{ 0.0000, 0.0040, 0.0080, 0.0120, 0.0160, 0.0199, 0.0239, 0.0279, 0.0319, 0.0359 },
			{ 0.0398, 0.0438, 0.0478, 0.0517, 0.0557, 0.0596, 0.0636, 0.0675, 0.0714, 0.0753 },
			{ 0.0793, 0.0832, 0.0871, 0.0910, 0.0948, 0.0987, 0.1026, 0.1064, 0.1103, 0.1141 },
			{ 0.1179, 0.1217, 0.1255, 0.1293, 0.1331, 0.1368, 0.1406, 0.1443, 0.1480, 0.1517 },
			{ 0.1554, 0.1591, 0.1628, 0.1664, 0.1700, 0.1736, 0.1772, 0.1808, 0.1844, 0.1879 },
			{ 0.1915, 0.1950, 0.1985, 0.2019, 0.2054, 0.2088, 0.2123, 0.2157, 0.2190, 0.2224 },
			{ 0.2257, 0.2291, 0.2324, 0.2357, 0.2389, 0.2422, 0.2454, 0.2486, 0.2517, 0.2549 },
			{ 0.2580, 0.2611, 0.2642, 0.2673, 0.2704, 0.2734, 0.2764, 0.2794, 0.2823, 0.2852 },
			{ 0.2881, 0.2910, 0.2939, 0.2967, 0.2995, 0.3023, 0.3051, 0.3078, 0.3106, 0.3133 },
			{ 0.3159, 0.3186, 0.3212, 0.3238, 0.3264, 0.3289, 0.3315, 0.3340, 0.3365, 0.3389 },
			{ 0.3413, 0.3438, 0.3461, 0.3485, 0.3508, 0.3531, 0.3554, 0.3577, 0.3599, 0.3621 },
			{ 0.3643, 0.3665, 0.3686, 0.3708, 0.3729, 0.3749, 0.3770, 0.3790, 0.3810, 0.3830 },
			{ 0.3849, 0.3869, 0.3888, 0.3907, 0.3925, 0.3944, 0.3962, 0.3980, 0.3997, 0.4015 },
			{ 0.4032, 0.4049, 0.4066, 0.4082, 0.4099, 0.4115, 0.4131, 0.4147, 0.4162, 0.4177 },
			{ 0.4192, 0.4207, 0.4222, 0.4236, 0.4251, 0.4265, 0.4279, 0.4292, 0.4306, 0.4319 },
			{ 0.4332, 0.4345, 0.4357, 0.4370, 0.4382, 0.4394, 0.4406, 0.4418, 0.4429, 0.4441 },
			{ 0.4452, 0.4463, 0.4474, 0.4484, 0.4495, 0.4505, 0.4515, 0.4525, 0.4535, 0.4545 },
			{ 0.4554, 0.4564, 0.4573, 0.4582, 0.4591, 0.4599, 0.4608, 0.4616, 0.4625, 0.4633 },
			{ 0.4641, 0.4649, 0.4656, 0.4664, 0.4671, 0.4678, 0.4686, 0.4693, 0.4699, 0.4706 },
			{ 0.4713, 0.4719, 0.4726, 0.4732, 0.4738, 0.4744, 0.4750, 0.4756, 0.4761, 0.4767 },
			{ 0.4772, 0.4778, 0.4783, 0.4788, 0.4793, 0.4798, 0.4803, 0.4808, 0.4812, 0.4817 },
			{ 0.4821, 0.4826, 0.4830, 0.4834, 0.4838, 0.4842, 0.4846, 0.4850, 0.4854, 0.4857 },
			{ 0.4861, 0.4864, 0.4868, 0.4871, 0.4875, 0.4878, 0.4881, 0.4884, 0.4887, 0.4890 },
			{ 0.4893, 0.4896, 0.4898, 0.4901, 0.4904, 0.4906, 0.4909, 0.4911, 0.4913, 0.4916 },
			{ 0.4918, 0.4920, 0.4922, 0.4925, 0.4927, 0.4929, 0.4931, 0.4932, 0.4934, 0.4936 },
			{ 0.4938, 0.4940, 0.4941, 0.4943, 0.4945, 0.4946, 0.4948, 0.4949, 0.4951, 0.4952 },
			{ 0.4953, 0.4955, 0.4956, 0.4957, 0.4959, 0.4960, 0.4961, 0.4962, 0.4963, 0.4964 },
			{ 0.4965, 0.4966, 0.4967, 0.4968, 0.4969, 0.4970, 0.4971, 0.4972, 0.4973, 0.4974 },
			{ 0.4974, 0.4975, 0.4976, 0.4977, 0.4977, 0.4978, 0.4979, 0.4979, 0.4980, 0.4981 },
			{ 0.4981, 0.4982, 0.4982, 0.4983, 0.4984, 0.4984, 0.4985, 0.4985, 0.4986, 0.4986 },
			{ 0.4987, 0.4987, 0.4987, 0.4988, 0.4988, 0.4989, 0.4989, 0.4989, 0.4990, 0.4990 },
	};

	unsigned int label = (int)std::abs(_z * 100 + 0.5);
	unsigned int i = label / 10;
	unsigned int j = label % 10;

	double probability = 0.5;
	if (i <= 30)
	{
		probability = table[i][j];
	}

	if (_z > 0)	probability = 0.5 + probability;
	else		probability = 0.5 - probability;

	assert(probability >= 0 && probability <= 1);
	return probability;
}