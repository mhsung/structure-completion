#include "MeshCuboidTrainer.h"

#include <Eigen/Core>
#include <Eigen/LU>
#include <QFileInfo>

#include <fstream>
#include <iostream>
#include <sstream>


MeshCuboidTrainer::MeshCuboidTrainer()
{

}

void MeshCuboidTrainer::clear()
{
	for (std::vector< std::list<MeshCuboidFeatures *> >::iterator f_it = feature_list_.begin();
		f_it != feature_list_.end(); ++f_it)
	{
		for (std::list<MeshCuboidFeatures *>::iterator f_jt = (*f_it).begin();
			f_jt != (*f_it).end(); ++f_jt)
			delete (*f_jt);
		(*f_it).clear();
	}
	feature_list_.clear();

	for (std::vector< std::list<MeshCuboidTransformation *> >::iterator t_it = transformation_list_.begin();
		t_it != transformation_list_.end(); ++t_it)
	{
		for (std::list<MeshCuboidTransformation *>::iterator t_jt = (*t_it).begin();
			t_jt != (*t_it).end(); ++t_jt)
			delete (*t_jt);
		(*t_it).clear();
	}
	transformation_list_.clear();
}

bool MeshCuboidTrainer::load_object_list(const std::string &_filename)
{
	std::ifstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't load file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	object_list_.clear();
	std::string buffer;

	while (!file.eof())
	{
		std::getline(file, buffer);
		if (buffer == "") break;
		object_list_.push_back(buffer);
	}

	return true;
}

bool MeshCuboidTrainer::load_features(const std::string &_filename_prefix)
{
	for (std::vector< std::list<MeshCuboidFeatures *> >::iterator f_it = feature_list_.begin();
		f_it != feature_list_.end(); ++f_it)
	{
		for (std::list<MeshCuboidFeatures *>::iterator f_jt = (*f_it).begin();
			f_jt != (*f_it).end(); ++f_jt)
			delete (*f_jt);
		(*f_it).clear();
	}
	feature_list_.clear();


	for (unsigned int cuboid_index = 0; true; ++cuboid_index)
	{
		std::stringstream sstr;
		sstr << _filename_prefix << std::to_string(cuboid_index) << std::string(".csv");
		std::string attributes_filename = sstr.str();

		QFileInfo attributes_file(attributes_filename.c_str());
		if (!attributes_file.exists())
			break;

		std::cout << "Loading '" << attributes_filename << "'..." << std::endl;

		std::list<MeshCuboidFeatures *> stats;
		MeshCuboidFeatures::load_feature_collection(
			attributes_filename.c_str(), stats);

		feature_list_.push_back(stats);
	}

	return true;
}

bool MeshCuboidTrainer::load_transformations(const std::string &_filename_prefix)
{
	for (std::vector< std::list<MeshCuboidTransformation *> >::iterator t_it = transformation_list_.begin();
		t_it != transformation_list_.end(); ++t_it)
	{
		for (std::list<MeshCuboidTransformation *>::iterator t_jt = (*t_it).begin();
			t_jt != (*t_it).end(); ++t_jt)
			delete (*t_jt);
		(*t_it).clear();
	}
	transformation_list_.clear();


	for (unsigned int cuboid_index = 0; true; ++cuboid_index)
	{
		std::stringstream  sstr;
		sstr << _filename_prefix << std::to_string(cuboid_index) << std::string(".csv");
		std::string transformation_filename = sstr.str();

		QFileInfo transformation_file(transformation_filename.c_str());
		if (!transformation_file.exists())
			break;

		std::cout << "Loading '" << transformation_filename << "'..." << std::endl;

		std::list<MeshCuboidTransformation *> stats;
		MeshCuboidTransformation::load_transformation_collection(
			transformation_filename.c_str(), stats);

		transformation_list_.push_back(stats);
	}

	return true;
}

MeshCuboidJointNormalRelationTrainer::MeshCuboidJointNormalRelationTrainer()
	: MeshCuboidTrainer()
{

}

void MeshCuboidJointNormalRelationTrainer::compute_relations(
	std::vector< std::vector<MeshCuboidJointNormalRelations *> > &_relations,
	const std::list<std::string> *_ignored_object_list) const
{
	unsigned int num_labels = feature_list_.size();
	assert(transformation_list_.size() == num_labels);

	_relations.clear();
	_relations.resize(num_labels);
	for (unsigned int cuboid_index = 0; cuboid_index < num_labels; ++cuboid_index)
		_relations[cuboid_index].resize(num_labels, NULL);


	for (unsigned int label_index_1 = 0; label_index_1 < num_labels; ++label_index_1)
	{
		for (unsigned int label_index_2 = 0; label_index_2 < num_labels; ++label_index_2)
		{
			if (label_index_1 == label_index_2) continue;

			// NOTE:
			// 'object_list_' should contain all object names.
			assert(object_list_.size() == feature_list_[label_index_1].size());
			assert(object_list_.size() == feature_list_[label_index_2].size());
			assert(object_list_.size() == transformation_list_[label_index_1].size());
			assert(object_list_.size() == transformation_list_[label_index_2].size());

			std::vector<MeshCuboidFeatures *> feature_1;
			std::vector<MeshCuboidFeatures *> feature_2;
			std::vector<MeshCuboidTransformation *> transformation_1;
			std::vector<MeshCuboidTransformation *> transformation_2;

			feature_1.reserve(feature_list_[label_index_1].size());
			feature_2.reserve(feature_list_[label_index_2].size());
			transformation_1.reserve(transformation_list_[label_index_1].size());
			transformation_2.reserve(transformation_list_[label_index_2].size());

			std::list<std::string>::const_iterator o_it = object_list_.begin();
			std::list<MeshCuboidFeatures *>::const_iterator f_it_1 = feature_list_[label_index_1].begin();
			std::list<MeshCuboidFeatures *>::const_iterator f_it_2 = feature_list_[label_index_2].begin();
			std::list<MeshCuboidTransformation *>::const_iterator t_it_1 = transformation_list_[label_index_1].begin();
			std::list<MeshCuboidTransformation *>::const_iterator t_it_2 = transformation_list_[label_index_2].begin();

			int num_objects = 0;
			while (true)
			{
				if (f_it_1 == feature_list_[label_index_1].end()
					|| f_it_2 == feature_list_[label_index_2].end()
					|| t_it_1 == transformation_list_[label_index_1].end()
					|| t_it_2 == transformation_list_[label_index_2].end())
					break;

				bool has_values = (!(*f_it_1)->has_nan() && !(*f_it_2)->has_nan());

				if (has_values && _ignored_object_list)
				{
					// Check whether the current object should be ignored.
					for (std::list<std::string>::const_iterator io_it = _ignored_object_list->begin();
						io_it != _ignored_object_list->end(); ++io_it)
					{
						if ((*o_it) == (*io_it))
						{
							std::cout << "Mesh [" << (*o_it) << "] is ignored." << std::endl;
							has_values = false;
							break;
						}
					}
				}

				if (has_values)
				{
					feature_1.push_back(*f_it_1);
					feature_2.push_back(*f_it_2);
					transformation_1.push_back(*t_it_1);
					transformation_2.push_back(*t_it_2);
					++num_objects;
				}

				++o_it;
				++f_it_1;
				++f_it_2;
				++t_it_1;
				++t_it_2;
			}

			if (num_objects == 0) continue;


			_relations[label_index_1][label_index_2] = new MeshCuboidJointNormalRelations();
			MeshCuboidJointNormalRelations *relation_12 = _relations[label_index_1][label_index_2];
			assert(relation_12);

			assert(feature_1.size() == num_objects);
			assert(feature_2.size() == num_objects);
			assert(transformation_1.size() == num_objects);
			assert(transformation_2.size() == num_objects);

			Eigen::MatrixXd X_1(num_objects, MeshCuboidFeatures::k_num_features);
			Eigen::MatrixXd X_2(num_objects, MeshCuboidFeatures::k_num_features);

			for (int object_index = 0; object_index < num_objects; ++object_index)
			{
				assert(feature_1[object_index]);
				assert(feature_2[object_index]);
				assert(transformation_1[object_index]);
				assert(transformation_2[object_index]);

				Eigen::VectorXd transformed_feature_1 = transformation_2[object_index]->get_transformed_features(
					(*feature_1[object_index]));

				Eigen::VectorXd transformed_feature_2 = transformation_1[object_index]->get_transformed_features(
					(*feature_2[object_index]));

				assert(transformed_feature_1.size() == MeshCuboidFeatures::k_num_features);
				assert(transformed_feature_2.size() == MeshCuboidFeatures::k_num_features);
				
				X_1.row(object_index) = transformed_feature_1.transpose();
				X_2.row(object_index) = transformed_feature_2.transpose();
			}

			Eigen::MatrixXd X(num_objects, 2 * MeshCuboidFeatures::k_num_features);
			X << X_1, X_2;

			Eigen::RowVectorXd mean_X = X.colwise().mean();
			Eigen::MatrixXd centered_X = X.rowwise() - mean_X;

			Eigen::MatrixXd cov = (centered_X.transpose() * centered_X) / static_cast<double>(num_objects);
			cov = cov + 1.0E-3 * Eigen::MatrixXd::Identity(cov.rows(), cov.cols());
			Eigen::MatrixXd inv_cov = cov.inverse();

			relation_12->set_mean(mean_X.transpose());
			relation_12->set_inv_cov(inv_cov);

			Eigen::MatrixXd diff = (X.rowwise() - mean_X).transpose();
			Eigen::VectorXd error = (diff.transpose() * inv_cov * diff).diagonal();
			std::cout << "(" << label_index_1 << ", " << label_index_2 << "): max_error = " << error.maxCoeff() << std::endl;
		}
	}
}
