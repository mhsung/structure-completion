#include "MeshCuboidEvaluator.h"

#include <fstream>
#include <iostream>


MeshCuboidEvaluator::MeshCuboidEvaluator(const MeshCuboidStructure &_cuboid_structure,
	const std::string _mesh_name,
	const std::string _cuboid_structure_name)
	: cuboid_structure_(_cuboid_structure)
	, mesh_name_(_mesh_name)
	, cuboid_structure_name_(_cuboid_structure_name)
{
}

bool MeshCuboidEvaluator::save_evaluate_results(const std::string _filename, bool _verbose)
{
	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Saving " << _filename << "..." << std::endl;


	evaluate_all();

	file << "Mesh name," << mesh_name_ << std::endl;
	file << "Cuboid structure name," << cuboid_structure_name_ << std::endl;

	if (_verbose)
	{
		std::cout << "Mesh name: " << mesh_name_ << std::endl;
		std::cout << "Cuboid structure name: " << cuboid_structure_name_ << std::endl;
	}

	for (std::map<std::string, Real>::const_iterator it = evaluation_results_.begin();
		it != evaluation_results_.end(); ++it)
	{
		const std::string key = (*it).first;
		Real value = (*it).second;

		file << key << "," << value << std::endl;

		if (_verbose)
			std::cout << key << ": " << value << std::endl;
	}

	file.close();
	return true;
}

void MeshCuboidEvaluator::evaluate_all()
{
	evaluate_segmentation();
	evaluate_cuboid_distance();
}

void MeshCuboidEvaluator::evaluate_segmentation()
{
	const MyMesh *mesh = cuboid_structure_.mesh_;
	assert(mesh);

	const unsigned int num_labels = cuboid_structure_.num_labels();
	assert(cuboid_structure_.label_cuboids_.size() == num_labels);

	unsigned int num_sample_points = cuboid_structure_.num_sample_points();
	unsigned int num_cuboid_sample_points = 0;
	unsigned int num_correct_sample_point_labels = 0;
	//int num_changed_sample_point_labels = 0;

	for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
	{
		Label label = cuboid_structure_.get_label(label_index);

		for (std::vector<MeshCuboid *>::const_iterator c_it = cuboid_structure_.label_cuboids_[label_index].begin();
			c_it != cuboid_structure_.label_cuboids_[label_index].end(); ++c_it)
		{
			const MeshCuboid *cuboid = (*c_it);
			assert(cuboid);
			assert(cuboid->get_label_index() == label_index);

			const std::vector<MeshSamplePoint *> cuboid_sample_points = cuboid->get_sample_points();
			for (std::vector<MeshSamplePoint *>::const_iterator s_it = cuboid_sample_points.begin();
				s_it != cuboid_sample_points.end(); ++s_it)
			{
				const MeshSamplePoint *sample_point = (*s_it);
				assert(sample_point);

				//SamplePointIndex sample_point_index = sample_point->sample_point_index_;
				//assert(sample_point_index < _sample_point_label_indices.size());
				//assert(_cuboid_structure.sample_points_[sample_point_index] == sample_point);

				//if (label_index != _sample_point_label_indices[sample_point_index])
				//{
				//	_sample_point_label_indices[sample_point_index] = label_index;
				//	++num_changed_sample_point_labels;
				//}

				assert(sample_point->corr_fid_ < mesh->n_faces());
				MyMesh::FaceHandle fh = mesh->face_handle(sample_point->corr_fid_);
				Label face_label = mesh->property(mesh->face_label_, fh);

				if (face_label == label)
					++num_correct_sample_point_labels;

				++num_cuboid_sample_points;
			}
		}
	}

	evaluation_results_["num_sample_points"] = static_cast<Real>(num_sample_points);

	Real segmented_sample_point_ratio = static_cast<Real>(num_cuboid_sample_points) / num_sample_points;
	evaluation_results_["segmented_sample_point_ratio"] = segmented_sample_point_ratio;

	Real correct_sample_point_label_ratio = static_cast<Real>(num_correct_sample_point_labels) / num_cuboid_sample_points;
	evaluation_results_["correct_sample_point_label_ratio"] = correct_sample_point_label_ratio;
}

void MeshCuboidEvaluator::evaluate_cuboid_distance()
{
	const MyMesh *mesh = cuboid_structure_.mesh_;
	assert(mesh);

	const unsigned int num_labels = cuboid_structure_.num_labels();

	// FIXME:
	// The cuboid structure should not deep copy all sample points.
	MeshCuboidStructure ground_truth_cuboid_structure = cuboid_structure_;
	ground_truth_cuboid_structure.clear_cuboids();
	ground_truth_cuboid_structure.get_mesh_face_label_cuboids();

	
	unsigned int num_missing_labels = 0;
	Real max_cuboid_distance = 0;

	for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
	{
		// NOTE:
		// The current implementation assumes that there is only one part for each label.
		assert(ground_truth_cuboid_structure.label_cuboids_[label_index].size() <= 1);
		if (ground_truth_cuboid_structure.label_cuboids_[label_index].empty())
			continue;
		
		MeshCuboid *ground_truth_cuboid = ground_truth_cuboid_structure.label_cuboids_[label_index].front();
		assert(ground_truth_cuboid);

		if (cuboid_structure_.label_cuboids_[label_index].empty())
		{
			++num_missing_labels;
		}
		else
		{
			MeshCuboid *cuboid = cuboid_structure_.label_cuboids_[label_index].front();
			assert(cuboid);

			Real cuboid_distance = MeshCuboid::distance_between_cuboids(ground_truth_cuboid, cuboid);
			max_cuboid_distance = std::max(max_cuboid_distance, cuboid_distance);
		}
	}

	evaluation_results_["num_missing_labels"] = static_cast<Real>(num_missing_labels);
	evaluation_results_["max_cuboid_distance"] = max_cuboid_distance;
}
