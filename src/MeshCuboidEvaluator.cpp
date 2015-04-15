#include "MeshCuboidEvaluator.h"

#include "ICP.h"

#include <fstream>
#include <iostream>


MeshCuboidEvaluator::MeshCuboidEvaluator(
	const MeshCuboidStructure *_ground_truth_cuboid_structure)
	//const std::string _mesh_name,
	//const std::string _cuboid_structure_name)
	: ground_truth_cuboid_structure_(_ground_truth_cuboid_structure)
	//, mesh_name_(_mesh_name)
	//, cuboid_structure_name_(_cuboid_structure_name)
	, mesh_name_("")
	, cuboid_structure_name_("")
{
	assert(ground_truth_cuboid_structure_);
}

bool MeshCuboidEvaluator::save_evaluate_results(
	const MeshCuboidStructure *_test_cuboid_structure,
	const char *_filename, bool _verbose)
{
	assert(_test_cuboid_structure);

	std::ofstream file(_filename);
	if (!file)
	{
		std::cerr << "Can't save file: \"" << _filename << "\"" << std::endl;
		return false;
	}

	if (_verbose)
		std::cout << "Saving " << _filename << "..." << std::endl;

	evaluate_all(_test_cuboid_structure);
	
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

void MeshCuboidEvaluator::evaluate_all(
	const MeshCuboidStructure *_test_cuboid_structure)
{
	assert(_test_cuboid_structure);

	//evaluate_segmentation(_test_cuboid_structure);
	//evaluate_cuboid_distance(_test_cuboid_structure);
}

/*
void MeshCuboidEvaluator::evaluate_segmentation(
	const MeshCuboidStructure *_test_cuboid_structure)
{
	assert(_test_cuboid_structure);

	const MyMesh *mesh = _test_cuboid_structure->mesh_;
	assert(mesh);

	const unsigned int num_labels = _test_cuboid_structure->num_labels();
	assert(_test_cuboid_structure->label_cuboids_.size() == num_labels);

	unsigned int num_sample_points = _test_cuboid_structure->num_sample_points();
	unsigned int num_cuboid_sample_points = 0;
	unsigned int num_correct_sample_point_labels = 0;
	//int num_changed_sample_point_labels = 0;

	for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
	{
		Label label = _test_cuboid_structure->get_label(label_index);

		for (std::vector<MeshCuboid *>::const_iterator c_it = _test_cuboid_structure->label_cuboids_[label_index].begin();
			c_it != _test_cuboid_structure->label_cuboids_[label_index].end(); ++c_it)
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

void MeshCuboidEvaluator::evaluate_cuboid_distance(
	const MeshCuboidStructure *_test_cuboid_structure)
{
	assert(_test_cuboid_structure);

	const MyMesh *mesh = _test_cuboid_structure->mesh_;
	assert(mesh);

	const unsigned int num_labels = _test_cuboid_structure->num_labels();

	// FIXME:
	// The cuboid structure should not deep copy all sample points.
	MeshCuboidStructure ground_truth_cuboid_structure = _test_cuboid_structure;
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

		if (_test_cuboid_structure->label_cuboids_[label_index].empty())
		{
			++num_missing_labels;
		}
		else
		{
			MeshCuboid *cuboid = _test_cuboid_structure->label_cuboids_[label_index].front();
			assert(cuboid);

			Real cuboid_distance = MeshCuboid::distance_between_cuboids(ground_truth_cuboid, cuboid);
			max_cuboid_distance = std::max(max_cuboid_distance, cuboid_distance);
		}
	}

	evaluation_results_["num_missing_labels"] = static_cast<Real>(num_missing_labels);
	evaluation_results_["max_cuboid_distance"] = max_cuboid_distance;
}
*/

void MeshCuboidEvaluator::evaluate_point_to_point_distances(
	const MeshCuboidStructure *_test_cuboid_structure,
	const char *_filename)
{
	assert(_test_cuboid_structure);

	const unsigned int num_neighbor_ranges = 1000;
	Real max_neighbor_range = 0.1;
	Eigen::VectorXd neighbor_ranges;
	neighbor_ranges.setLinSpaced(num_neighbor_ranges, 0.0, max_neighbor_range);

	const unsigned int num_labels = ground_truth_cuboid_structure_->num_labels();

	for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
	{
		MeshCuboid *ground_truth_cuboid = NULL, *test_cuboid = NULL;

		// NOTE:
		// The current implementation assumes that there is only one part for each label.
		assert(ground_truth_cuboid_structure_->label_cuboids_[label_index].size() <= 1);
		if (!ground_truth_cuboid_structure_->label_cuboids_[label_index].empty())
			ground_truth_cuboid = ground_truth_cuboid_structure_->label_cuboids_[label_index].front();

		assert(_test_cuboid_structure->label_cuboids_[label_index].size() <= 1);
		if (!_test_cuboid_structure->label_cuboids_[label_index].empty())
			test_cuboid = _test_cuboid_structure->label_cuboids_[label_index].front();

		if (!ground_truth_cuboid || !test_cuboid)
			continue;

		const unsigned int num_ground_truth_sample_points = ground_truth_cuboid->num_sample_points();
		const unsigned int num_test_sample_points = test_cuboid->num_sample_points();

		if (num_ground_truth_sample_points == 0 || num_test_sample_points == 0)
			continue;


		// Create a ground truth sample point KD-tree.
		Eigen::MatrixXd ground_truth_sample_points(3, num_ground_truth_sample_points);
		for (SamplePointIndex sample_point_index = 0; sample_point_index < num_ground_truth_sample_points;
			++sample_point_index)
		{
			assert(ground_truth_cuboid->get_sample_point(sample_point_index));
			for (unsigned int i = 0; i < 3; ++i)
				ground_truth_sample_points.col(sample_point_index)(i) =
				ground_truth_cuboid->get_sample_point(sample_point_index)->point_[i];
		}

		ANNpointArray ground_truth_sample_ann_points;
		ANNkd_tree *ground_truth_sample_ann_kd_tree = ICP::create_kd_tree(ground_truth_sample_points,
			ground_truth_sample_ann_points);
		assert(ground_truth_sample_ann_points);
		assert(ground_truth_sample_ann_kd_tree);


		// Create a test sample point KD-tree.
		Eigen::MatrixXd test_sample_points(3, num_test_sample_points);
		for (SamplePointIndex sample_point_index = 0; sample_point_index < num_test_sample_points;
			++sample_point_index)
		{
			assert(test_cuboid->get_sample_point(sample_point_index));
			for (unsigned int i = 0; i < 3; ++i)
				test_sample_points.col(sample_point_index)(i) =
				test_cuboid->get_sample_point(sample_point_index)->point_[i];
		}

		ANNpointArray test_sample_ann_points;
		ANNkd_tree *test_sample_ann_kd_tree = ICP::create_kd_tree(test_sample_points,
			test_sample_ann_points);
		assert(test_sample_ann_points);
		assert(test_sample_ann_kd_tree);


		// Ground truth -> test.
		Eigen::VectorXd ground_truth_to_test_distances;
		ICP::get_closest_points(test_sample_ann_kd_tree, ground_truth_sample_points,
			ground_truth_to_test_distances);
		assert(ground_truth_to_test_distances.rows() == num_ground_truth_sample_points);

		// Test -> ground truth.
		Eigen::VectorXd test_to_ground_truth_distances;
		ICP::get_closest_points(ground_truth_sample_ann_kd_tree, test_sample_points,
			test_to_ground_truth_distances);
		assert(test_to_ground_truth_distances.rows() == num_test_sample_points);


		Real distance_error = (ground_truth_to_test_distances.squaredNorm() / num_ground_truth_sample_points)
			+ (test_to_ground_truth_distances.squaredNorm() / num_test_sample_points);
		assert(distance_error >= 0);


		Eigen::MatrixXd precision_recall(2, num_neighbor_ranges);

		for (unsigned int i = 0; i < num_neighbor_ranges; ++i)
		{
			// Precision.
			precision_recall.row(0)(i) = static_cast<double>(
				(test_to_ground_truth_distances.array() <= neighbor_ranges[i]).count())
				/ num_test_sample_points;

			// Recall.
			precision_recall.row(1)(i) = static_cast<double>(
				(ground_truth_to_test_distances.array() <= neighbor_ranges[i]).count())
				/ num_ground_truth_sample_points;
		}

		std::stringstream sstr;
		sstr << _filename << "_" << label_index << ".csv";
		std::ofstream file(sstr.str());

		if (!file.good())
		{
			do {
				std::cout << "Error: Cannot save file: '" << sstr.str() << "'." << std::endl;
				std::cout << '\n' << "Press the Enter key to continue.";
			} while (std::cin.get() != '\n');
		}

		Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ",");
		file << precision_recall.format(csv_format) << std::endl;
		file.close();
		

		annDeallocPts(ground_truth_sample_ann_points);
		annDeallocPts(test_sample_ann_points);
		delete ground_truth_sample_ann_kd_tree;
		delete test_sample_ann_kd_tree;
	}
}
