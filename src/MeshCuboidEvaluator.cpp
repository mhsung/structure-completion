#include "MeshCuboidEvaluator.h"

#include "MeshCuboidParameters.h"
#include "ICP.h"

#include <fstream>
#include <iostream>


MeshCuboidEvaluator::MeshCuboidEvaluator(
	MeshCuboidStructure *_ground_truth_cuboid_structure)
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

void MeshCuboidEvaluator::evaluate_point_to_point_distances(
	const std::vector<MeshSamplePoint *> _ground_truth_sample_points,
	const std::vector<MeshSamplePoint *> _test_sample_points,
	const char *_filename, bool _record_error)
{
	const unsigned int num_ground_truth_sample_points = _ground_truth_sample_points.size();
	const unsigned int num_test_sample_points = _test_sample_points.size();

	if (num_ground_truth_sample_points == 0 || num_test_sample_points == 0)
		return;


	Eigen::VectorXd neighbor_ranges;
	neighbor_ranges.setLinSpaced(FLAGS_param_eval_num_neighbor_range_samples,
		0.0, FLAGS_param_eval_max_neighbor_distance);


	// Create a ground truth sample point KD-tree.
	Eigen::MatrixXd ground_truth_sample_points(3, num_ground_truth_sample_points);
	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_ground_truth_sample_points;
		++sample_point_index)
	{
		assert(_ground_truth_sample_points[sample_point_index]);
		for (unsigned int i = 0; i < 3; ++i)
			ground_truth_sample_points.col(sample_point_index)(i) =
			_ground_truth_sample_points[sample_point_index]->point_[i];
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
		assert(_test_sample_points[sample_point_index]);
		for (unsigned int i = 0; i < 3; ++i)
			test_sample_points.col(sample_point_index)(i) =
			_test_sample_points[sample_point_index]->point_[i];
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


	if (_record_error)
	{
		double eval_neighbor_distance_min_max_range = FLAGS_param_eval_max_neighbor_distance
			- FLAGS_param_eval_min_neighbor_distance;
		assert(eval_neighbor_distance_min_max_range);

		// Accuracy.
		for (SamplePointIndex sample_point_index = 0; sample_point_index < num_test_sample_points;
			++sample_point_index)
		{
			_test_sample_points[sample_point_index]->error_ = 0.0;

			double distance = test_to_ground_truth_distances[sample_point_index];
			if (distance >= FLAGS_param_eval_min_neighbor_distance)
			{
				double error = (distance - FLAGS_param_eval_min_neighbor_distance) / eval_neighbor_distance_min_max_range;
				error = std::min(error, 1.0);

				// Map (min, max) to 0.5, 1.0.
				error = 0.5 + error * 0.5;
				_test_sample_points[sample_point_index]->error_ = error;
			}
		}

		// Completeness.
		for (SamplePointIndex sample_point_index = 0; sample_point_index < num_ground_truth_sample_points;
			++sample_point_index)
		{
			_ground_truth_sample_points[sample_point_index]->error_ = 0.0;

			double distance = ground_truth_to_test_distances[sample_point_index];
			if (distance >= FLAGS_param_eval_min_neighbor_distance)
			{
				double error = (distance - FLAGS_param_eval_min_neighbor_distance) / eval_neighbor_distance_min_max_range;
				error = std::min(error, 1.0);

				// Map (min, max) to 0.5, 1.0.
				error = 0.5 + error * 0.5;
				_ground_truth_sample_points[sample_point_index]->error_ = error;
			}
		}
	}


	Eigen::VectorXd accuracy(FLAGS_param_eval_num_neighbor_range_samples);
	Eigen::VectorXd completeness(FLAGS_param_eval_num_neighbor_range_samples);

	for (unsigned int i = 0; i < FLAGS_param_eval_num_neighbor_range_samples; ++i)
	{
		accuracy[i] = static_cast<double>(
			(test_to_ground_truth_distances.array() <= neighbor_ranges[i]).count())
			/ num_test_sample_points;

		completeness[i] = static_cast<double>(
			(ground_truth_to_test_distances.array() <= neighbor_ranges[i]).count())
			/ num_ground_truth_sample_points;
	}

	std::ofstream file(_filename);
	if (!file.good())
	{
		do {
			std::cout << "Error: Cannot save file: '" << _filename << "'." << std::endl;
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	std::cout << "Saving '" << _filename << "'..." << std::endl;

	Eigen::IOFormat csv_format(Eigen::StreamPrecision, 0, ",");
	file << neighbor_ranges.transpose().format(csv_format) << std::endl;
	file << accuracy.transpose().format(csv_format) << std::endl;
	file << completeness.transpose().format(csv_format) << std::endl;
	file.close();


	annDeallocPts(ground_truth_sample_ann_points);
	annDeallocPts(test_sample_ann_points);
	delete ground_truth_sample_ann_kd_tree;
	delete test_sample_ann_kd_tree;
}

void MeshCuboidEvaluator::evaluate_point_to_point_distances(
	const MeshCuboidStructure *_test_cuboid_structure,
	const char *_filename)
{
	assert(_test_cuboid_structure);

	std::stringstream output_filename_sstr;


	// For entire object.
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << _filename << ".csv";

	evaluate_point_to_point_distances(
		ground_truth_cuboid_structure_->sample_points_,
		_test_cuboid_structure->sample_points_,
		output_filename_sstr.str().c_str(), true);


	// For each symmetric part sets.
	const unsigned int num_symmetric_label_sets
		= ground_truth_cuboid_structure_->label_symmetries_.size();

	
	for (unsigned int symmetric_label_set_index = 0;
		symmetric_label_set_index < num_symmetric_label_sets; ++symmetric_label_set_index)
	{
		const std::list<LabelIndex> &label_symmetry
			= ground_truth_cuboid_structure_->label_symmetries_[symmetric_label_set_index];

		std::vector<MeshSamplePoint *> ground_truth_sample_points;
		std::vector<MeshSamplePoint *> test_sample_points;

		for (std::list<LabelIndex>::const_iterator it = label_symmetry.begin(); it != label_symmetry.end(); ++it)
		{
			LabelIndex label_index = (*it);
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

		
			const std::vector<MeshSamplePoint *>& ground_truth_cuboid_sample_points
				= ground_truth_cuboid->get_sample_points();
			const std::vector<MeshSamplePoint *>& test_cuboid_sample_points
				= test_cuboid->get_sample_points();

			ground_truth_sample_points.insert(ground_truth_sample_points.end(),
				ground_truth_cuboid_sample_points.begin(), ground_truth_cuboid_sample_points.end());
			test_sample_points.insert(test_sample_points.end(),
				test_cuboid_sample_points.begin(), test_cuboid_sample_points.end());
		}

		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << _filename << "_" << symmetric_label_set_index << ".csv";

		evaluate_point_to_point_distances( ground_truth_sample_points, test_sample_points,
			output_filename_sstr.str().c_str());
	}
}

void MeshCuboidEvaluator::evaluate_point_labeling(
	const MeshCuboidStructure *_test_cuboid_structure,
	const char *_filename)
{
	assert(_test_cuboid_structure);

	const MyMesh *mesh = _test_cuboid_structure->mesh_;
	assert(mesh);

	std::ofstream file(_filename);
	if (!file.good())
	{
		do {
			std::cout << "Error: Cannot save file: '" << _filename << "'." << std::endl;
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	std::cout << "Saving '" << _filename << "'..." << std::endl;

	// For entire object.
	unsigned int num_total_labeled_samples = 0;
	unsigned int num_total_correctly_labeled_samples = 0;

	const unsigned int num_symmetric_label_sets
		= ground_truth_cuboid_structure_->label_symmetries_.size();

	for (unsigned int symmetric_label_set_index = 0;
		symmetric_label_set_index < num_symmetric_label_sets; ++symmetric_label_set_index)
	{
		// For each symmetric part sets.
		unsigned int num_labeled_samples = 0;
		unsigned int num_correctly_labeled_samples = 0;

		const std::list<LabelIndex> &label_symmetry
			= ground_truth_cuboid_structure_->label_symmetries_[symmetric_label_set_index];

		for (std::list<LabelIndex>::const_iterator it = label_symmetry.begin(); it != label_symmetry.end(); ++it)
		{
			LabelIndex label_index = (*it);
			MeshCuboid *test_cuboid = NULL;

			// NOTE:
			// The current implementation assumes that there is only one part for each label.
			assert(_test_cuboid_structure->label_cuboids_[label_index].size() <= 1);
			if (!_test_cuboid_structure->label_cuboids_[label_index].empty())
				test_cuboid = _test_cuboid_structure->label_cuboids_[label_index].front();

			if (!test_cuboid)
				continue;

			Label label = _test_cuboid_structure->get_label(label_index);

			const std::vector<MeshSamplePoint *>& test_cuboid_sample_points
				= test_cuboid->get_sample_points();
			for (std::vector<MeshSamplePoint *>::const_iterator s_it = test_cuboid_sample_points.begin();
				s_it != test_cuboid_sample_points.end(); ++s_it)
			{
				const MeshSamplePoint *sample_point = (*s_it);
				assert(sample_point);
				assert(sample_point->corr_fid_ < mesh->n_faces());
				MyMesh::FaceHandle fh = mesh->face_handle(sample_point->corr_fid_);
				Label face_label = mesh->property(mesh->face_label_, fh);

				if (face_label == label)
					++num_correctly_labeled_samples;

				++num_labeled_samples;

			}
		}

		file << symmetric_label_set_index << "," <<
			num_correctly_labeled_samples << "," <<
			num_labeled_samples << ",";

		if (num_labeled_samples == 0)
		{
			file << "none" << std::endl;
		}
		else
		{
			Real label_accuracy = static_cast<Real>(num_correctly_labeled_samples) / num_labeled_samples;
			file << label_accuracy << std::endl;
		}

		num_total_labeled_samples += num_labeled_samples;
		num_total_correctly_labeled_samples += num_correctly_labeled_samples;
	}

	file << "all," <<
		num_total_correctly_labeled_samples << "," <<
		num_total_labeled_samples << ",";
	
	if (num_total_labeled_samples == 0)
	{
		file << "none" << std::endl;
	}
	else
	{
		Real total_label_accuracy = static_cast<Real>(num_total_correctly_labeled_samples) /
			num_total_labeled_samples;
		file << total_label_accuracy << std::endl;
	}

	file.close();
}

void MeshCuboidEvaluator::evaluate_cuboid_distance(
	const MeshCuboidStructure *_test_cuboid_structure, const char *_filename)
{
	assert(_test_cuboid_structure);

	const MyMesh *mesh = _test_cuboid_structure->mesh_;
	assert(mesh);

	std::ofstream file(_filename);
	if (!file.good())
	{
		do {
			std::cout << "Error: Cannot save file: '" << _filename << "'." << std::endl;
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	std::cout << "Saving '" << _filename << "'..." << std::endl;

	// For entire object.
	Real total_max_cuboid_distance = -1;

	const unsigned int num_symmetric_label_sets
		= ground_truth_cuboid_structure_->label_symmetries_.size();

	for (unsigned int symmetric_label_set_index = 0;
		symmetric_label_set_index < num_symmetric_label_sets; ++symmetric_label_set_index)
	{
		// For each symmetric part sets.
		Real max_cuboid_distance = -1;

		const std::list<LabelIndex> &label_symmetry
			= ground_truth_cuboid_structure_->label_symmetries_[symmetric_label_set_index];

		for (std::list<LabelIndex>::const_iterator it = label_symmetry.begin(); it != label_symmetry.end(); ++it)
		{
			LabelIndex label_index = (*it);
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

			// NOTE: Cuboid surface points are replaced.
			ground_truth_cuboid->create_grid_points_on_cuboid_surface(FLAGS_param_num_cuboid_surface_points);
			test_cuboid->create_grid_points_on_cuboid_surface(FLAGS_param_num_cuboid_surface_points);

			Real cuboid_distance = MeshCuboid::distance_between_cuboids(ground_truth_cuboid, test_cuboid);
			max_cuboid_distance = std::max(cuboid_distance, max_cuboid_distance);
		}

		file << symmetric_label_set_index << ",";
		if (max_cuboid_distance < 0)
			file << "none" << std::endl;
		else
			file << max_cuboid_distance << std::endl;

		total_max_cuboid_distance = std::max(max_cuboid_distance, total_max_cuboid_distance);
	}

	file << "all,";
	if (total_max_cuboid_distance < 0)
		file << "none" << std::endl;
	else
		file << total_max_cuboid_distance << std::endl;

	file.close();
}