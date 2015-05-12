#include "MeshViewerCore.h"
#include "MeshCuboidEvaluator.h"
#include "MeshCuboidFusion.h"
#include "MeshCuboidParameters.h"
#include "MeshCuboidTrainer.h"

#include <Eigen/Core>
#include <Eigen/Geometry> 
#include <QDir>
#include <QFileInfo>


void get_bounding_cylinder(const MeshCuboidStructure &_cuboid_structure,
	Eigen::MatrixXd &_sample_points, Eigen::VectorXd &_bbox_center,
	Real &_xy_size, Real &_z_size)
{
	const unsigned int num_input_points = _cuboid_structure.sample_points_.size();
	assert(num_input_points > 0);
	_sample_points = Eigen::MatrixXd(3, num_input_points);

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_input_points;
		++sample_point_index)
	{
		const MeshSamplePoint *sample_point = _cuboid_structure.sample_points_[sample_point_index];
		assert(sample_point);
		for (int i = 0; i < 3; ++i)
			_sample_points.col(sample_point_index)[i] = sample_point->point_[i];
	}

	Eigen::VectorXd input_bbox_min = _sample_points.rowwise().minCoeff();
	Eigen::VectorXd input_bbox_max = _sample_points.rowwise().maxCoeff();
	Eigen::VectorXd input_bbox_size = input_bbox_max - input_bbox_min;

	_bbox_center = 0.5 * (input_bbox_min + input_bbox_max);
	_xy_size = std::sqrt(input_bbox_size[0] * input_bbox_size[0] + input_bbox_size[1] * input_bbox_size[1]);
	_z_size = input_bbox_size[2];
}

void get_transformed_sample_points(const MeshCuboidStructure &_cuboid_structure,
	const Real _xy_size, const Real _z_size, const Real _angle, Eigen::MatrixXd &_transformed_sample_points)
{
	Eigen::MatrixXd sample_points;
	Eigen::VectorXd bbox_center;
	Real xy_size, z_size;
	get_bounding_cylinder(_cuboid_structure, sample_points, bbox_center, xy_size, z_size);

	_transformed_sample_points = sample_points.colwise() - bbox_center;
	for (SamplePointIndex sample_point_index = 0; sample_point_index < _transformed_sample_points.cols();
		++sample_point_index)
	{
		_transformed_sample_points.col(sample_point_index)[0] *= (_xy_size / xy_size);
		_transformed_sample_points.col(sample_point_index)[1] *= (_xy_size / xy_size);
		_transformed_sample_points.col(sample_point_index)[2] *= (_z_size / z_size);
	}

	if (_angle != 0)
	{
		Eigen::AngleAxisd axis_rotation(_angle, Eigen::Vector3d::UnitZ());
		_transformed_sample_points = axis_rotation.toRotationMatrix() * _transformed_sample_points;
	}

	_transformed_sample_points = _transformed_sample_points.colwise() + bbox_center;
}

void MeshViewerCore::run_part_assembly_align_database(const std::string _mesh_filepath,
	Real &_xy_size, Real &_z_size, Real &_angle)
{
	Eigen::MatrixXd input_points;
	Eigen::VectorXd input_bbox_center;
	get_bounding_cylinder(cuboid_structure_, input_points, input_bbox_center,
		_xy_size, _z_size);

	ANNpointArray input_ann_points;
	ANNkd_tree* input_ann_kd_tree = ICP::create_kd_tree(input_points, input_ann_points);


	const unsigned int num_angles = 360;
	Real angle_unit = 2 * M_PI / static_cast<Real>(num_angles);
	std::vector<Real> angle_scores(num_angles, 0);


	MyMesh example_mesh;
	MeshCuboidStructure example_cuboid_structure(&example_mesh);


	// For every file in the base path.
	QDir input_dir((FLAGS_data_root_path + FLAGS_mesh_path).c_str());
	assert(input_dir.exists());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
	input_dir.setSorting(QDir::Name);


	QFileInfoList dir_list = input_dir.entryInfoList();
	for (int i = 0; i < dir_list.size(); i++)
	{
		QFileInfo example_file_info = dir_list.at(i);
		if (example_file_info.exists() &&
			(example_file_info.suffix().compare("obj") == 0
			|| example_file_info.suffix().compare("off") == 0))
		{
			std::string example_mesh_filepath = std::string(example_file_info.filePath().toLocal8Bit());

			// Skip if the mesh is the same with the input mesh.
			if (example_mesh_filepath.compare(_mesh_filepath) == 0)
				continue;

			std::string other_mesh_name(example_file_info.baseName().toLocal8Bit());

			bool ret = load_object_info(example_mesh, example_cuboid_structure,
				example_mesh_filepath.c_str(), LoadSamplePoints, NULL, false);
			if (!ret) continue;


			Eigen::MatrixXd scaled_example_points;
			get_transformed_sample_points(example_cuboid_structure, _xy_size, _z_size, 0, scaled_example_points);

			for (unsigned int angle_index = 0; angle_index < num_angles; ++angle_index)
			{
				Eigen::AngleAxisd axis_rotation(angle_index * angle_unit, Eigen::Vector3d::UnitZ());
				Eigen::MatrixXd rotated_example_points = axis_rotation.toRotationMatrix() * scaled_example_points;

				ANNpointArray rotated_example_ann_points;
				ANNkd_tree* rotated_example_ann_kd_tree = ICP::create_kd_tree(rotated_example_points,
					rotated_example_ann_points);

				Eigen::VectorXd input_to_rotated_example_distances;
				ICP::get_closest_points(rotated_example_ann_kd_tree, input_points, input_to_rotated_example_distances);

				Eigen::VectorXd rotated_example_to_input_distances;
				ICP::get_closest_points(input_ann_kd_tree, rotated_example_points, rotated_example_to_input_distances);

				Real score = std::max(input_to_rotated_example_distances.maxCoeff(),
					rotated_example_to_input_distances.maxCoeff());
				angle_scores[angle_index] += score;

				annDeallocPts(rotated_example_ann_points);
				delete rotated_example_ann_kd_tree;
			}
		}
	}

	annDeallocPts(input_ann_points);
	delete input_ann_kd_tree;

	Real min_score = std::numeric_limits<Real>::max();
	unsigned int min_angle_index = 0;
	for (unsigned int angle_index = 0; angle_index < num_angles; ++angle_index)
	{
		if (angle_scores[angle_index] < min_score)
		{
			min_angle_index = angle_index;
			min_score = angle_scores[angle_index];
		}
	}

	std::cout << "[" << min_angle_index << "]: " << min_score << std::endl;
	_angle = min_angle_index * angle_unit;
}

void MeshViewerCore::run_part_assembly_render_alignment(const std::string _mesh_filepath,
	const Real _xy_size, const Real _z_size, const Real _angle, const std::string _output_filename)
{
	MyMesh example_mesh;
	MeshCuboidStructure example_cuboid_structure(&example_mesh);

	// Load the same object.
	bool ret = load_object_info(example_mesh, example_cuboid_structure,
		_mesh_filepath.c_str(), LoadDenseSamplePoints, NULL, false);
	assert(ret);

	// Render partial scanned input points with error = 0 (for coloring).
	Eigen::MatrixXd aligned_example_points;
	get_transformed_sample_points(example_cuboid_structure, _xy_size, _z_size, _angle, aligned_example_points);

	for (SamplePointIndex sample_point_index = 0; sample_point_index < cuboid_structure_.num_sample_points();
		++sample_point_index)
	{
		MeshSamplePoint *sample_point = cuboid_structure_.sample_points_[sample_point_index];
		sample_point->error_ = 0.0;
	}

	// Render transformed original points with error = 1 (for coloring).
	assert(aligned_example_points.cols() == example_cuboid_structure.num_sample_points());
	for (SamplePointIndex sample_point_index = 0; sample_point_index < example_cuboid_structure.num_sample_points();
		++sample_point_index)
	{
		MeshSamplePoint *sample_point = example_cuboid_structure.sample_points_[sample_point_index];
		MyMesh::Normal normal = sample_point->normal_;
		MyMesh::Point point;
		for (int i = 0; i < 3; ++i)
			point[i] = aligned_example_points.col(sample_point_index)[i];
		MeshSamplePoint *new_sample_point = cuboid_structure_.add_sample_point(point, normal);
		new_sample_point->error_ = 1.0;
	}

	setDrawMode(COLORED_POINT_SAMPLES);
	updateGL();
	snapshot(_output_filename);

	setDrawMode(CUSTOM_VIEW);
}

void get_consistent_matching_parts(
	const MeshCuboidStructure &_cuboid_structure,
	const MeshCuboidTrainer &_trainer,
	std::vector< std::pair<Real, std::string> > &_label_matched_object_scores,
	std::vector<std::string> &_label_matched_objects)
{
	unsigned int num_labels = _cuboid_structure.num_labels();
	assert(_label_matched_object_scores.size() == num_labels);

	_label_matched_objects.clear();
	_label_matched_objects.resize(num_labels, "");


	// Remove conflicted labels.
	std::vector< std::list<LabelIndex> > conflicted_labels;
	_trainer.get_conflicted_labels(conflicted_labels);
	assert(conflicted_labels.size() == num_labels);

	bool *is_label_index_visited = new bool[num_labels];
	memset(is_label_index_visited, false, num_labels * sizeof(bool));

	while (true)
	{
		Real max_score = 0;
		LabelIndex max_score_label_index = 0;

		for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
		{
			Real score = _label_matched_object_scores[label_index].first;
			if (!is_label_index_visited[label_index] && score > max_score)
			{
				max_score = _label_matched_object_scores[label_index].first;
				max_score_label_index = label_index;
			}
		}

		if (max_score == 0)
			break;

		std::string other_mesh_filepath = _label_matched_object_scores[max_score_label_index].second;
		_label_matched_objects[max_score_label_index] = other_mesh_filepath;
		is_label_index_visited[max_score_label_index] = true;

		for (std::list<LabelIndex>::iterator it = conflicted_labels[max_score_label_index].begin();
			it != conflicted_labels[max_score_label_index].end(); ++it)
		{
			assert(*it < num_labels);
			is_label_index_visited[*it] = true;
		}
	}

	delete[] is_label_index_visited;


	// NOTE:
	// Select the same 3D model for symmetric parts.
	for (std::vector< MeshCuboidSymmetryGroupInfo >::const_iterator it = _cuboid_structure.symmetry_group_info_.begin();
		it != _cuboid_structure.symmetry_group_info_.end(); ++it)
	{
		const MeshCuboidSymmetryGroupInfo &symmetry_group = (*it);
		for (std::vector< std::pair<LabelIndex, LabelIndex> >::const_iterator jt = symmetry_group.pair_label_indices_.begin();
			jt != symmetry_group.pair_label_indices_.end(); ++jt)
		{
			LabelIndex label_index_1 = (*jt).first;
			LabelIndex label_index_2 = (*jt).second;

			// If either one of symmetric cuboids exists.
			if (_label_matched_objects[label_index_1] != ""
				|| _label_matched_objects[label_index_2] != "")
			{
				if (_label_matched_object_scores[label_index_1].first >= _label_matched_object_scores[label_index_2].first)
					_label_matched_objects[label_index_2] = _label_matched_objects[label_index_1];
				else
					_label_matched_objects[label_index_1] = _label_matched_objects[label_index_2];
			}
		}
	}
}

void MeshViewerCore::run_part_assembly_match_parts(const std::string _mesh_filepath,
	const Real _xy_size, const Real _z_size, const Real _angle,
	const MeshCuboidTrainer &_trainer, std::vector<std::string> &_label_matched_objects)
{
	// Parameters.
	const Real part_assembly_window_size = FLAGS_param_part_assembly_window_size *
		cuboid_structure_.mesh_->get_object_diameter();
	const Real part_assembly_voxel_size = FLAGS_param_part_assembly_voxel_size *
		cuboid_structure_.mesh_->get_object_diameter();
	const Real part_assembly_voxel_variance = FLAGS_param_part_assembly_voxel_variance *
		cuboid_structure_.mesh_->get_object_diameter();
	const Real distance_param = 2 * part_assembly_voxel_variance * part_assembly_voxel_variance;
	assert(distance_param > 0);


	// Create input point KD tree.
	unsigned int num_labels = cuboid_structure_.num_labels();

	Eigen::MatrixXd input_sample_points(3, cuboid_structure_.num_sample_points());
	for (SamplePointIndex sample_point_index = 0; sample_point_index < cuboid_structure_.num_sample_points();
		++sample_point_index)
	{
		MeshSamplePoint *sample_point = cuboid_structure_.sample_points_[sample_point_index];
		assert(sample_point);
		for (unsigned int i = 0; i < 3; ++i)
			input_sample_points.col(sample_point_index)(i) = sample_point->point_[i];
	}

	ANNpointArray input_ann_points;
	ANNkd_tree* input_ann_kd_tree = ICP::create_kd_tree(input_sample_points, input_ann_points);


	// Load database.
	MyMesh example_mesh;
	MeshCuboidStructure example_cuboid_structure(&example_mesh);

	bool ret = true;
	ret = ret & example_cuboid_structure.load_labels((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_info_filename).c_str());
	ret = ret & example_cuboid_structure.load_label_symmetries((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_symmetry_info_filename).c_str());

	// Load symmetry groups.
	ret = ret & example_cuboid_structure.load_symmetry_groups((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_symmetry_group_info_filename).c_str());

	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open label information files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}


	// For every file in the base path.
	QDir input_dir((FLAGS_data_root_path + FLAGS_mesh_path).c_str());
	assert(input_dir.exists());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
	input_dir.setSorting(QDir::Name);


	// (Score, mesh_filepath)
	std::vector< std::pair<Real, std::string> > label_matched_object_scores(num_labels);
	for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
		label_matched_object_scores[label_index].first = 0;


	QFileInfoList dir_list = input_dir.entryInfoList();
	for (int i = 0; i < dir_list.size(); i++)
	{
		QFileInfo example_file_info = dir_list.at(i);
		if (example_file_info.exists() &&
			(example_file_info.suffix().compare("obj") == 0
			|| example_file_info.suffix().compare("off") == 0))
		{
			std::string example_mesh_filepath = std::string(example_file_info.filePath().toLocal8Bit());

			// Skip if the mesh is the same with the input mesh.
			if (example_mesh_filepath.compare(_mesh_filepath) == 0)
				continue;

			std::string example_mesh_name(example_file_info.baseName().toLocal8Bit());
			std::string example_cuboid_filepath = FLAGS_training_dir + std::string("/")
				+ example_mesh_name + std::string(".arff");

			bool ret = load_object_info(example_mesh, example_cuboid_structure,
				example_mesh_filepath.c_str(), LoadDenseSamplePoints, example_cuboid_filepath.c_str(), false);
			if (!ret) continue;

			std::cout << "mesh: " << example_mesh_name << std::endl;

			Eigen::MatrixXd transformed_example_points;
			get_transformed_sample_points(example_cuboid_structure, _xy_size, _z_size, _angle, transformed_example_points);

			// Assign transformed points.
			for (SamplePointIndex sample_point_index = 0; sample_point_index < example_cuboid_structure.num_sample_points();
				++sample_point_index)
			{
				MeshSamplePoint *sample_point = example_cuboid_structure.sample_points_[sample_point_index];
				for (int i = 0; i < 3; ++i)
					sample_point->point_[i] = transformed_example_points.col(sample_point_index)[i];
			}

			for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
			{
				MeshCuboid *example_cuboid = NULL;

				// NOTE:
				// The current implementation assumes that there is only one part for each label.
				assert(example_cuboid_structure.label_cuboids_[label_index].size() <= 1);
				if (!example_cuboid_structure.label_cuboids_[label_index].empty())
					example_cuboid = example_cuboid_structure.label_cuboids_[label_index].front();

				if (!example_cuboid) continue;
				else if (example_cuboid->num_sample_points() == 0) continue;

				// Create cuboid KD-tree.
				std::vector<MyMesh::Point> example_sample_points;
				example_cuboid->get_sample_points(example_sample_points);

				Real max_bbox_size = 0;
				for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
					max_bbox_size = std::max(max_bbox_size, example_cuboid->get_bbox_size()[axis_index]);
				assert(max_bbox_size > 0);

				MyMesh::Point local_coord_bbox_min = example_cuboid->get_bbox_center() - MyMesh::Point(max_bbox_size);
				MyMesh::Point local_coord_bbox_max = example_cuboid->get_bbox_center() + MyMesh::Point(max_bbox_size);
				
				MeshCuboidVoxelGrid local_coord_voxels(local_coord_bbox_min, local_coord_bbox_max,
					part_assembly_voxel_size);
				int num_voxels = local_coord_voxels.n_voxels();


				// Compute distance maps.
				Eigen::VectorXd input_distance_map;
				local_coord_voxels.get_distance_map(input_ann_points, input_ann_kd_tree, input_distance_map);
				assert(input_distance_map.rows() == num_voxels);
				for (unsigned int i = 0; i < num_voxels; ++i)
					input_distance_map[i] = std::exp(-input_distance_map[i] * input_distance_map[i] / distance_param);

				Eigen::VectorXd example_voxel_occupancies;
				local_coord_voxels.get_voxel_occupancies(example_sample_points, example_voxel_occupancies);
				assert(example_voxel_occupancies.rows() == num_voxels);

				int num_input_occupied_voxels = (input_distance_map.array() > 0.95).count();
				int num_example_occupied_voxels = (example_voxel_occupancies.array() >= 1).count();
				if (num_input_occupied_voxels == 0 || num_example_occupied_voxels == 0)
					continue;

				// Equation (1) ~ (3).
				Real score = input_distance_map.dot(example_voxel_occupancies);
				assert(score >= 0);
				score = (1 - 0.7) * (score / num_example_occupied_voxels)
					+ (0.7) * (score / num_input_occupied_voxels);

				// DEBUG.
				printf("[%s] (%d): %lf\n", example_mesh_name.c_str(), label_index, score);

				if (score > label_matched_object_scores[label_index].first)
				{
					label_matched_object_scores[label_index].first = score;
					label_matched_object_scores[label_index].second = example_mesh_filepath;
				}
			}
		}
	}

	annDeallocPts(input_ann_points);
	delete input_ann_kd_tree;

	get_consistent_matching_parts(cuboid_structure_, _trainer,
		label_matched_object_scores, _label_matched_objects);
}

void MeshViewerCore::run_part_assembly_reconstruction(const std::string _mesh_filepath,
	const Real _xy_size, const Real _z_size, const Real _angle,
	const std::vector<std::string> &_label_matched_objects)
{
	unsigned int num_labels = cuboid_structure_.num_labels();
	assert(_label_matched_objects.size() == num_labels);

	MyMesh example_mesh;
	MeshCuboidStructure example_cuboid_structure(&example_mesh);

	bool ret = example_cuboid_structure.load_labels((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_info_filename).c_str());

	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open label information files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}


	cuboid_structure_.clear_sample_points();
	cuboid_structure_.clear_cuboids();

	for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
	{
		std::string mesh_filepath = _label_matched_objects[label_index];
		if (mesh_filepath == "")
			continue;

		QFileInfo file_info(mesh_filepath.c_str());
		std::string mesh_name(file_info.baseName().toLocal8Bit());
		std::string cuboid_filepath = FLAGS_training_dir + std::string("/") + mesh_name + std::string(".arff");

		std::cout << "--------" << std::endl;
		std::cout << "Label (" << label_index << "):" << std::endl;
		std::cout << "Mesh name: " << mesh_name << std::endl;
		
		bool ret = load_object_info(example_mesh, example_cuboid_structure,
			mesh_filepath.c_str(), LoadDenseSamplePoints, cuboid_filepath.c_str());
		assert(ret);

		Eigen::MatrixXd transformed_example_points;
		get_transformed_sample_points(example_cuboid_structure, _xy_size, _z_size, _angle, transformed_example_points);

		// Assign transformed points.
		for (SamplePointIndex sample_point_index = 0; sample_point_index < example_cuboid_structure.num_sample_points();
			++sample_point_index)
		{
			MeshSamplePoint *sample_point = example_cuboid_structure.sample_points_[sample_point_index];
			for (int i = 0; i < 3; ++i)
				sample_point->point_[i] = transformed_example_points.col(sample_point_index)[i];
		}


		assert(example_cuboid_structure.label_cuboids_[label_index].size() <= 1);
		MeshCuboid *example_cuboid = NULL;
		if (!example_cuboid_structure.label_cuboids_[label_index].empty())
			example_cuboid = example_cuboid_structure.label_cuboids_[label_index].front();
		assert(example_cuboid);

		const int num_points = example_cuboid->num_sample_points();
		std::cout << "# of sample points: " << num_points << std::endl;
		std::cout << "Copying... ";

		MeshCuboid *cuboid = new MeshCuboid(label_index);
		cuboid_structure_.label_cuboids_[label_index].push_back(cuboid);

		for (int point_index = 0; point_index < num_points; ++point_index)
		{
			const MeshSamplePoint* sample_point = example_cuboid->get_sample_point(point_index);

			MyMesh::Point point = sample_point->point_;
			MyMesh::Normal normal = sample_point->normal_;
			MeshSamplePoint *new_sample_point = cuboid_structure_.add_sample_point(point, normal);

			new_sample_point->label_index_confidence_.resize(num_labels, 0.0);
			new_sample_point->label_index_confidence_[label_index] = 1.0;

			cuboid->add_sample_point(new_sample_point);
		}
	}
}

void MeshViewerCore::run_part_assembly()
{
	// Load basic information.
	bool ret = true;

	ret = ret & cuboid_structure_.load_labels((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_info_filename).c_str());
	ret = ret & cuboid_structure_.load_label_symmetries((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_symmetry_info_filename).c_str());

	// Load symmetry groups.
	ret = ret & cuboid_structure_.load_symmetry_groups((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_symmetry_group_info_filename).c_str());

	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open label information files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}


	ret = true;
	MeshCuboidTrainer trainer;
	ret = ret & trainer.load_object_list(FLAGS_training_dir + std::string("/") + FLAGS_object_list_filename);
	ret = ret & trainer.load_features(FLAGS_training_dir + std::string("/") + FLAGS_feature_filename_prefix);
	ret = ret & trainer.load_transformations(FLAGS_training_dir + std::string("/") + FLAGS_transformation_filename_prefix);

	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open training files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}


	// Check file paths.
	setDrawMode(CUSTOM_VIEW);

	std::string mesh_filepath = FLAGS_data_root_path + FLAGS_mesh_path + std::string("/") + FLAGS_mesh_filename;
	QFileInfo mesh_file(mesh_filepath.c_str());
	std::string mesh_name = std::string(mesh_file.baseName().toLocal8Bit());

	if (!mesh_file.exists())
	{
		std::cerr << "Error: The mesh file does not exist (" << mesh_filepath << ")." << std::endl;
		return;
	}
	else if (!open_mesh(mesh_filepath.c_str()))
	{
		std::cerr << "Error: The mesh file cannot be opened (" << mesh_filepath << ")." << std::endl;
		return;
	}

	std::string filename_prefix = std::string("/") + mesh_name + std::string("_");
	std::stringstream output_filename_sstr;

	QDir output_dir;
	std::string mesh_output_path = FLAGS_part_assembly_dir + std::string("/") + mesh_name;
	output_dir.mkpath(FLAGS_part_assembly_dir.c_str());
	output_dir.mkpath(mesh_output_path.c_str());

	std::list<std::string> ignored_object_list;
	ignored_object_list.push_back(mesh_name);

	std::vector< std::vector<MeshCuboidJointNormalRelations *> > joint_normal_relations;
	//MeshCuboidTrainer::load_joint_normal_relations(num_labels, "joint_normal_", joint_normal_relations);
	trainer.get_joint_normal_relations(joint_normal_relations, &ignored_object_list);
	MeshCuboidJointNormalRelationPredictor joint_normal_predictor(joint_normal_relations);


	// Initialize basic information.
	unsigned int num_labels = cuboid_structure_.num_labels();
	cuboid_structure_.clear_cuboids();
	cuboid_structure_.clear_sample_points();


	// Load files.
	double snapshot_modelview_matrix[16];
	double occlusion_modelview_matrix[16];

	open_modelview_matrix_file(FLAGS_pose_filename.c_str());
	memcpy(snapshot_modelview_matrix, modelview_matrix(), 16 * sizeof(double));

	if (FLAGS_occlusion_pose_filename == "")
		set_random_view_direction(true);
	else
	{
		ret = open_modelview_matrix_file(FLAGS_occlusion_pose_filename.c_str());
		assert(ret);
	}
	memcpy(occlusion_modelview_matrix, modelview_matrix(), 16 * sizeof(double));
	//save_modelview_matrix_file((mesh_output_path + std::string("/occlusion_pose.txt")).c_str());


	ret = load_object_info(mesh_, cuboid_structure_, mesh_filepath.c_str(), LoadSamplePoints);
	if (!ret) return;
	set_modelview_matrix(occlusion_modelview_matrix, false);
	remove_occluded_points();
	set_modelview_matrix(snapshot_modelview_matrix);


	std::cout << "Align database... " << std::endl;
	Real xy_size, z_size, angle;
	run_part_assembly_align_database(mesh_filepath, xy_size, z_size, angle);

	std::cout << "Match parts... " << std::endl;
	std::vector<std::string> label_matched_objects;
	run_part_assembly_match_parts(mesh_filepath, xy_size, z_size, angle,
		trainer, label_matched_objects);


	ret = load_object_info(mesh_, cuboid_structure_, mesh_filepath.c_str(), LoadDenseSamplePoints);
	if (!ret) return;
	set_modelview_matrix(occlusion_modelview_matrix, false);
	remove_occluded_points();

	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << mesh_output_path << filename_prefix << std::string("view");
	snapshot(output_filename_sstr.str().c_str());
	set_modelview_matrix(snapshot_modelview_matrix);
	updateGL();


	std::cout << "Render alignment... " << std::endl;
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << mesh_output_path << filename_prefix << std::string("assembly");
	run_part_assembly_render_alignment(mesh_filepath, xy_size, z_size, angle, output_filename_sstr.str());


	std::cout << "Assemble parts... " << std::endl;
	run_part_assembly_reconstruction(mesh_filepath, xy_size, z_size, angle, label_matched_objects);


	// Evaluation.
	MyMesh ground_truth_mesh;
	MeshCuboidStructure ground_truth_cuboid_structure(&ground_truth_mesh);

	ret = true;
	ret = ret & ground_truth_cuboid_structure.load_labels((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_info_filename).c_str());
	ret = ret & ground_truth_cuboid_structure.load_label_symmetries((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_symmetry_info_filename).c_str());

	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open label information files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	ret = load_object_info(ground_truth_mesh, ground_truth_cuboid_structure,
		mesh_filepath.c_str(), LoadGroundTruthData);
	assert(ret);
	MeshCuboidEvaluator evaluator(&ground_truth_cuboid_structure);

	evaluator.evaluate_point_to_point_distances(&cuboid_structure_, output_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points_to_ply(output_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points((output_filename_sstr.str() + std::string(".pts")).c_str());
	cuboid_structure_.save_sample_point_labels((output_filename_sstr.str() + std::string("_label.arff")).c_str());

	setDrawMode(COLORED_POINT_SAMPLES);

	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << mesh_output_path << filename_prefix << std::string("assembly_accuracy");
	snapshot(output_filename_sstr.str().c_str());

	cuboid_structure_ = ground_truth_cuboid_structure;
	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << mesh_output_path << filename_prefix << std::string("assembly_completeness");
	snapshot(output_filename_sstr.str().c_str()); updateGL();
	snapshot(output_filename_sstr.str().c_str());

	setDrawMode(CUSTOM_VIEW);
}
