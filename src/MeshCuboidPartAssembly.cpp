#include "MeshViewerCore.h"
#include "MeshCuboidParameters.h"

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
	const Real _scale_xy_size, const Real _scale_z_size, const Real _angle, Eigen::MatrixXd &_transformed_sample_points)
{
	Eigen::MatrixXd sample_points;
	Eigen::VectorXd bbox_center;
	Real xy_size, z_size;
	get_bounding_cylinder(_cuboid_structure, sample_points, bbox_center, xy_size, z_size);

	_transformed_sample_points = sample_points.colwise() - bbox_center;
	for (SamplePointIndex sample_point_index = 0; sample_point_index < _transformed_sample_points.cols();
		++sample_point_index)
	{
		_transformed_sample_points.col(sample_point_index)[0] *= (_scale_xy_size / xy_size);
		_transformed_sample_points.col(sample_point_index)[1] *= (_scale_xy_size / xy_size);
		_transformed_sample_points.col(sample_point_index)[2] *= (_scale_z_size / z_size);
	}
	_transformed_sample_points = _transformed_sample_points.colwise() + bbox_center;

	if (_angle != 0)
	{
		Eigen::AngleAxisd axis_rotation(_angle, Eigen::Vector3d::UnitZ());
		_transformed_sample_points = axis_rotation.toRotationMatrix() * _transformed_sample_points;
	}
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


	const unsigned int num_angles = 180;
	Real angle_unit = 2 * M_PI / static_cast<Real>(num_angles);
	std::vector<Real> angle_scores(180, 0);


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
		QFileInfo other_file_info = dir_list.at(i);
		if (other_file_info.exists() &&
			(other_file_info.suffix().compare("obj") == 0
			|| other_file_info.suffix().compare("off") == 0))
		{
			std::string other_mesh_filepath = std::string(other_file_info.filePath().toLocal8Bit());

			// Skip if the mesh is the same with the input mesh.
			if (other_mesh_filepath.compare(_mesh_filepath) == 0)
				continue;

			std::string other_mesh_name(other_file_info.baseName().toLocal8Bit());

			bool ret = load_object_info(example_mesh, example_cuboid_structure,
				other_mesh_filepath.c_str(), LoadSamplePoints);
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

				Real score = (input_to_rotated_example_distances.sum() + rotated_example_to_input_distances.sum())
					/ (input_to_rotated_example_distances.rows() + rotated_example_to_input_distances.rows());
				angle_scores[angle_index] += score;

				annDeallocPts(rotated_example_ann_points);
				delete rotated_example_ann_kd_tree;
			}
		}
	}

	annDeallocPts(input_ann_points); delete input_ann_kd_tree;

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

	bool ret = load_object_info(example_mesh, example_cuboid_structure,
		_mesh_filepath.c_str(), LoadDenseSamplePoints);
	assert(ret);

	Eigen::MatrixXd aligned_example_points;
	get_transformed_sample_points(example_cuboid_structure, _xy_size, _z_size, _angle, aligned_example_points);

	for (SamplePointIndex sample_point_index = 0; sample_point_index < cuboid_structure_.num_sample_points();
		++sample_point_index)
	{
		MeshSamplePoint *sample_point = cuboid_structure_.sample_points_[sample_point_index];
		sample_point->error_ = 0.0;
	}

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

void MeshViewerCore::run_part_assembly_match_parts(const std::string _mesh_filepath, const Real _angle)
{
	// ---- //
	unsigned int num_input_points = cuboid_structure_.sample_points_.size();
	assert(num_input_points > 0);
	Eigen::MatrixXd input_points(3, num_input_points);

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_input_points;
		++sample_point_index)
	{
		MeshSamplePoint *sample_point = cuboid_structure_.sample_points_[sample_point_index];
		for (int i = 0; i < 3; ++i)
			input_points.col(sample_point_index)[i] = sample_point->point_[i];
	}

	Eigen::VectorXd input_bbox_min = input_points.rowwise().minCoeff();
	Eigen::VectorXd input_bbox_max = input_points.rowwise().maxCoeff();
	Eigen::VectorXd input_bbox_size = input_bbox_max - input_bbox_min;
	Eigen::VectorXd input_bbox_center = 0.5 * (input_bbox_min + input_bbox_max);

	Real input_z_size = input_bbox_size[2];
	Real input_xy_size = std::sqrt(input_bbox_size[0] * input_bbox_size[0]
		+ input_bbox_size[1] * input_bbox_size[1]);

	ANNpointArray input_ann_points;
	ANNkd_tree* input_ann_kd_tree = ICP::create_kd_tree(input_points, input_ann_points);


	const unsigned int num_angles = 180;
	Real angle_unit = 2 * M_PI / static_cast<Real>(num_angles);
	std::vector<Real> angle_scores(180, 0);
	// ---- //


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


	QFileInfoList dir_list = input_dir.entryInfoList();
	for (int i = 0; i < dir_list.size(); i++)
	{
		QFileInfo other_file_info = dir_list.at(i);
		if (other_file_info.exists() &&
			(other_file_info.suffix().compare("obj") == 0
			|| other_file_info.suffix().compare("off") == 0))
		{
			std::string other_mesh_filepath = std::string(other_file_info.filePath().toLocal8Bit());

			// Skip if the mesh is the same with the input mesh.
			if (other_mesh_filepath.compare(_mesh_filepath) == 0)
				continue;

			std::string other_mesh_name(other_file_info.baseName().toLocal8Bit());

			bool ret = load_object_info(example_mesh, example_cuboid_structure,
				other_mesh_filepath.c_str(), LoadDenseTestData);
			if (!ret) continue;

			// ---- //
			unsigned int num_example_points = example_cuboid_structure.sample_points_.size();
			assert(num_example_points > 0);
			Eigen::MatrixXd example_points(3, num_example_points);

			for (SamplePointIndex sample_point_index = 0; sample_point_index < num_example_points;
				++sample_point_index)
			{
				MeshSamplePoint *sample_point = example_cuboid_structure.sample_points_[sample_point_index];
				for (int i = 0; i < 3; ++i)
					example_points.col(sample_point_index)[i] = sample_point->point_[i];
			}

			Eigen::VectorXd example_bbox_min = example_points.rowwise().minCoeff();
			Eigen::VectorXd example_bbox_max = example_points.rowwise().maxCoeff();
			Eigen::VectorXd example_bbox_size = example_bbox_max - example_bbox_min;
			Eigen::VectorXd example_bbox_center = 0.5 * (example_bbox_min + example_bbox_max);

			Real example_z_size = example_bbox_size[2];
			Real example_xy_size = std::sqrt(example_bbox_size[0] * example_bbox_size[0]
				+ example_bbox_size[1] * example_bbox_size[1]);

			Eigen::MatrixXd transformed_example_points = example_points.colwise() - example_bbox_center;
			for (SamplePointIndex sample_point_index = 0; sample_point_index < num_example_points;
				++sample_point_index)
			{
				transformed_example_points.col(sample_point_index)[0] *= (input_xy_size / example_xy_size);
				transformed_example_points.col(sample_point_index)[1] *= (input_xy_size / example_xy_size);
				transformed_example_points.col(sample_point_index)[2] *= (input_z_size / example_z_size);
			}
			transformed_example_points = transformed_example_points.colwise() + example_bbox_center;



			Eigen::AngleAxisd axis_rotation(_angle, Eigen::Vector3d::UnitZ());
			Eigen::MatrixXd rotated_example_points = axis_rotation.toRotationMatrix() * transformed_example_points;

			ANNpointArray example_ann_points;
			ANNkd_tree* example_ann_kd_tree = ICP::create_kd_tree(rotated_example_points, example_ann_points);

			Eigen::VectorXd input_to_example_distances;
			ICP::get_closest_points(example_ann_kd_tree, input_points, input_to_example_distances);

			Eigen::VectorXd example_to_input_distances;
			ICP::get_closest_points(input_ann_kd_tree, example_points, example_to_input_distances);

			annDeallocPts(example_ann_points); delete example_ann_kd_tree;
			// ---- //
		}
	}

	annDeallocPts(input_ann_points); delete input_ann_kd_tree;

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
	std::string part_assembly_output_path = FLAGS_output_dir + std::string("/part_assembly/");
	output_dir.mkpath(FLAGS_output_dir.c_str());
	output_dir.mkpath(part_assembly_output_path.c_str());


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

	Real xy_size, z_size, angle;
	run_part_assembly_align_database(mesh_filepath, xy_size, z_size, angle);


	ret = load_object_info(mesh_, cuboid_structure_, mesh_filepath.c_str(), LoadDenseSamplePoints);
	if (!ret) return;
	set_modelview_matrix(occlusion_modelview_matrix, false);
	remove_occluded_points();
	set_modelview_matrix(snapshot_modelview_matrix);

	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << part_assembly_output_path << filename_prefix << std::string("align");
	run_part_assembly_render_alignment(mesh_filepath, xy_size, z_size, angle, output_filename_sstr.str());

	//run_part_assembly_match_parts(mesh_filepath, angle);
}