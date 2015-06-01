#include "MeshViewerCore.h"
#include "MeshCuboidEvaluator.h"
#include "MeshCuboidFusion.h"
#include "MeshCuboidParameters.h"
#include "MeshCuboidPredictor.h"
#include "MeshCuboidRelation.h"
#include "MeshCuboidTrainer.h"
#include "MeshCuboidSolver.h"
#include "simplerandom.h"
#include "SymmetryDetection.h"
//#include "QGLOcculsionTestWidget.h"

#include <sstream>
#include <Eigen/Core>
#include <gflags/gflags.h>
#include <QDir>
#include <QFileInfo>
#include "MRFEnergy.h"


void MeshViewerCore::parse_arguments()
{
	std::cout << "data_root_path = " << FLAGS_data_root_path << std::endl;
	std::cout << "label_info_path = " << FLAGS_data_root_path + FLAGS_label_info_path << std::endl;
	std::cout << "mesh_path = " << FLAGS_data_root_path + FLAGS_mesh_path << std::endl;
	std::cout << "sample_path = " << FLAGS_data_root_path + FLAGS_sample_path << std::endl;
	std::cout << "sample_label_path = " << FLAGS_data_root_path + FLAGS_sample_label_path << std::endl;
	std::cout << "mesh_label_path = " << FLAGS_data_root_path + FLAGS_mesh_label_path << std::endl;
	std::cout << "output_dir = " << FLAGS_output_dir << std::endl;

	if (FLAGS_run_ground_truth_cuboids)
	{
		std::cout << "mesh_filename = " << FLAGS_mesh_filename << std::endl;
		compute_ground_truth_cuboids();
		exit(EXIT_FAILURE);
	}
	else if (FLAGS_run_training)
	{
		train();
		exit(EXIT_FAILURE);
	}
	else if (FLAGS_run_prediction)
	{
		std::cout << "mesh_filename = " << FLAGS_mesh_filename << std::endl;
		predict();
		exit(EXIT_FAILURE);
	}
	else if (FLAGS_run_part_assembly)
	{
		std::cout << "mesh_filename = " << FLAGS_mesh_filename << std::endl;
		run_part_assembly();
		exit(EXIT_FAILURE);
	}
	else if (FLAGS_run_symmetry_detection)
	{
		std::cout << "mesh_filename = " << FLAGS_mesh_filename << std::endl;
		run_symmetry_detection_msh2pln();
		exit(EXIT_FAILURE);
	}
	else if (FLAGS_run_baseline)
	{
		std::cout << "mesh_filename = " << FLAGS_mesh_filename << std::endl;
		run_baseline_stats();
		exit(EXIT_FAILURE);
	}
	else if (FLAGS_run_render_assembly)
	{
		std::cout << "mesh_filename = " << FLAGS_mesh_filename << std::endl;
		render_part_assembly_cuboids();
		exit(EXIT_FAILURE);
	}
	else if (FLAGS_run_render_evaluation)
	{
		std::cout << "mesh_filename = " << FLAGS_mesh_filename << std::endl;
		run_render_evaluation();
		exit(EXIT_FAILURE);
	}
	else if (FLAGS_run_extract_symmetry_info)
	{
		std::cout << "mesh_filename = " << FLAGS_mesh_filename << std::endl;
		run_extract_symmetry_info();
		exit(EXIT_FAILURE);
	}
}

bool MeshViewerCore::load_object_info(
	MyMesh &_mesh,
	MeshCuboidStructure &_cuboid_structure,
	const char* _mesh_filepath,
	const LoadObjectInfoOption _option,
	const char* _cuboid_filepath,
	bool _verbose)
{
	bool ret;
	QFileInfo file_info(_mesh_filepath);

	std::string mesh_name(file_info.baseName().toLocal8Bit());
	std::string mesh_label_filepath = FLAGS_data_root_path + FLAGS_mesh_label_path + std::string("/") + mesh_name + std::string(".seg");
	std::string sample_filepath = FLAGS_data_root_path + FLAGS_sample_path + std::string("/") + mesh_name + std::string(".pts");
	std::string sample_label_filepath = FLAGS_data_root_path + FLAGS_sample_label_path + std::string("/") + mesh_name + std::string(".arff");
	std::string dense_sample_filepath = FLAGS_data_root_path + FLAGS_dense_sample_path + std::string("/") + mesh_name + std::string(".pts");

	QFileInfo mesh_file(_mesh_filepath);
	QFileInfo mesh_label_file(mesh_label_filepath.c_str());
	QFileInfo sample_file(sample_filepath.c_str());
	QFileInfo sample_label_file(sample_label_filepath.c_str());
	QFileInfo dense_sample_file(dense_sample_filepath.c_str());
	
	if (!mesh_file.exists())
	{
		std::cerr << "Error: The mesh file does not exist (" << _mesh_filepath << ")." << std::endl;
		return false;
	}
	if (!(_option == LoadTestData || _option == LoadDenseTestData) && !mesh_label_file.exists())
	{
		std::cerr << "Error: The mesh label file does not exist (" << mesh_label_filepath << ")." << std::endl;
		return false;
	}
	if (!sample_file.exists())
	{
		std::cerr << "Error: The sample file does not exist (" << sample_filepath << ")." << std::endl;
		return false;
	}
	if ((_option == LoadTestData || _option == LoadDenseTestData) && !sample_label_file.exists())
	{
		std::cerr << "Error: The sample label file does not exist (" << sample_label_filepath << ")." << std::endl;
		return false;
	}
	if ((_option == LoadDenseSamplePoints || _option == LoadDenseTestData) && !dense_sample_file.exists())
	{
		std::cerr << "Error: The dense sample file does not exist (" << dense_sample_filepath << ")." << std::endl;
		return false;
	}

	_cuboid_structure.clear_cuboids();
	_cuboid_structure.clear_sample_points();

	if ((&_mesh) == (&mesh_))
	{
		ret = open_mesh(_mesh_filepath);
		assert(ret);
	}
	else
	{
		ret = _mesh.open_mesh(_mesh_filepath, _verbose);
		assert(ret);
	}

	if (_option == LoadSamplePoints)
	{
		if (_verbose) std::cout << " - Load mesh face labels." << std::endl;
		ret = _mesh.load_face_label_simple(mesh_label_filepath.c_str(), false);
		assert(ret);

		if (_verbose) std::cout << " - Load sample points." << std::endl;
		ret = _cuboid_structure.load_sample_points(sample_filepath.c_str(), false);
		assert(ret);

		_cuboid_structure.apply_mesh_face_labels_to_sample_points();
	}
	else if (_option == LoadDenseSamplePoints)
	{
		if (_verbose) std::cout << " - Load mesh face labels." << std::endl;
		ret = _mesh.load_face_label_simple(mesh_label_filepath.c_str(), false);
		assert(ret);

		// NOTE:
		// Load dense sample points as base sample points.
		// Do not use load_dense_sample_points() function,
		// which requires to load sparse sample point first.
		if (_verbose) std::cout << " - Load dense sample points." << std::endl;
		ret = _cuboid_structure.load_sample_points(dense_sample_filepath.c_str(), false);
		assert(ret);

		_cuboid_structure.apply_mesh_face_labels_to_sample_points();
	}
	else if (_option == LoadGroundTruthData)
	{
		if (_verbose) std::cout << " - Load mesh face labels." << std::endl;
		ret = _mesh.load_face_label_simple(mesh_label_filepath.c_str(), false);
		assert(ret);

		// NOTE:
		// Load dense sample points as base sample points.
		// Do not use load_dense_sample_points() function,
		// which requires to load sparse sample point first.
		if (_verbose) std::cout << " - Load dense sample points." << std::endl;
		ret = _cuboid_structure.load_sample_points(dense_sample_filepath.c_str(), false);
		assert(ret);

		if (_verbose) std::cout << " - Compute ground truth cuboids." << std::endl;
		_cuboid_structure.get_mesh_face_label_cuboids();
	}
	else if (_option == LoadTestData || _option == LoadDenseTestData)
	{
		if (_verbose) std::cout << " - Load sample points." << std::endl;
		ret = _cuboid_structure.load_sample_points(sample_filepath.c_str(), false);
		assert(ret);
		assert(_cuboid_structure.num_sample_points() > 0);

		if (_verbose) std::cout << " - Load sample point labels." << std::endl;
		ret = _cuboid_structure.load_sample_point_labels(sample_label_filepath.c_str());
		assert(ret);

		if (_option == LoadDenseTestData)
		{
			if (_verbose) std::cout << " - Load dense sample points." << std::endl;
			ret = _cuboid_structure.load_dense_sample_points(dense_sample_filepath.c_str(), false);
			assert(ret);
		}
	}

	if (_cuboid_filepath)
	{
		assert(_option != LoadGroundTruthData);

		ret = _cuboid_structure.load_cuboids(_cuboid_filepath, _verbose);
		if (!ret) return false;

		std::vector<MeshCuboid *> all_cuboids = _cuboid_structure.get_all_cuboids();
		for (std::vector<MeshCuboid *>::iterator it = all_cuboids.begin(); it != all_cuboids.end(); ++it)
		{
			MeshCuboid *cuboid = (*it);
			cuboid->clear_sample_points();
		}

		for (SamplePointIndex sample_point_index = 0; sample_point_index < _cuboid_structure.num_sample_points();
			++sample_point_index)
		{
			MeshSamplePoint* sample_point = _cuboid_structure.sample_points_[sample_point_index];
			assert(sample_point);

			FaceIndex fid = sample_point->corr_fid_;
			assert(fid < _cuboid_structure.mesh_->n_faces());

			MyMesh::FaceHandle fh = _cuboid_structure.mesh_->face_handle(fid);
			Label label = _cuboid_structure.mesh_->property(_cuboid_structure.mesh_->face_label_, fh);
			LabelIndex label_index;
			bool ret = _cuboid_structure.exist_label(label, &label_index);

			if (!ret || label_index >= _cuboid_structure.num_labels())
				continue;

			MeshCuboid *cuboid = NULL;
			// NOTE:
			// The current implementation assumes that there is only one part for each label.
			assert(_cuboid_structure.label_cuboids_[label_index].size() <= 1);
			if (!_cuboid_structure.label_cuboids_[label_index].empty())
				cuboid = _cuboid_structure.label_cuboids_[label_index].front();

			if (cuboid)
				cuboid->add_sample_point(sample_point);
		}
	}

	mesh_.clear_colors();

	return true;
}

void MeshViewerCore::open_cuboid_file(const char *_filename)
{
	assert(_filename);

	// Load basic information.
	bool ret = true;

	ret = ret & cuboid_structure_.load_labels((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_info_filename).c_str());
	ret = ret & cuboid_structure_.load_label_symmetries((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_symmetry_info_filename).c_str());

	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open label information files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	// To be removed.
	//if (FLAGS_use_symmetric_group_cuboids)
	//{
	//	cuboid_structure_.add_symmetric_group_labels();
	//}

	cuboid_structure_.load_cuboids(_filename);
	draw_cuboid_axes_ = true;
	setDrawMode(CUSTOM_VIEW);
	updateGL();
}

void MeshViewerCore::set_view_direction()
{
	Eigen::Vector3d view_direction;
	// -Z direction.
	view_direction << -modelview_matrix()[2], -modelview_matrix()[6], -modelview_matrix()[10];
	Eigen::Matrix3d rotation_mat;
	rotation_mat << modelview_matrix()[0], modelview_matrix()[4], modelview_matrix()[8],
		modelview_matrix()[1], modelview_matrix()[5], modelview_matrix()[9],
		modelview_matrix()[2], modelview_matrix()[6], modelview_matrix()[10];
	Eigen::Vector3d translation_vec;
	translation_vec << modelview_matrix()[12], modelview_matrix()[13], modelview_matrix()[14];
	Eigen::Vector3d view_point = -rotation_mat.transpose() * translation_vec;

	for (unsigned int i = 0; i < 3; ++i)
	{
		view_point_[i] = view_point[i];
		view_direction_[i] = view_direction[i];
	}

	update_cuboid_surface_points(cuboid_structure_, modelview_matrix());
	updateGL();
}

void MeshViewerCore::set_random_view_direction(bool _set_modelview_matrix)
{
	Eigen::Matrix4d centering_mat = Eigen::Matrix4d::Identity();
	for (int i = 0; i < 3; ++i)
		centering_mat(i, 3) = -mesh().get_bbox_center()[i];

	// Sample points on each face.
	static SimpleRandomCong_t rng_cong;
	simplerandom_cong_seed(&rng_cong, FLAGS_random_view_seed);

	double x_angle_min = 0;
	double x_angle_max = 2 * M_PI;

	double x_angle = static_cast<double>(simplerandom_cong_next(&rng_cong))
		/ std::numeric_limits<uint32_t>::max();
	x_angle = x_angle * (x_angle_max - x_angle_min) + x_angle_min;
	Eigen::Matrix4d x_axis_random_rotation_mat = Eigen::Matrix4d::Identity();
	x_axis_random_rotation_mat(1, 1) = cos(x_angle);
	x_axis_random_rotation_mat(1, 2) = -sin(x_angle);
	x_axis_random_rotation_mat(2, 1) = sin(x_angle);
	x_axis_random_rotation_mat(2, 2) = cos(x_angle);

	double y_angle_min = 0;
	double y_angle_max = 2 * M_PI;

	double y_angle = static_cast<double>(simplerandom_cong_next(&rng_cong))
		/ std::numeric_limits<uint32_t>::max();
	y_angle = y_angle * (y_angle_max - y_angle_min) + y_angle_min;
	Eigen::Matrix4d y_axis_random_rotation_mat = Eigen::Matrix4d::Identity();
	y_axis_random_rotation_mat(2, 2) = cos(y_angle);
	y_axis_random_rotation_mat(2, 0) = -sin(y_angle);
	y_axis_random_rotation_mat(0, 2) = sin(y_angle);
	y_axis_random_rotation_mat(0, 0) = cos(y_angle);

	// longitude.
	double z_angle_min = 0;
	double z_angle_max = 2 * M_PI;

	double z_angle = static_cast<double>(simplerandom_cong_next(&rng_cong))
		/ std::numeric_limits<uint32_t>::max();
	z_angle = z_angle * (z_angle_max - z_angle_min) + z_angle_min;
	Eigen::Matrix4d z_axis_random_rotation_mat = Eigen::Matrix4d::Identity();
	z_axis_random_rotation_mat(0, 0) = cos(z_angle);
	z_axis_random_rotation_mat(0, 1) = -sin(z_angle);
	z_axis_random_rotation_mat(1, 0) = sin(z_angle);
	z_axis_random_rotation_mat(1, 1) = cos(z_angle);


	// No XY-plane translation.
	//double x_offset = static_cast<double>(simplerandom_cong_next(&rng_cong))
	//	/ std::numeric_limits<uint32_t>::max();
	//x_offset = (x_offset - 0.5) * mesh().get_object_diameter();

	//double y_offset = static_cast<double>(simplerandom_cong_next(&rng_cong))
	//	/ std::numeric_limits<uint32_t>::max();
	//y_offset = (y_offset - 0.5) * mesh().get_object_diameter();

	const double z_offset = -1.5 * mesh().get_object_diameter();

	Eigen::Matrix4d translation_mat = Eigen::Matrix4d::Identity();
	//translation_mat(0, 3) = x_offset;
	//translation_mat(1, 3) = y_offset;
	translation_mat(2, 3) = z_offset;

	Eigen::Matrix4d transformation_mat = translation_mat
		* x_axis_random_rotation_mat * y_axis_random_rotation_mat * z_axis_random_rotation_mat * centering_mat;

	// -Z direction.
	Eigen::Vector3d view_direction = -transformation_mat.row(2).transpose().segment(0, 3);
	Eigen::Vector3d view_point = transformation_mat.col(3).segment(0, 3);

	for (unsigned int i = 0; i < 3; ++i)
	{
		view_point_[i] = view_point[i];
		view_direction_[i] = view_direction[i];
	}

	if (_set_modelview_matrix)
	{
		GLdouble matrix[16];
		for (int col = 0; col < 4; ++col)
			for (int row = 0; row < 4; ++row)
				matrix[4 * col + row] = transformation_mat(row, col);

		set_modelview_matrix(matrix);
	}

	updateGL();
}

void MeshViewerCore::compute_ground_truth_cuboids()
{
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

	unsigned int num_labels = cuboid_structure_.num_labels();

	std::stringstream output_filename_sstr;
	QDir output_dir;
	output_dir.mkpath(FLAGS_training_dir.c_str());

	setDrawMode(CUSTOM_VIEW);
	open_modelview_matrix_file(FLAGS_pose_filename.c_str());
	
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

	ret = load_object_info(mesh_, cuboid_structure_, mesh_filepath.c_str(), LoadGroundTruthData);
	if (!ret) return;

	mesh_.clear_colors();
	open_modelview_matrix_file(FLAGS_pose_filename.c_str());
	updateGL();

	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << FLAGS_training_dir << std::string("/") << mesh_name;
	snapshot(output_filename_sstr.str().c_str());

	if (FLAGS_param_optimize_training_cuboids)
	{
		cuboid_structure_.compute_symmetry_groups();

		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << FLAGS_training_dir << std::string("/") << mesh_name << std::string("_log.txt");

		// Iterate only once.
		MeshCuboidPredictor predictor(num_labels);
		optimize_attributes(cuboid_structure_, NULL, predictor,
			FLAGS_param_opt_single_energy_term_weight, FLAGS_param_opt_symmetry_energy_term_weight,
			1, output_filename_sstr.str(), this,
			FLAGS_param_optimize_with_non_linear_constraints);

		updateGL();
		output_filename_sstr.clear(); output_filename_sstr.str("");
		output_filename_sstr << FLAGS_training_dir << std::string("/") << mesh_name << std::string("_optimized");
		snapshot(output_filename_sstr.str().c_str());
	}

	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << FLAGS_training_dir << std::string("/") << mesh_name << std::string(".arff");
	cuboid_structure_.save_cuboids(output_filename_sstr.str());
}

void MeshViewerCore::train()
{
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

	// To be removed.
	//if (FLAGS_use_symmetric_group_cuboids)
	//{
	//	cuboid_structure_.add_symmetric_group_labels();
	//}

	unsigned int num_labels = cuboid_structure_.num_labels();

	std::vector< std::list<MeshCuboidFeatures *> > feature_list(num_labels);
	std::vector< std::list<MeshCuboidTransformation *> > transformation_list(num_labels);
	//std::vector< std::list<MeshCuboidAttributes *> > attributes_list(num_labels);
	//std::vector< std::list<MeshCuboidManualFeatures *> > manual_single_feature_list(num_labels);
	//std::vector< std::vector< std::list<MeshCuboidManualFeatures *> > > manual_pair_feature_list(num_labels);

	//for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
	//	manual_pair_feature_list[label_index].resize(num_labels);
	

	setDrawMode(CUSTOM_VIEW);
	open_modelview_matrix_file(FLAGS_pose_filename.c_str());

	// For every file in the base path.
	QDir input_dir((FLAGS_data_root_path + FLAGS_mesh_path).c_str());
	assert(input_dir.exists());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
	input_dir.setSorting(QDir::Name);

	std::stringstream output_filename_sstr;
	QDir output_dir;
	output_dir.mkpath(FLAGS_training_dir.c_str());

	std::ofstream mesh_name_list_file(FLAGS_training_dir + std::string("/") + FLAGS_object_list_filename);
	assert(mesh_name_list_file);


	QFileInfoList dir_list = input_dir.entryInfoList();
	for (int i = 0; i < dir_list.size(); i++)
	{
		QFileInfo file_info = dir_list.at(i);
		if (file_info.exists() &&
			(file_info.suffix().compare("obj") == 0
			|| file_info.suffix().compare("off") == 0))
		{
			std::string mesh_filepath = std::string(file_info.filePath().toLocal8Bit());
			std::string mesh_name = std::string(file_info.baseName().toLocal8Bit());

			output_filename_sstr.clear(); output_filename_sstr.str("");
			output_filename_sstr << FLAGS_training_dir << std::string("/") << mesh_name << std::string(".arff");

			bool ret = load_object_info(mesh_, cuboid_structure_, mesh_filepath.c_str(), LoadMesh,
				output_filename_sstr.str().c_str(), false);
			if (!ret) continue;

			mesh_name_list_file << mesh_name << std::endl;

			// NOTE:
			// Compute ground truth cuboids first.
			/*
			bool ret = load_object_info(mesh_, cuboid_structure_, mesh_filepath.c_str(), LoadGroundTruthData);
			if (!ret) continue;

			mesh_name_list_file << mesh_name << std::endl;

			mesh_.clear_colors();
			open_modelview_matrix_file(FLAGS_pose_filename.c_str());
			updateGL();
			output_filename_sstr.clear(); output_filename_sstr.str("");
			output_filename_sstr << FLAGS_training_dir << std::string("/") << mesh_name;
			snapshot(output_filename_sstr.str().c_str());

			if (FLAGS_param_optimize_training_cuboids)
			{
				cuboid_structure_.compute_symmetry_groups();

				output_filename_sstr.clear(); output_filename_sstr.str("");
				output_filename_sstr << FLAGS_training_dir << std::string("/") << mesh_name << std::string("_log.txt");

				// Iterate only once.
				MeshCuboidPredictor predictor(num_labels);
				optimize_attributes(cuboid_structure_, NULL, predictor,
					FLAGS_param_opt_single_energy_term_weight, FLAGS_param_opt_symmetry_energy_term_weight,
					1, output_filename_sstr.str(), this,
					FLAGS_param_optimize_with_non_linear_constraints);

				updateGL();
				output_filename_sstr.clear(); output_filename_sstr.str("");
				output_filename_sstr << FLAGS_training_dir << std::string("/") << mesh_name << std::string("_optimized");
				snapshot(output_filename_sstr.str().c_str());
			}

			output_filename_sstr.clear(); output_filename_sstr.str("");
			output_filename_sstr << FLAGS_training_dir << std::string("/") << mesh_name << std::string(".arff");
			cuboid_structure_.save_cuboids(output_filename_sstr.str());
			*/

			assert(cuboid_structure_.num_labels() == num_labels);
			for (LabelIndex label_index_1 = 0; label_index_1 < num_labels; ++label_index_1)
			{
				MeshCuboidTransformation *transformation = new MeshCuboidTransformation(mesh_name);
				//MeshCuboidAttributes *attributes = new MeshCuboidAttributes(mesh_name);
				MeshCuboidFeatures *features = new MeshCuboidFeatures(mesh_name);
				//MeshSingleCuboidManualFeatures *manual_single_features = new MeshSingleCuboidManualFeatures(mesh_name);

				Label label_1 = cuboid_structure_.labels_[label_index_1];
				MeshCuboid *cuboid_1 = NULL;
				
				// NOTE:
				// The current implementation assumes that there is only one part for each label.
				assert(cuboid_structure_.label_cuboids_[label_index_1].size() <= 1);
				if (!cuboid_structure_.label_cuboids_[label_index_1].empty())
					cuboid_1 = cuboid_structure_.label_cuboids_[label_index_1].front();

				//for (LabelIndex label_index_2 = label_index_1; label_index_2 < num_labels; ++label_index_2)
				//{
				//	MeshPairCuboidManualFeatures *manual_pair_features = new MeshPairCuboidManualFeatures(mesh_name);

				//	Label label_2 = cuboid_structure_.labels_[label_index_2];
				//	MeshCuboid *cuboid_2 = NULL;
				//	
				//	// NOTE:
				//	// The current implementation assumes that there is only one part for each label.
				//	assert(cuboid_structure_.label_cuboids_[label_index_2].size() <= 1);
				//	if (!cuboid_structure_.label_cuboids_[label_index_2].empty())
				//		cuboid_2 = cuboid_structure_.label_cuboids_[label_index_2].front();

				//	if (cuboid_1 && cuboid_2)
				//	{
				//		assert(label_index_1 <= label_index_2);
				//		manual_pair_features->compute_values(cuboid_1, cuboid_2);
				//		manual_pair_feature_list[label_index_1][label_index_2].push_back(manual_pair_features);
				//	}
				//}
				
				if (cuboid_1)
				{
					transformation->compute_transformation(cuboid_1);
					features->compute_features(cuboid_1);
					//manual_single_features->compute_values(cuboid_1);
				}

				transformation_list[label_index_1].push_back(transformation);
				feature_list[label_index_1].push_back(features);
				//manual_single_feature_list[label_index_1].push_back(manual_single_features);
			}
		}
	}


	for (LabelIndex label_index_1 = 0; label_index_1 < num_labels; ++label_index_1)
	{
		//for (LabelIndex label_index_2 = label_index_1; label_index_2 < num_labels; ++label_index_2)
		//{
		//	std::stringstream pair_features_filename_sstr;
		//	pair_features_filename_sstr << FLAGS_pair_feature_filename_prefix
		//		<< label_index_1 << std::string("_")
		//		<< label_index_2 << std::string(".csv");
		//	MeshCuboidManualFeatures::save_keys_and_values(
		//		manual_pair_feature_list[label_index_1][label_index_2],
		//		pair_features_filename_sstr.str().c_str());

		//	std::stringstream pair_stats_filename_sstr;
		//	pair_stats_filename_sstr << FLAGS_pair_stats_filename_prefix
		//		<< label_index_1 << std::string("_")
		//		<< label_index_2 << std::string(".csv");
		//	MeshCuboidManualFeatures::save_stats(
		//		manual_pair_feature_list[label_index_1][label_index_2],
		//		pair_stats_filename_sstr.str().c_str());
		//}

		std::stringstream transformation_filename_sstr;
		transformation_filename_sstr << FLAGS_training_dir << std::string("/") << FLAGS_transformation_filename_prefix
			<< label_index_1 << std::string(".csv");
		MeshCuboidTransformation::save_transformation_collection(transformation_filename_sstr.str().c_str(),
			transformation_list[label_index_1]);

		std::stringstream feature_filename_sstr;
		feature_filename_sstr << FLAGS_training_dir << std::string("/") << FLAGS_feature_filename_prefix
			<< label_index_1 << std::string(".csv");
		MeshCuboidFeatures::save_feature_collection(feature_filename_sstr.str().c_str(),
			feature_list[label_index_1]);

		//MeshCuboidAttributes::save_values(attributes_list[label_index_1],
		//	attributes_filename.c_str());

		//std::stringstream single_features_filename_sstr;
		//single_features_filename_sstr << FLAGS_single_feature_filename_prefix
		//	<< label_index_1 << std::string(".csv");
		//MeshCuboidManualFeatures::save_keys_and_values(manual_single_feature_list[label_index_1],
		//	single_features_filename_sstr.str().c_str());

		//std::stringstream single_stats_filename_sstr;
		//single_stats_filename_sstr << FLAGS_single_stats_filename_prefix
		//	<< label_index_1 << std::string(".csv");
		//MeshCuboidManualFeatures::save_stats(manual_single_feature_list[label_index_1],
		//	single_stats_filename_sstr.str().c_str());
	}


	// Deallocate.
	for (LabelIndex label_index_1 = 0; label_index_1 < num_labels; ++label_index_1)
	{
		//for (LabelIndex label_index_2 = label_index_1; label_index_2 < num_labels; ++label_index_2)
		//{
		//	for (std::list<MeshCuboidManualFeatures *>::iterator it = manual_pair_feature_list[label_index_1][label_index_2].begin();
		//		it != manual_pair_feature_list[label_index_1][label_index_2].end(); ++it)
		//		delete (*it);
		//}

		for (std::list<MeshCuboidTransformation *>::iterator it = transformation_list[label_index_1].begin();
			it != transformation_list[label_index_1].end(); ++it)
			delete (*it);

		//for (std::list<MeshCuboidAttributes *>::iterator it = attributes_list[label_index_1].begin();
		//	it != attributes_list[label_index_1].end(); ++it)
		//	delete (*it);

		for (std::list<MeshCuboidFeatures *>::iterator it = feature_list[label_index_1].begin();
			it != feature_list[label_index_1].end(); ++it)
			delete (*it);

		//for (std::list<MeshCuboidManualFeatures *>::iterator it = manual_single_feature_list[label_index_1].begin();
		//	it != manual_single_feature_list[label_index_1].end(); ++it)
		//	delete (*it);
	}

	mesh_name_list_file.close();

	std::cout << std::endl;
	std::cout << " -- Batch Completed. -- " << std::endl;
}

/*
void MeshViewerCore::train_file_files()
{
	unsigned int num_cuboids = 0;

	MeshCuboidTrainer trainer;
	trainer.load_object_list(FLAGS_training_dir + std::string("/") + FLAGS_object_list_filename);
	trainer.load_features(FLAGS_training_dir + std::string("/") + FLAGS_feature_filename_prefix);
	trainer.load_transformations(FLAGS_training_dir + std::string("/") + FLAGS_transformation_filename_prefix);

	//
	std::vector< std::vector<MeshCuboidJointNormalRelations *> > joint_normal_relations;
	trainer.get_joint_normal_relations(joint_normal_relations);

	num_cuboids = joint_normal_relations.size();
	for (unsigned int cuboid_index_1 = 0; cuboid_index_1 < num_cuboids; ++cuboid_index_1)
	{
		assert(joint_normal_relations[cuboid_index_1].size() == num_cuboids);
		for (unsigned int cuboid_index_2 = 0; cuboid_index_2 < num_cuboids; ++cuboid_index_2)
		{
			if (cuboid_index_1 == cuboid_index_2) continue;

			const MeshCuboidJointNormalRelations *relation_12 = joint_normal_relations[cuboid_index_1][cuboid_index_2];
			if (!relation_12) continue;
			
			std::stringstream relation_filename_sstr;
			relation_filename_sstr << FLAGS_joint_normal_relation_filename_prefix << cuboid_index_1
				<< "_" << cuboid_index_2 << ".csv";

			std::cout << "Saving '" << relation_filename_sstr.str() << "'..." << std::endl;
			relation_12->save_joint_normal_csv(relation_filename_sstr.str().c_str());
		}
	}
	//

	////
	//std::vector< std::vector<MeshCuboidCondNormalRelations *> > cond_normal_relations;
	//trainer.get_cond_normal_relations(cond_normal_relations);

	//num_cuboids = cond_normal_relations.size();
	//for (unsigned int cuboid_index_1 = 0; cuboid_index_1 < num_cuboids; ++cuboid_index_1)
	//{
	//	assert(cond_normal_relations[cuboid_index_1].size() == num_cuboids);
	//	for (unsigned int cuboid_index_2 = 0; cuboid_index_2 < num_cuboids; ++cuboid_index_2)
	//	{
	//		if (cuboid_index_1 == cuboid_index_2) continue;

	//		const MeshCuboidCondNormalRelations *relation_12 = cond_normal_relations[cuboid_index_1][cuboid_index_2];
	//		if (!relation_12) continue;

	//		std::stringstream relation_filename_sstr;
	//		relation_filename_sstr << FLAGS_cond_normal_relation_filename_prefix << cuboid_index_1
	//			<< "_" << cuboid_index_2 << ".csv";

	//		std::cout << "Saving '" << relation_filename_sstr.str() << "'..." << std::endl;
	//		relation_12->save_cond_normal_csv(relation_filename_sstr.str().c_str());
	//	}
	//}
	////
}
*/

void MeshViewerCore::batch_predict()
{
	// For every file in the base path.
	QDir input_dir((FLAGS_data_root_path + FLAGS_mesh_path).c_str());
	assert(input_dir.exists());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
	input_dir.setSorting(QDir::Name);

	QFileInfoList dir_list = input_dir.entryInfoList();
	for (int file_index = 0; file_index < dir_list.size(); file_index++)
	{
		QFileInfo file_info = dir_list.at(file_index);
		if (file_info.exists() &&
			(file_info.suffix().compare("obj") == 0
			|| file_info.suffix().compare("off") == 0))
		{
			FLAGS_mesh_filename = std::string(file_info.fileName().toLocal8Bit());
			std::cout << "mesh_filename = " << FLAGS_mesh_filename << std::endl;
			predict();
		}
	}

	std::cout << " -- Batch Completed. -- " << std::endl;
}

void MeshViewerCore::predict()
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

	// To be removed.
	//if (FLAGS_use_symmetric_group_cuboids)
	//{
	//	cuboid_structure_.add_symmetric_group_labels();
	//}

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
	std::stringstream snapshot_filename_sstr;
	std::stringstream log_filename_sstr;

	QDir output_dir;
	std::string mesh_output_path = FLAGS_output_dir + std::string("/") + mesh_name;
	std::string mesh_intermediate_path = FLAGS_output_dir + std::string("/") + mesh_name + std::string("/temp");
	output_dir.mkpath(FLAGS_output_dir.c_str());
	output_dir.mkpath(mesh_output_path.c_str());
	output_dir.mkpath(mesh_intermediate_path.c_str());


	// Initialize basic information.
	unsigned int num_labels = cuboid_structure_.num_labels();
	cuboid_structure_.clear_cuboids();
	cuboid_structure_.clear_sample_points();

	std::list<std::string> ignored_object_list;
	ignored_object_list.push_back(mesh_name);

	std::vector< std::vector<MeshCuboidJointNormalRelations *> > joint_normal_relations;
	//MeshCuboidTrainer::load_joint_normal_relations(num_labels, "joint_normal_", joint_normal_relations);
	trainer.get_joint_normal_relations(joint_normal_relations, &ignored_object_list);
	MeshCuboidJointNormalRelationPredictor joint_normal_predictor(joint_normal_relations);

	//std::vector< std::vector<MeshCuboidCondNormalRelations *> > cond_normal_relations;
	//MeshCuboidTrainer::load_cond_normal_relations(num_labels, "conditional_normal_", cond_normal_relations);
	//trainer.get_cond_normal_relations(cond_normal_relations, &ignored_object_list);
	//MeshCuboidCondNormalRelationPredictor cond_normal_predictor(cond_normal_relations);


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
	save_modelview_matrix_file((mesh_output_path + std::string("/occlusion_pose.txt")).c_str());


	//
	ret = load_object_info(mesh_, cuboid_structure_, mesh_filepath.c_str(), LoadDenseTestData);
	if (!ret) return;
	set_modelview_matrix(occlusion_modelview_matrix, false);

	// View plane mask.
	if (FLAGS_use_view_plane_mask)
	{
		std::string view_plane_mask_filename = mesh_output_path + std::string("/view_mask.txt");
		QFileInfo view_plane_mask_file_info(view_plane_mask_filename.c_str());

		if (!view_plane_mask_file_info.exists())
		{
			compute_view_plane_mask_range(occlusion_modelview_matrix);

			// Save range.
			std::ofstream view_plane_mask_file(view_plane_mask_filename.c_str());
			assert(view_plane_mask_file.is_open());
			view_plane_mask_file << FLAGS_param_view_plane_mask_min_x << std::endl;
			view_plane_mask_file << FLAGS_param_view_plane_mask_min_y << std::endl;
			view_plane_mask_file << FLAGS_param_view_plane_mask_max_x << std::endl;
			view_plane_mask_file << FLAGS_param_view_plane_mask_max_y << std::endl;
			view_plane_mask_file.close();
		}

		std::ifstream view_plane_mask_file(view_plane_mask_filename.c_str());
		assert(view_plane_mask_file.is_open());
		std::string buffer;
		std::getline(view_plane_mask_file, buffer); FLAGS_param_view_plane_mask_min_x = std::atof(buffer.c_str());
		std::getline(view_plane_mask_file, buffer); FLAGS_param_view_plane_mask_min_y = std::atof(buffer.c_str());
		std::getline(view_plane_mask_file, buffer); FLAGS_param_view_plane_mask_max_x = std::atof(buffer.c_str());
		std::getline(view_plane_mask_file, buffer); FLAGS_param_view_plane_mask_max_y = std::atof(buffer.c_str());
		view_plane_mask_file.close();
	}

	remove_occluded_points();
	
	updateGL();
	snapshot_filename_sstr.clear(); snapshot_filename_sstr.str("");
	snapshot_filename_sstr << mesh_output_path << filename_prefix << std::string("view");
	snapshot(snapshot_filename_sstr.str().c_str());

	set_modelview_matrix(snapshot_modelview_matrix);
	updateGL();
	snapshot_filename_sstr.clear(); snapshot_filename_sstr.str("");
	snapshot_filename_sstr << mesh_output_path << filename_prefix << std::string("input");
	snapshot(snapshot_filename_sstr.str().c_str());

	cuboid_structure_.save_sample_points_to_ply(snapshot_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points((snapshot_filename_sstr.str() + std::string(".pts")).c_str());
	//


	ret = load_object_info(mesh_, cuboid_structure_, mesh_filepath.c_str(), LoadTestData);
	if (!ret) return;

	std::cout << " - Remove occluded points." << std::endl;
	set_modelview_matrix(occlusion_modelview_matrix, false);
	remove_occluded_points();
	set_modelview_matrix(snapshot_modelview_matrix);


	std::cout << " - Cluster points and construct initial cuboids." << std::endl;
	draw_cuboid_axes_ = false;

	cuboid_structure_.compute_label_cuboids();

	// Split cuboids if sample points are far away each other.
	cuboid_structure_.split_label_cuboids();

	// Remove cuboids in symmetric labels.
	cuboid_structure_.remove_symmetric_cuboids();

	if (cuboid_structure_.get_all_cuboids().empty())
		return;

	update_cuboid_surface_points(cuboid_structure_, occlusion_modelview_matrix);


	// Sub-routine.
	bool first_iteration = true;
	unsigned int num_final_cuboid_structure_candidates = 0;

	std::list< std::pair<std::string, MeshCuboidStructure> > cuboid_structure_candidates;
	cuboid_structure_candidates.push_back(std::make_pair(std::string("0"), cuboid_structure_));
	std::set<LabelIndex> ignored_label_indices;

	while (!cuboid_structure_candidates.empty())
	{
		// FIXME:
		// The cuboid structure should not deep copy all sample points.
		// Use smart pointers for sample points.
		std::string cuboid_structure_name = cuboid_structure_candidates.front().first;
		cuboid_structure_ = cuboid_structure_candidates.front().second;
		cuboid_structure_candidates.pop_front();

		unsigned int snapshot_index = 0;
		log_filename_sstr.clear(); log_filename_sstr.str("");
		log_filename_sstr << mesh_intermediate_path << filename_prefix
			<< std::string("c_") << cuboid_structure_name << std::string("_log.txt");
		std::ofstream log_file(log_filename_sstr.str());
		log_file.clear(); log_file.close();


		updateGL();
		snapshot_filename_sstr.clear(); snapshot_filename_sstr.str("");
		snapshot_filename_sstr << mesh_intermediate_path << filename_prefix
			<< std::string("c_") << cuboid_structure_name << std::string("_")
			<< std::string("s_") << snapshot_index;
		snapshot(snapshot_filename_sstr.str().c_str());
		++snapshot_index;
		draw_cuboid_axes_ = true;
		

		std::cout << "\n1. Recognize labels and axes configurations." << std::endl;
		// NOTE:
		// Use symmetric label information only at the first time of the iteration.
		recognize_labels_and_axes_configurations(cuboid_structure_,
			joint_normal_predictor, log_filename_sstr.str(), first_iteration,
			true);
		first_iteration = false;

		//
		cuboid_structure_.compute_symmetry_groups();
		//

		updateGL();
		snapshot_filename_sstr.clear(); snapshot_filename_sstr.str("");
		snapshot_filename_sstr << mesh_intermediate_path << filename_prefix
			<< std::string("c_") << cuboid_structure_name << std::string("_")
			<< std::string("s_") << snapshot_index;
		snapshot(snapshot_filename_sstr.str().c_str());
		++snapshot_index;


		std::cout << "\n2. Segment sample points." << std::endl;
		segment_sample_points(cuboid_structure_);

		updateGL();
		snapshot_filename_sstr.clear(); snapshot_filename_sstr.str("");
		snapshot_filename_sstr << mesh_intermediate_path << filename_prefix
			<< std::string("c_") << cuboid_structure_name << std::string("_")
			<< std::string("s_") << snapshot_index;
		snapshot(snapshot_filename_sstr.str().c_str());
		++snapshot_index;


		std::cout << "\n3. Optimize cuboid attributes." << std::endl;

		optimize_attributes(cuboid_structure_, occlusion_modelview_matrix, joint_normal_predictor,
			FLAGS_param_opt_single_energy_term_weight, FLAGS_param_opt_symmetry_energy_term_weight,
			FLAGS_param_opt_max_iterations, log_filename_sstr.str(), this, false);

		if (FLAGS_param_optimize_with_non_linear_constraints)
		{
			cuboid_structure_.compute_symmetry_groups();

			optimize_attributes(cuboid_structure_, occlusion_modelview_matrix, joint_normal_predictor,
				FLAGS_param_opt_single_energy_term_weight, FLAGS_param_opt_symmetry_energy_term_weight,
				FLAGS_param_opt_max_iterations, log_filename_sstr.str(), this, true);
		}

		updateGL();
		snapshot_filename_sstr.clear(); snapshot_filename_sstr.str("");
		snapshot_filename_sstr << FLAGS_output_dir + std::string("/Temp") << filename_prefix
			<< std::string("c_") << cuboid_structure_name << std::string("_")
			<< std::string("s_") << snapshot_index;
		snapshot(snapshot_filename_sstr.str().c_str());
		++snapshot_index;


		std::cout << "\n4. Add missing cuboids." << std::endl;
		assert(cuboid_structure_.num_labels() == num_labels);
		std::list<LabelIndex> given_label_indices;
		for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
			if (!cuboid_structure_.label_cuboids_[label_index].empty())
				given_label_indices.push_back(label_index);

		std::list< std::list<LabelIndex> > missing_label_index_groups;
		trainer.get_missing_label_index_groups(given_label_indices, missing_label_index_groups,
			&ignored_label_indices);
		

		bool is_cuboid_added = (!missing_label_index_groups.empty());

		if (!missing_label_index_groups.empty())
		{
			unsigned int missing_label_index_group_index = 0;

			for (std::list< std::list<LabelIndex> >::iterator it = missing_label_index_groups.begin();
				it != missing_label_index_groups.end(); ++it)
			{
				std::list<LabelIndex> &missing_label_indices = (*it);
				MeshCuboidStructure new_cuboid_structure = cuboid_structure_;

				// FIXME:
				// Any missing cuboid may not be added.
				// Then, you should escapt the loop.
				ret = add_missing_cuboids(new_cuboid_structure, occlusion_modelview_matrix,
					missing_label_indices, joint_normal_predictor, ignored_label_indices);
				//ret = add_missing_cuboids(new_cuboid_structure, occlusion_modelview_matrix,
				//	missing_label_indices, joint_normal_relations, ignored_label_indices);

				if (!ret)
				{
					is_cuboid_added = false;
				}
				else
				{
					std::stringstream new_cuboid_structure_name;
					new_cuboid_structure_name << cuboid_structure_name << missing_label_index_group_index;
					cuboid_structure_candidates.push_front(
						std::make_pair(new_cuboid_structure_name.str(), new_cuboid_structure));
					++missing_label_index_group_index;
				}
			}
		}

		// If there was a case when no cuboid is added, reconstruct using the current cuboid structure.
		if (!is_cuboid_added)
		{
			// Escape loop.
			snapshot_filename_sstr.clear(); snapshot_filename_sstr.str("");
			snapshot_filename_sstr << mesh_output_path << filename_prefix << num_final_cuboid_structure_candidates;

			draw_point_correspondences_ = false;
			reconstruct(
				mesh_filepath.c_str(),
				snapshot_modelview_matrix,
				occlusion_modelview_matrix,
				snapshot_filename_sstr.str().c_str());
			draw_point_correspondences_ = true;

			ignored_label_indices.clear();
			++num_final_cuboid_structure_candidates;
		}
	}


	//annDeallocPts(occlusion_test_ann_points);
	//delete occlusion_test_points_kd_tree;

	for (LabelIndex label_index_1 = 0; label_index_1 < joint_normal_relations.size(); ++label_index_1)
		for (LabelIndex label_index_2 = 0; label_index_2 < joint_normal_relations[label_index_1].size(); ++label_index_2)
			delete joint_normal_relations[label_index_1][label_index_2];

	//for (LabelIndex label_index_1 = 0; label_index_1 < cond_normal_relations.size(); ++label_index_1)
	//	for (LabelIndex label_index_2 = 0; label_index_2 < cond_normal_relations[label_index_1].size(); ++label_index_2)
	//		delete cond_normal_relations[label_index_1][label_index_2];
}

void MeshViewerCore::run_symmetry_detection()
{
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

	std::stringstream output_filename_sstr;

	std::string mesh_output_path = FLAGS_symmetry_detection_dir + std::string("/") + mesh_name;
	QDir output_dir;
	output_dir.mkpath(FLAGS_symmetry_detection_dir.c_str());
	output_dir.mkpath(mesh_output_path.c_str());

	setDrawMode(CUSTOM_VIEW);
	open_modelview_matrix_file(FLAGS_pose_filename.c_str());

	bool ret = load_object_info(mesh_, cuboid_structure_, mesh_filepath.c_str(), LoadSamplePoints);
	if (!ret) return;

	mesh_.clear_colors();
	open_modelview_matrix_file(FLAGS_pose_filename.c_str());
	updateGL();

	const unsigned int num_input_points = cuboid_structure_.sample_points_.size();
	assert(num_input_points > 0);
	Eigen::MatrixXd sample_points = Eigen::MatrixXd(3, num_input_points);

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_input_points;
		++sample_point_index)
	{
		const MeshSamplePoint *sample_point = cuboid_structure_.sample_points_[sample_point_index];
		assert(sample_point);
		for (int i = 0; i < 3; ++i)
			sample_points.col(sample_point_index)[i] = sample_point->point_[i];
	}

	std::list<SymmetryDetection::ReflectionPlane> reflection_planes;
	SymmetryDetection::detect_reflectional_symmetry(sample_points,
		FLAGS_param_sparse_neighbor_distance, 0.7, reflection_planes);

	//
	for (std::vector< MeshCuboidReflectionSymmetryGroup* >::iterator it
		= cuboid_structure_.reflection_symmetry_groups_.begin();
		it != cuboid_structure_.reflection_symmetry_groups_.end(); ++it)
		delete (*it);
	cuboid_structure_.reflection_symmetry_groups_.clear();
	//

	for (std::list<SymmetryDetection::ReflectionPlane>::iterator it = reflection_planes.begin();
		it != reflection_planes.end(); ++it)
	{
		MyMesh::Normal normal;
		MyMesh::Point point;
		for (int i = 0; i < 3; ++i)
		{
			normal[i] = (*it).normal_[i];
			point[i] = (*it).point_[i];
		}

		MeshCuboidReflectionSymmetryGroup *symmetry_group = new MeshCuboidReflectionSymmetryGroup(
			normal, dot(normal, point));
		cuboid_structure_.reflection_symmetry_groups_.push_back(symmetry_group);
	}

	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << mesh_output_path << std::string("/") << mesh_name;
	snapshot(output_filename_sstr.str().c_str());
}


void MeshViewerCore::run_symmetry_detection_msh2pln()
{
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

	// View plane mask.
	if (FLAGS_use_view_plane_mask)
	{
		std::string mesh_output_path = FLAGS_output_dir + std::string("/") + mesh_name;
		std::string view_plane_mask_filename = mesh_output_path + std::string("/view_mask.txt");
		std::ifstream view_plane_mask_file(view_plane_mask_filename.c_str());

		if (!view_plane_mask_file.is_open())
		{
			std::cerr << "Error: The view plane mask file is not opened (" << view_plane_mask_filename << ")." << std::endl;
			return;
		}

		std::string buffer;
		std::getline(view_plane_mask_file, buffer); FLAGS_param_view_plane_mask_min_x = std::atof(buffer.c_str());
		std::getline(view_plane_mask_file, buffer); FLAGS_param_view_plane_mask_min_y = std::atof(buffer.c_str());
		std::getline(view_plane_mask_file, buffer); FLAGS_param_view_plane_mask_max_x = std::atof(buffer.c_str());
		std::getline(view_plane_mask_file, buffer); FLAGS_param_view_plane_mask_max_y = std::atof(buffer.c_str());
		view_plane_mask_file.close();
	}


	std::string filename_prefix = std::string("/") + mesh_name + std::string("_");
	std::stringstream output_filename_sstr;

	QDir output_dir;
	std::string mesh_output_path = FLAGS_symmetry_detection_dir + std::string("/") + mesh_name;
	output_dir.mkpath(FLAGS_output_dir.c_str());
	output_dir.mkpath(mesh_output_path.c_str());


	// Initialize basic information.
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
		bool ret = open_modelview_matrix_file(FLAGS_occlusion_pose_filename.c_str());
		assert(ret);
	}
	memcpy(occlusion_modelview_matrix, modelview_matrix(), 16 * sizeof(double));
	//save_modelview_matrix_file((mesh_output_path + std::string("/occlusion_pose.txt")).c_str());


	//
	bool ret = load_object_info(mesh_, cuboid_structure_, mesh_filepath.c_str(), LoadDenseSamplePoints);
	if (!ret) return;
	set_modelview_matrix(occlusion_modelview_matrix, false);
	remove_occluded_points();
	updateGL();

	set_modelview_matrix(snapshot_modelview_matrix);
	updateGL();
	//


	//
	unsigned int num_sample_points = cuboid_structure_.num_sample_points();
	std::vector<MyMesh::Point> sample_points(num_sample_points);
	Eigen::MatrixXd sample_points_mat(3, num_sample_points);

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points; ++sample_point_index)
	{
		MeshSamplePoint *sample_point = cuboid_structure_.sample_points_[sample_point_index];
		assert(sample_point);
		sample_points[sample_point_index] = sample_point->point_;
		for (unsigned int i = 0; i < 3; ++i)
			sample_points_mat.col(sample_point_index)[i] = sample_point->point_[i];
	}

	ANNpointArray sample_ann_points;
	ANNkd_tree *sample_ann_kd_tree = ICP::create_kd_tree(sample_points_mat, sample_ann_points);
	//


	const std::string msh2pln_output_filename = mesh_output_path + std::string("/") + mesh_name + std::string(".txt");

	std::ifstream file(msh2pln_output_filename);
	if (!file)
	{
		std::cerr << "Can't open file: \"" << msh2pln_output_filename << "\"" << std::endl;
		return;
	}

	std::string buffer;
	std::stringstream strstr;
	std::string token;
	
	const Real squared_neighbor_distance = FLAGS_param_sparse_neighbor_distance *
		FLAGS_param_sparse_neighbor_distance * mesh_.get_object_diameter();

	unsigned int max_num_symmetric_point_pairs = 0;
	assert(cuboid_structure_.reflection_symmetry_groups_.empty());


	while (!file.eof())
	{
		std::getline(file, buffer);
		if (buffer == "") continue;

		strstr.str(std::string());
		strstr.clear();
		strstr.str(buffer);

		do{ std::getline(strstr, token, ' '); } while (token == "" && !strstr.eof());
		Real center_x = std::stof(token);

		do{ std::getline(strstr, token, ' '); } while (token == "" && !strstr.eof());
		Real center_y = std::stof(token);

		do{ std::getline(strstr, token, ' '); } while (token == "" && !strstr.eof());
		Real center_z = std::stof(token);


		do{ std::getline(strstr, token, ' '); } while (token == "" && !strstr.eof());
		Real normal_x = std::stof(token);

		do{ std::getline(strstr, token, ' '); } while (token == "" && !strstr.eof());
		Real normal_y = std::stof(token);

		do{ std::getline(strstr, token, ' '); } while (token == "" && !strstr.eof());
		Real normal_z = std::stof(token);

		//printf("center = (%lf, %lf, %lf)\n", center_x, center_y, center_z);
		//printf("normal = (%lf, %lf, %lf)\n", normal_x, normal_y, normal_z);

		MyMesh::Point center(center_x, center_y, center_z);
		MyMesh::Normal normal(normal_x, normal_y, normal_z);
		MeshCuboidReflectionSymmetryGroup *new_symmetry_group = new MeshCuboidReflectionSymmetryGroup(
			normal, dot(normal, center));
		

		std::list<MeshCuboidSymmetryGroup::WeightedPointPair> sample_point_pairs;
		new_symmetry_group->get_symmetric_sample_point_pairs(
			sample_points, sample_ann_points, sample_ann_kd_tree,
			squared_neighbor_distance, sample_point_pairs);

		if (sample_point_pairs.size() > max_num_symmetric_point_pairs)
		{
			max_num_symmetric_point_pairs = sample_point_pairs.size();

			if (!cuboid_structure_.reflection_symmetry_groups_.empty())
			{
				assert(cuboid_structure_.reflection_symmetry_groups_.size() == 1);
				delete cuboid_structure_.reflection_symmetry_groups_.front();
			}

			cuboid_structure_.reflection_symmetry_groups_.clear();
			cuboid_structure_.reflection_symmetry_groups_.push_back(new_symmetry_group);
		}
	}

	file.close();
	annDeallocPts(sample_ann_points);
	delete sample_ann_kd_tree;

	
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << mesh_output_path << filename_prefix << std::string("symm_detection");
	updateGL();
	snapshot(output_filename_sstr.str().c_str());

	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << mesh_output_path << filename_prefix << std::string("symmetry_info.txt");
	cuboid_structure_.save_symmetry_groups(output_filename_sstr.str());


	if (!cuboid_structure_.reflection_symmetry_groups_.empty())
	{
		MeshCuboidReflectionSymmetryGroup *symmetry_group = cuboid_structure_.reflection_symmetry_groups_.front();
		assert(symmetry_group);

		for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points; ++sample_point_index)
		{
			MeshSamplePoint *sample_point = cuboid_structure_.sample_points_[sample_point_index];
			assert(sample_point);

			MyMesh::Point symmetric_point = symmetry_group->get_symmetric_point(sample_point->point_);
			MyMesh::Point symmetric_normal = symmetry_group->get_symmetric_normal(sample_point->normal_);
			cuboid_structure_.add_sample_point(symmetric_point, symmetric_normal);
		}
	}

	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << mesh_output_path << filename_prefix << std::string("symm_detection_recon");
	updateGL();
	snapshot(output_filename_sstr.str().c_str());

	
	//
	MyMesh ground_truth_mesh;
	MeshCuboidStructure ground_truth_cuboid_structure(&ground_truth_mesh);

	ret = load_object_info(ground_truth_mesh, ground_truth_cuboid_structure,
		mesh_filepath.c_str(), LoadDenseSamplePoints);
	assert(ret);
	MeshCuboidEvaluator evaluator(&ground_truth_cuboid_structure);

	setDrawMode(COLORED_POINT_SAMPLES);
	//

	evaluator.evaluate_point_to_point_distances(&cuboid_structure_, output_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points_to_ply(output_filename_sstr.str().c_str());
	cuboid_structure_.save_sample_points((output_filename_sstr.str() + std::string(".pts")).c_str());
	MeshCuboidStructure symmetry_reconstruction(cuboid_structure_);

	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << mesh_output_path << filename_prefix << std::string("symm_detection_accuracy");
	snapshot(output_filename_sstr.str().c_str());

	cuboid_structure_ = ground_truth_cuboid_structure;
	updateGL();
	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << mesh_output_path << filename_prefix << std::string("symm_detection_completeness");
	snapshot(output_filename_sstr.str().c_str());


	setDrawMode(CUSTOM_VIEW);
}

void MeshViewerCore::run_baseline_stats()
{
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
	std::string mesh_output_path = FLAGS_baseline_dir + std::string("/") + mesh_name;
	output_dir.mkpath(FLAGS_output_dir.c_str());
	output_dir.mkpath(mesh_output_path.c_str());


	// Initialize basic information.
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
		bool ret = open_modelview_matrix_file(FLAGS_occlusion_pose_filename.c_str());
		assert(ret);
	}
	memcpy(occlusion_modelview_matrix, modelview_matrix(), 16 * sizeof(double));
	//save_modelview_matrix_file((mesh_output_path + std::string("/occlusion_pose.txt")).c_str());


	//
	bool ret = load_object_info(mesh_, cuboid_structure_, mesh_filepath.c_str(), LoadDenseSamplePoints);
	if (!ret) return;
	set_modelview_matrix(occlusion_modelview_matrix, false);
	remove_occluded_points();
	updateGL();

	set_modelview_matrix(snapshot_modelview_matrix);
	updateGL();
	//


	output_filename_sstr.clear(); output_filename_sstr.str("");
	output_filename_sstr << mesh_output_path << filename_prefix << std::string("baseline");
	updateGL();
	snapshot(output_filename_sstr.str().c_str());

	//
	MyMesh ground_truth_mesh;
	MeshCuboidStructure ground_truth_cuboid_structure(&ground_truth_mesh);

	ret = load_object_info(ground_truth_mesh, ground_truth_cuboid_structure,
		mesh_filepath.c_str(), LoadDenseSamplePoints);
	assert(ret);
	MeshCuboidEvaluator evaluator(&ground_truth_cuboid_structure);
	evaluator.evaluate_point_to_point_distances(&cuboid_structure_, output_filename_sstr.str().c_str());


	setDrawMode(CUSTOM_VIEW);
}

void MeshViewerCore::batch_render_point_clusters()
{
	bool ret = true;
	ret = ret & cuboid_structure_.load_labels((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_info_filename).c_str());
	ret = ret & cuboid_structure_.load_label_symmetries((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_symmetry_info_filename).c_str());
	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open label information files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}


	// To be removed.
	//if (FLAGS_use_symmetric_group_cuboids)
	//{
	//	cuboid_structure_.add_symmetric_group_labels();
	//}


	unsigned int num_labels = cuboid_structure_.num_labels();


	setDrawMode(CUSTOM_VIEW);
	open_modelview_matrix_file(FLAGS_pose_filename.c_str());

	// For every file in the base path.
	QDir input_dir((FLAGS_data_root_path + FLAGS_mesh_path).c_str());
	assert(input_dir.exists());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
	input_dir.setSorting(QDir::Name);

	QDir output_dir;
	output_dir.mkpath(FLAGS_output_dir.c_str());


	QFileInfoList dir_list = input_dir.entryInfoList();
	for (int i = 0; i < dir_list.size(); i++)
	{
		QFileInfo file_info = dir_list.at(i);
		if (file_info.exists() &&
			(file_info.suffix().compare("obj") == 0
			|| file_info.suffix().compare("off") == 0))
		{
			std::string mesh_name = std::string(file_info.baseName().toLocal8Bit());
			std::string mesh_filename = std::string(file_info.filePath().toLocal8Bit());

			std::string sample_filename = FLAGS_data_root_path + FLAGS_sample_path + std::string("/") + mesh_name + std::string(".pts");
			std::string sample_label_filename = FLAGS_data_root_path + FLAGS_sample_label_path + std::string("/") + mesh_name + std::string(".arff");
			std::string mesh_label_filename = FLAGS_data_root_path + FLAGS_mesh_label_path + std::string("/") + mesh_name + std::string(".seg");
			std::string snapshot_filename = FLAGS_output_dir + std::string("/") + mesh_name + std::string("_in");

			QFileInfo mesh_file(mesh_filename.c_str());
			QFileInfo sample_file(sample_filename.c_str());
			QFileInfo sample_label_file(sample_label_filename.c_str());
			QFileInfo mesh_label_file(mesh_label_filename.c_str());

			if (!mesh_file.exists()
				|| !sample_file.exists()
				|| !sample_label_file.exists()
				|| !mesh_label_file.exists())
				continue;


			cuboid_structure_.clear_cuboids();
			cuboid_structure_.clear_sample_points();

			open_mesh(mesh_filename.c_str());
			open_sample_point_file(sample_filename.c_str());
			open_sample_point_label_file(sample_label_filename.c_str());

			cuboid_structure_.compute_label_cuboids();

			mesh_.clear_colors();
			draw_cuboid_axes_ = false;
			updateGL();
			snapshot(snapshot_filename.c_str());
		}
	}

	draw_cuboid_axes_ = true;
}

void MeshViewerCore::batch_render_cuboids()
{
	bool ret = true;
	ret = ret & cuboid_structure_.load_labels((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_info_filename).c_str());
	ret = ret & cuboid_structure_.load_label_symmetries((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_symmetry_info_filename).c_str());
	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open label information files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}


	// To be removed.
	//if (FLAGS_use_symmetric_group_cuboids)
	//{
	//	cuboid_structure_.add_symmetric_group_labels();
	//}


	unsigned int num_labels = cuboid_structure_.num_labels();


	setDrawMode(CUSTOM_VIEW);
	open_modelview_matrix_file(FLAGS_pose_filename.c_str());

	// For every file in the base path.
	QDir input_dir(FLAGS_output_dir.c_str());
	assert(input_dir.exists());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
	input_dir.setSorting(QDir::Name);

	QDir output_dir;
	output_dir.mkpath((FLAGS_output_dir + std::string("/cuboid_snapshots/")).c_str());


	QFileInfoList dir_list = input_dir.entryInfoList();
	for (int i = 0; i < dir_list.size(); i++)
	{
		QFileInfo file_info = dir_list.at(i);
		if (file_info.exists() &&
			file_info.suffix().compare("arff") == 0)
		{
			std::string cuboid_name = std::string(file_info.baseName().toLocal8Bit());
			std::string cuboid_filename = std::string(file_info.filePath().toLocal8Bit());
			std::string snapshot_filename = FLAGS_output_dir + std::string("/cuboid_snapshots/") + cuboid_name;

			QFileInfo cuboid_file(cuboid_filename.c_str());
			if (!cuboid_file.exists())
				continue;

			cuboid_structure_.clear_cuboids();
			cuboid_structure_.clear_sample_points();
			cuboid_structure_.load_cuboids(cuboid_filename.c_str());

			draw_cuboid_axes_ = true;
			updateGL();
			snapshot(snapshot_filename.c_str());
		}
	}
}

void MeshViewerCore::batch_render_points()
{
	bool ret = true;
	ret = ret & cuboid_structure_.load_labels((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_info_filename).c_str());
	ret = ret & cuboid_structure_.load_label_symmetries((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_symmetry_info_filename).c_str());
	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open label information files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}


	// To be removed.
	//if (FLAGS_use_symmetric_group_cuboids)
	//{
	//	cuboid_structure_.add_symmetric_group_labels();
	//}


	unsigned int num_labels = cuboid_structure_.num_labels();


	setDrawMode(CUSTOM_VIEW);
	open_modelview_matrix_file(FLAGS_pose_filename.c_str());

	// For every file in the base path.
	QDir input_dir((FLAGS_data_root_path + FLAGS_mesh_path).c_str());
	assert(input_dir.exists());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
	input_dir.setSorting(QDir::Name);

	QDir output_dir;
	output_dir.mkpath((FLAGS_output_dir + std::string("/reconstruction_snapshots/")).c_str());


	const std::string sample_filename_postfix("_0_reconstructed.pts");

	QFileInfoList dir_list = input_dir.entryInfoList();
	for (int i = 0; i < dir_list.size(); i++)
	{
		QFileInfo file_info = dir_list.at(i);
		if (file_info.exists() &&
			(file_info.suffix().compare("obj") == 0
			|| file_info.suffix().compare("off") == 0))
		{
			std::string mesh_name = std::string(file_info.baseName().toLocal8Bit());
			std::string mesh_filename = std::string(file_info.filePath().toLocal8Bit());

			//std::string sample_filename = FLAGS_data_root_path + FLAGS_sample_path + std::string("/") + mesh_name + std::string(".pts");
			std::string sample_filename = FLAGS_output_dir + std::string("/") + mesh_name + sample_filename_postfix;
			std::string snapshot_filename = FLAGS_output_dir + std::string("/reconstruction_snapshots/") + mesh_name;

			QFileInfo mesh_file(mesh_filename.c_str());
			QFileInfo sample_file(sample_filename.c_str());

			if (!mesh_file.exists()
				|| !sample_file.exists())
				continue;


			cuboid_structure_.clear_cuboids();
			cuboid_structure_.clear_sample_points();

			open_mesh(mesh_filename.c_str());
			open_sample_point_file(sample_filename.c_str());

			mesh_.clear_colors();
			draw_cuboid_axes_ = false;
			updateGL();
			snapshot(snapshot_filename.c_str());
		}
	}

	draw_cuboid_axes_ = true;
}


void MeshViewerCore::test_figure_1()
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


	// Initialize basic information.
	unsigned int num_labels = cuboid_structure_.num_labels();
	cuboid_structure_.clear_cuboids();
	cuboid_structure_.clear_sample_points();

	std::list<std::string> ignored_object_list;
	std::vector< std::vector<MeshCuboidJointNormalRelations *> > joint_normal_relations;
	//MeshCuboidTrainer::load_joint_normal_relations(num_labels, "joint_normal_", joint_normal_relations);
	trainer.get_joint_normal_relations(joint_normal_relations, &ignored_object_list);
	MeshCuboidJointNormalRelationPredictor joint_normal_predictor(joint_normal_relations);


	assert(cuboid_structure_.label_cuboids_.size() == num_labels);
	MeshCuboidStructure temp_cuboid_structure(cuboid_structure_);

	// Set label indices.
	LabelIndex main_label_index = 0;
	std::vector<LabelIndex> target_labels;
	target_labels.push_back(1);

	const double variance = 20;
	Eigen::VectorXd pairwise_features_vec;
	double error = 0.0;


	//
	// Ignore height attributes.
	Eigen::VectorXd mask = Eigen::VectorXd::Ones(MeshCuboidJointNormalRelations::k_mat_size);
	unsigned int start_index = 0;

	start_index += 3 * (MeshCuboid::k_num_corners);
	mask.segment(start_index, (MeshCuboid::k_num_corners + 1)).setZero();

	start_index += ((MeshCuboid::k_num_corners + 1) + 3 * (MeshCuboid::k_num_corners + 1));
	mask.segment(start_index, (MeshCuboid::k_num_corners + 1)).setZero();

	start_index += ((MeshCuboid::k_num_corners + 1) + 3 * (MeshCuboid::k_num_corners));
	mask.segment(start_index, (MeshCuboid::k_num_corners + 1)).setZero();

	start_index += ((MeshCuboid::k_num_corners + 1) + 3 * (MeshCuboid::k_num_corners + 1));
	mask.segment(start_index, (MeshCuboid::k_num_corners + 1)).setZero();



	for (std::vector<LabelIndex>::iterator it = target_labels.begin(); it != target_labels.end(); ++it)
	{
		LabelIndex target_label_index = (*it);
		assert(main_label_index < num_labels);
		assert(target_label_index < num_labels);

		MeshCuboidJointNormalRelations *joint_normal_relation = joint_normal_relations[main_label_index][target_label_index];
		assert(joint_normal_relation);
		Eigen::VectorXd mean = joint_normal_relation->get_mean();

		std::array<MyMesh::Point, MeshCuboid::k_num_corners> cuboid_corners_11;
		std::array<MyMesh::Point, MeshCuboid::k_num_corners> cuboid_corners_12;
		for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
		{
			for (unsigned int i = 0; i < 3; ++i)
			{
				unsigned int offset_1 = 0;
				unsigned int offset_2 = (MeshCuboidFeatures::k_num_features - MeshCuboidFeatures::k_corner_index)
					+ MeshCuboidFeatures::k_corner_index;
				cuboid_corners_11[corner_index][i] = mean[offset_1 + 3 * corner_index + i];
				cuboid_corners_12[corner_index][i] = mean[offset_2 + 3 * corner_index + i];
			}
		}

		MeshCuboid *mean_main_cuboid = new MeshCuboid(main_label_index);
		mean_main_cuboid->set_bbox_corners(cuboid_corners_11);
		mean_main_cuboid->update_axes_center_size_corner_points();

		MeshCuboid *mean_target_cuboid = new MeshCuboid(target_label_index);
		mean_target_cuboid->set_bbox_corners(cuboid_corners_12);
		mean_target_cuboid->update_axes_center_size_corner_points();

		MeshCuboidTransformation mean_transformation_1;
		MeshCuboidTransformation mean_transformation_2;
		mean_transformation_1.compute_transformation(mean_main_cuboid);
		mean_transformation_2.compute_transformation(mean_target_cuboid);


		joint_normal_relation->get_pairwise_cuboid_features(mean_main_cuboid, mean_target_cuboid,
			&mean_transformation_1, &mean_transformation_2, pairwise_features_vec);
		Eigen::VectorXd diff = pairwise_features_vec - mean;



		//diff = diff.cwiseProduct(mask);
		//error = diff.transpose() * joint_normal_relation->get_inv_cov() * diff;
		//error = std::exp(-error / variance);
		//std::cout << ">>> error = " << error << std::endl;

		// Add 'error_' member variable to 'MeshCuboid'.
		//mean_main_cuboid->error_ = 1.0;
		//mean_target_cuboid->error_ = 1.0;

		temp_cuboid_structure.label_cuboids_[main_label_index].push_back(mean_main_cuboid);
		temp_cuboid_structure.label_cuboids_[target_label_index].push_back(mean_target_cuboid);
	}
	//


	// For every file in the base path.
	QDir input_dir((FLAGS_data_root_path + FLAGS_mesh_path).c_str());
	assert(input_dir.exists());
	input_dir.setFilter(QDir::Files | QDir::Hidden | QDir::NoSymLinks);
	input_dir.setSorting(QDir::Name);

	std::stringstream output_filename_sstr;
	QDir output_dir;
	output_dir.mkpath(FLAGS_training_dir.c_str());



	std::list<MeshCuboid *> temp_cuboid_2_list;

	QFileInfoList dir_list = input_dir.entryInfoList();
	for (int i = 0; i < dir_list.size(); i++)
	{
		QFileInfo file_info = dir_list.at(i);
		if (file_info.exists() &&
			(file_info.suffix().compare("obj") == 0
			|| file_info.suffix().compare("off") == 0))
		{
			std::string mesh_filepath = std::string(file_info.filePath().toLocal8Bit());
			std::string mesh_name = std::string(file_info.baseName().toLocal8Bit());

			output_filename_sstr.clear(); output_filename_sstr.str("");
			output_filename_sstr << FLAGS_training_dir << std::string("/") << mesh_name << std::string(".arff");

			bool ret = load_object_info(mesh_, cuboid_structure_, mesh_filepath.c_str(), LoadMesh,
				output_filename_sstr.str().c_str(), false);
			if (!ret) continue;


			for (std::vector<LabelIndex>::iterator it = target_labels.begin(); it != target_labels.end(); ++it)
			{
				LabelIndex target_label_index = (*it);

				MeshCuboidJointNormalRelations *joint_normal_relation = joint_normal_relations[main_label_index][target_label_index];
				assert(joint_normal_relation);
				Eigen::VectorXd mean = joint_normal_relation->get_mean();

				if (!cuboid_structure_.label_cuboids_[main_label_index].empty()
					&& !cuboid_structure_.label_cuboids_[target_label_index].empty())
				{
					MeshCuboid *main_cuboid = cuboid_structure_.label_cuboids_[main_label_index].front();
					MeshCuboid *target_cuboid = cuboid_structure_.label_cuboids_[target_label_index].front();

					MeshCuboidTransformation transformation_1;
					MeshCuboidTransformation transformation_2;
					transformation_1.compute_transformation(main_cuboid);
					transformation_2.compute_transformation(target_cuboid);

					MeshCuboidJointNormalRelations::get_pairwise_cuboid_features(main_cuboid, target_cuboid,
						&transformation_1, &transformation_2, pairwise_features_vec);

					std::array<MyMesh::Point, MeshCuboid::k_num_corners> cuboid_corners_11;
					std::array<MyMesh::Point, MeshCuboid::k_num_corners> cuboid_corners_12;
					for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
					{
						for (unsigned int i = 0; i < 3; ++i)
						{
							unsigned int offset_1 = 0;
							unsigned int offset_2 = (MeshCuboidFeatures::k_num_features - MeshCuboidFeatures::k_corner_index)
								+ MeshCuboidFeatures::k_corner_index;
							cuboid_corners_11[corner_index][i] = pairwise_features_vec[offset_1 + 3 * corner_index + i];
							cuboid_corners_12[corner_index][i] = pairwise_features_vec[offset_2 + 3 * corner_index + i];
						}
					}

					MeshCuboid *normalized_target_cuboid = new MeshCuboid(target_label_index);
					normalized_target_cuboid->set_bbox_corners(cuboid_corners_12);
					normalized_target_cuboid->update_axes_center_size_corner_points();

					Eigen::VectorXd diff = pairwise_features_vec - mean;
					diff = diff.cwiseProduct(mask);
					error = diff.transpose() * joint_normal_relation->get_inv_cov() * diff;
					error = std::exp(-error / variance);

					// Add 'error_' member variable to 'MeshCuboid'.
					//normalized_target_cuboid->error_ = error;
					assert(error >= 0);
					assert(error <= 1);
					std::cout << ">>> error = " << error << std::endl;

					temp_cuboid_structure.label_cuboids_[target_label_index].push_back(normalized_target_cuboid);
				}
			}
		}
	}

	cuboid_structure_.clear_cuboids();
	cuboid_structure_ = temp_cuboid_structure;

	draw_cuboid_axes_ = false;
	updateGL();
}

void MeshViewerCore::compute_view_plane_mask_range(const Real _modelview_matrix[16])
{
	unsigned int num_sample_points = cuboid_structure_.num_sample_points();
	assert(num_sample_points > 0);

	std::vector<MyMesh::Point> sample_points(num_sample_points);
	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points;
		++sample_point_index)
	{
		const MeshSamplePoint *sample_point = cuboid_structure_.sample_points_[sample_point_index];
		assert(sample_point);
		sample_points[sample_point_index] = sample_point->point_;
	}

	FLAGS_param_view_plane_mask_min_x = 0;
	FLAGS_param_view_plane_mask_min_y = 0;
	FLAGS_param_view_plane_mask_max_x = 0;
	FLAGS_param_view_plane_mask_max_y = 0;


	Eigen::Matrix3d rotation_mat;
	Eigen::Vector3d translation_vec;

	for (int row = 0; row < 3; ++row)
		for (int col = 0; col < 3; ++col)
			rotation_mat(row, col) = _modelview_matrix[4 * col + row];

	for (int row = 0; row < 3; ++row)
		translation_vec(row) = _modelview_matrix[4 * 3 + row];

	unsigned int num_points = sample_points.size();
	if (num_points == 0)
		return;

	Eigen::MatrixXd view_plane_points_mat = Eigen::MatrixXd(3, num_points);

	Eigen::Vector2d view_plane_bbox_min;
	Eigen::Vector2d view_plane_bbox_max;
	view_plane_bbox_min[0] = view_plane_bbox_min[1] = std::numeric_limits<double>::max();
	view_plane_bbox_max[0] = view_plane_bbox_max[1] = -std::numeric_limits<double>::max();

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_points;
		++sample_point_index)
	{
		Eigen::Vector3d point_vec;
		for (int i = 0; i < 3; ++i)
			point_vec(i) = sample_points[sample_point_index][i];

		// Transformation.
		point_vec = rotation_mat * point_vec + translation_vec;

		view_plane_points_mat.col(sample_point_index) = point_vec;
		for (int i = 0; i < 2; ++i)
		{
			view_plane_bbox_min[i] = std::min(view_plane_bbox_min[i], point_vec[i]);
			view_plane_bbox_max[i] = std::max(view_plane_bbox_max[i], point_vec[i]);
		}
	}

	Eigen::Vector2d view_plane_bbox_center = 0.5 * (view_plane_bbox_min + view_plane_bbox_max);
	Eigen::Vector2d view_plane_bbox_size = view_plane_bbox_max - view_plane_bbox_min;

	if (view_plane_bbox_size[0] == 0 || view_plane_bbox_size[1] == 0)
		return;

	Eigen::Vector2d view_plane_mask_range_center;
	Eigen::Vector2d view_plane_mask_range_initial_size;
	double scale_factor = 1.0;

	// NOTE:
	// Temporary share the same seed with random view point.
	static SimpleRandomCong_t rng_cong;
	simplerandom_cong_seed(&rng_cong, FLAGS_random_view_seed);

	double center_x = static_cast<double>(simplerandom_cong_next(&rng_cong)) / std::numeric_limits<uint32_t>::max();
	double center_y = static_cast<double>(simplerandom_cong_next(&rng_cong)) / std::numeric_limits<uint32_t>::max();
	double size_ratio_x = static_cast<double>(simplerandom_cong_next(&rng_cong)) / std::numeric_limits<uint32_t>::max();
	double size_ratio_y = static_cast<double>(simplerandom_cong_next(&rng_cong)) / std::numeric_limits<uint32_t>::max();

	center_x -= 0.5;
	center_y -= 0.5;
	size_ratio_x += 0.5;
	size_ratio_y += 0.5;

	view_plane_mask_range_center[0] = view_plane_bbox_center[0] - center_x * view_plane_bbox_size[0];
	view_plane_mask_range_center[1] = view_plane_bbox_center[1] - center_y * view_plane_bbox_size[1];
	view_plane_mask_range_initial_size[0] = size_ratio_x * view_plane_bbox_size[0];
	view_plane_mask_range_initial_size[1] = size_ratio_y * view_plane_bbox_size[1];


	for (unsigned int iter = 0; iter < 100; ++iter)
	{
		Eigen::Vector2d view_plane_mask_range_size = scale_factor * view_plane_mask_range_initial_size;
		Eigen::Vector2d view_plane_mask_range_min = view_plane_mask_range_center - view_plane_mask_range_size;
		Eigen::Vector2d view_plane_mask_range_max = view_plane_mask_range_center + view_plane_mask_range_size;
		unsigned int num_removed_points = 0;

		for (SamplePointIndex sample_point_index = 0; sample_point_index < num_points;
			++sample_point_index)
		{
			Eigen::Vector3d sample_point_vec = view_plane_points_mat.col(sample_point_index);
			if (sample_point_vec[0] >= view_plane_mask_range_min[0]
				&& sample_point_vec[1] >= view_plane_mask_range_min[1]
				&& sample_point_vec[0] <= view_plane_mask_range_max[0]
				&& sample_point_vec[1] <= view_plane_mask_range_max[1])
				++num_removed_points;
		}

		Real removed_point_proportion = static_cast<Real>(num_removed_points) / num_points;
		if (std::abs(removed_point_proportion - FLAGS_param_view_plane_mask_proportion) < 0.01)
		{
			break;
		}
		else if (removed_point_proportion < FLAGS_param_view_plane_mask_proportion)
		{
			// Increase the mask size.
			scale_factor += std::pow(0.5, iter + 1);
		}
		else
		{
			// Decrease the mask size.
			scale_factor -= std::pow(0.5, iter + 1);
		}
	}

	Eigen::Vector2d view_plane_mask_range_size = scale_factor * view_plane_mask_range_initial_size;
	Eigen::Vector2d view_plane_mask_range_min = view_plane_mask_range_center - view_plane_mask_range_size;
	Eigen::Vector2d view_plane_mask_range_max = view_plane_mask_range_center + view_plane_mask_range_size;

	FLAGS_param_view_plane_mask_min_x = view_plane_mask_range_min[0];
	FLAGS_param_view_plane_mask_min_y = view_plane_mask_range_min[1];
	FLAGS_param_view_plane_mask_max_x = view_plane_mask_range_max[0];
	FLAGS_param_view_plane_mask_max_y = view_plane_mask_range_max[1];
}

/*
void MeshViewerCore::run_reconstruction_test()
{
	std::string cuboid_filepath = "C:/project/app/cuboid-prediction/experiments/exp1_assembly_chairs/output/33_0.arff";

	// Load basic information.
	bool ret = true;

	ret = ret & cuboid_structure_.load_labels((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_info_filename).c_str());
	ret = ret & cuboid_structure_.load_label_symmetries((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_symmetry_info_filename).c_str());

	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open label information files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	cuboid_structure_.load_cuboids(cuboid_filepath.c_str());
	draw_cuboid_axes_ = true;
	setDrawMode(CUSTOM_VIEW);
	updateGL();

	reconstruct_using_database();

	std::cout << "# of total sample points = " << cuboid_structure_.num_sample_points() << std::endl;
	cuboid_structure_.save_sample_points_to_ply("test");
}

void MeshViewerCore::run_occlusion_test()
{
	assert(occlusion_test_widget_);

	std::vector<Eigen::Vector3f> sample_points;
	sample_points.reserve(cuboid_structure_.num_sample_points());
	for (unsigned point_index = 0; point_index < cuboid_structure_.num_sample_points();
		++point_index)
	{
		MyMesh::Point point = cuboid_structure_.sample_points_[point_index]->point_;
		Eigen::Vector3f point_vec;
		point_vec << point[0], point[1], point[2];
		sample_points.push_back(point_vec);
	}

	occlusion_test_widget_->set_sample_points(sample_points);
	occlusion_test_widget_->set_modelview_matrix(modelview_matrix());

	std::array<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, 4> occlusion_test_result;
	occlusion_test_widget_->get_occlusion_test_result(
		static_cast<float>(FLAGS_occlusion_test_radius), occlusion_test_result);

	occlusion_test_points_.clear();
	occlusion_test_points_.reserve(occlusion_test_result[0].rows() * occlusion_test_result[0].cols());

	for (unsigned int row = 0; row < occlusion_test_result[0].rows(); ++row)
	{
		for (unsigned int col = 0; col < occlusion_test_result[0].cols(); ++col)
		{
			Mesh::Point point;
			for (unsigned int i = 0; i < 3; ++i)
				point[i] = static_cast<Real>(occlusion_test_result[i](row, col));
			Real value = static_cast<Real>(occlusion_test_result[3](row, col));
			assert(value >= 0);
			assert(value <= 1);

			// NOTE:
			// Reverse values.
			// After reversing, the value '1' means that the voxel is perfectly occluded.
			occlusion_test_points_.push_back(std::make_pair(point, 1.0 - value));
		}
	}

	draw_occlusion_test_points_ = true;
	//updateGL();
}

ANNkd_tree* create_occlusion_test_points_kd_tree(
	const std::vector< std::pair<MyMesh::Point, Real> > &_occlusion_test_points,
	const MeshCuboidStructure &_cuboid_structure,
	Eigen::VectorXd &_values,
	ANNpointArray &_ann_points)
{
	Eigen::MatrixXd points = Eigen::MatrixXd(3, _occlusion_test_points.size());
	_values = Eigen::VectorXd(_occlusion_test_points.size());

	for (unsigned point_index = 0; point_index < _occlusion_test_points.size();
		++point_index)
	{
		MyMesh::Point point = _occlusion_test_points[point_index].first;
		Eigen::Vector3d point_vec;
		point_vec << point[0], point[1], point[2];
		points.col(point_index) = point_vec;
		_values(point_index) = _occlusion_test_points[point_index].second;
	}
	
	ANNkd_tree* kd_tree = ICP::create_kd_tree(points, _ann_points);
	return kd_tree;
}
*/
