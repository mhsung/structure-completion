#include "MeshViewerCore.h"
#include "MeshCuboidFusion.h"
#include "MeshCuboidParameters.h"
#include "MeshCuboidSolver.h"
//#include "MeshCuboidNonLinearSolver.h"

#include <sstream>
#include <Eigen/Geometry>
#include <QFileInfo>


void MeshViewerCore::test_initialize()
{
	setDrawMode(CUSTOM_VIEW);
	draw_cuboid_axes_ = true;

	cuboid_structure_.test_load_cuboids("cuboids.csv");

	unsigned int num_all_cuboid_labels = 2;


	for (LabelIndex label_index_1 = 0; label_index_1 < test_joint_normal_relations_.size(); ++label_index_1)
		for (LabelIndex label_index_2 = 0; label_index_2 < test_joint_normal_relations_[label_index_1].size(); ++label_index_2)
			delete test_joint_normal_relations_[label_index_1][label_index_2];

	test_joint_normal_relations_.clear();
	test_joint_normal_relations_.resize(num_all_cuboid_labels);
	for (LabelIndex label_index_1 = 0; label_index_1 < num_all_cuboid_labels; ++label_index_1)
	{
		test_joint_normal_relations_[label_index_1].resize(num_all_cuboid_labels, NULL);
		for (LabelIndex label_index_2 = 0; label_index_2 < num_all_cuboid_labels; ++label_index_2)
		{
			if (label_index_1 == label_index_2)
				continue;

			std::stringstream relation_filename_sstr;
			relation_filename_sstr << std::string("joint_normal_")
				<< label_index_1 << std::string("_")
				<< label_index_2 << std::string(".dat");

			QFileInfo relation_file(relation_filename_sstr.str().c_str());
			if (!relation_file.exists()) continue;

			test_joint_normal_relations_[label_index_1][label_index_2] = new MeshCuboidJointNormalRelations();
			bool ret = test_joint_normal_relations_[label_index_1][label_index_2]->load_joint_normal_dat(
				relation_filename_sstr.str().c_str());
			if (!ret)
			{
				do {
					std::cout << '\n' << "Press the Enter key to continue.";
				} while (std::cin.get() != '\n');
			}
		}
	}

	delete test_joint_normal_predictor_;
	test_joint_normal_predictor_ = new MeshCuboidJointNormalRelationPredictor(test_joint_normal_relations_);

	all_cuboids_.clear();
	for (std::vector< std::vector<MeshCuboid *> >::iterator it = cuboid_structure_.label_cuboids_.begin();
		it != cuboid_structure_.label_cuboids_.end(); ++it)
		all_cuboids_.insert(all_cuboids_.end(), (*it).begin(), (*it).end());

	open_modelview_matrix_file("occlusion_pose.txt");
	set_view_direction();
	memcpy(test_occlusion_modelview_matrix_, modelview_matrix(), 16 * sizeof(double));
	open_modelview_matrix_file("pose.txt");

	update_cuboid_surface_points(cuboid_structure_, test_occlusion_modelview_matrix_);
	updateGL();
}

void MeshViewerCore::test_translate(const MyMesh::Normal _translation)
{
	Eigen::Vector3d translation_vec;
	for (unsigned int i = 0; i < 3; ++i)
		translation_vec[i] = 0.001 * _translation[i];

	assert(!all_cuboids_.empty());
	all_cuboids_[0]->translate(translation_vec);

	update_cuboid_surface_points(cuboid_structure_, test_occlusion_modelview_matrix_);
	
	double single_total_energy = 0, pair_total_energy = 0;
	get_optimization_error(all_cuboids_, *test_joint_normal_predictor_, single_total_energy, pair_total_energy);
	std::cout << "Energy: (pair = " << pair_total_energy
		<< ", single = " << single_total_energy << ")" << std::endl;

	updateGL();
}

void MeshViewerCore::test_rotate(const Real _angle)
{
	Eigen::AngleAxisd rotation(_angle / 180.0 * M_PI, Eigen::Vector3d::UnitZ());

	assert(!all_cuboids_.empty());
	all_cuboids_[0]->rotate(rotation.toRotationMatrix(), false);
	all_cuboids_[0]->update_corner_points();

	Eigen::Matrix3d rotation_mat;
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		for (unsigned int i = 0; i < 3; ++i)
			rotation_mat.row(axis_index)(i) = all_cuboids_[0]->get_bbox_axis(axis_index)[i];

	rotation.fromRotationMatrix(rotation_mat);
	std::cout << "Angle = " << rotation.angle() / M_PI * 180.0 << std::endl;

	update_cuboid_surface_points(cuboid_structure_, test_occlusion_modelview_matrix_);

	double single_total_energy = 0, pair_total_energy = 0;
	get_optimization_error(all_cuboids_, *test_joint_normal_predictor_, single_total_energy, pair_total_energy);
	std::cout << "Energy: (pair = " << pair_total_energy
		<< ", single = " << single_total_energy << ")" << std::endl;

	updateGL();
}

void MeshViewerCore::test_scale(const Real _scale_x, const Real _scale_y)
{
	assert(!all_cuboids_.empty());

	MyMesh::Point bbox_center = all_cuboids_[0]->get_bbox_center();
	MyMesh::Normal bbox_size = all_cuboids_[0]->get_bbox_size();

	bbox_center[0] += (0.001 * _scale_x);
	bbox_center[1] += (0.001 * _scale_y);

	bbox_size[0] += (2 * 0.001 * _scale_x);
	bbox_size[1] += (2 * 0.001 * _scale_y);

	all_cuboids_[0]->set_bbox_center(bbox_center);
	all_cuboids_[0]->set_bbox_size(bbox_size, true);

	update_cuboid_surface_points(cuboid_structure_, test_occlusion_modelview_matrix_);

	double single_total_energy = 0, pair_total_energy = 0;
	get_optimization_error(all_cuboids_, *test_joint_normal_predictor_, single_total_energy, pair_total_energy);
	std::cout << "Energy: (pair = " << pair_total_energy
		<< ", single = " << single_total_energy << ")" << std::endl;

	updateGL();
}

void MeshViewerCore::test_optimize()
{
	optimize_attributes(cuboid_structure_, test_occlusion_modelview_matrix_, *test_joint_normal_predictor_,
		FLAGS_param_opt_single_energy_term_weight, FLAGS_param_opt_symmetry_energy_term_weight,
		FLAGS_param_opt_max_iterations, "log.txt", this, false);
}



bool MeshViewerCore::load_result_info(
	MyMesh &_mesh,
	MeshCuboidStructure &_cuboid_structure,
	const char* _mesh_filepath,
	const char* _sample_filepath,
	const char* _sample_label_filepath,
	const char* _cuboid_filepath)
{
	bool ret = true;
	ret = ret & _cuboid_structure.load_labels((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_info_filename).c_str());
	ret = ret & _cuboid_structure.load_label_symmetries((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_label_symmetry_info_filename).c_str());

	// Load symmetry groups.
	ret = ret & _cuboid_structure.load_symmetry_groups((FLAGS_data_root_path +
		FLAGS_label_info_path + FLAGS_symmetry_group_info_filename).c_str());

	if (!ret)
	{
		do {
			std::cout << "Error: Cannot open label information files.";
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	QFileInfo file_info(_mesh_filepath);
	std::string mesh_name(file_info.baseName().toLocal8Bit());

	QFileInfo mesh_file(_mesh_filepath);
	QFileInfo sample_file(_sample_filepath);

	QFileInfo cuboid_file(_cuboid_filepath);

	if (!mesh_file.exists())
	{
		std::cerr << "Error: The mesh file does not exist (" << _mesh_filepath << ")." << std::endl;
		return false;
	}
	if (!sample_file.exists())
	{
		std::cerr << "Error: The sample file does not exist (" << _sample_filepath << ")." << std::endl;
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
		ret = _mesh.open_mesh(_mesh_filepath);
		assert(ret);
	}

	ret = _cuboid_structure.load_sample_points(_sample_filepath, false);
	assert(ret);
	assert(_cuboid_structure.num_sample_points() > 0);

	if (_sample_label_filepath)
	{
		QFileInfo sample_label_file(_sample_label_filepath);
		if (!sample_label_file.exists())
		{
			std::cerr << "Error: The sample label file does not exist (" << _sample_label_filepath << ")." << std::endl;
			return false;
		}

		ret = _cuboid_structure.load_sample_point_labels(_sample_label_filepath);
		assert(ret);
	}

	if (_cuboid_filepath)
	{
		if (!cuboid_file.exists())
		{
			std::cerr << "Error: The cuboid file does not exist (" << _cuboid_filepath << ")." << std::endl;
			return false;
		}

		ret = _cuboid_structure.load_cuboids(_cuboid_filepath, false);
		assert(ret);

		std::vector<MeshCuboid *> all_cuboids = _cuboid_structure.get_all_cuboids();
		for (std::vector<MeshCuboid *>::iterator it = all_cuboids.begin(); it != all_cuboids.end(); ++it)
		{
			MeshCuboid *cuboid = (*it);
			cuboid->clear_sample_points();
		}
	}

	if (_sample_label_filepath && _cuboid_filepath)
	{
		std::vector<LabelIndex> sample_point_label_indices;
		_cuboid_structure.get_sample_point_label_indices_from_confidences(sample_point_label_indices);
		assert(sample_point_label_indices.size() == _cuboid_structure.num_sample_points());

		for (SamplePointIndex sample_point_index = 0; sample_point_index < _cuboid_structure.num_sample_points();
			++sample_point_index)
		{
			MeshSamplePoint* sample_point = _cuboid_structure.sample_points_[sample_point_index];
			assert(sample_point);
			LabelIndex label_index = sample_point_label_indices[sample_point_index];
			if (label_index >= _cuboid_structure.num_labels())
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

void MeshViewerCore::run_fusion_test()
{
	setDrawMode(CUSTOM_VIEW);

	bool ret;
	std::string mesh_name("3");
	std::string exp_root_path = "C:/project/app/cuboid-prediction/experiments/test_assembly_chairs/output/";

	std::string mesh_filepath = FLAGS_data_root_path + FLAGS_mesh_path + std::string("/") + mesh_name + std::string(".off");
	std::string sample_filepath = FLAGS_data_root_path + FLAGS_sample_path + std::string("/") + mesh_name + std::string(".pts");
	std::string sample_label_filepath = FLAGS_data_root_path + FLAGS_sample_label_path + std::string("/") + mesh_name + std::string(".arff");

	std::string cuboid_filepath = exp_root_path + std::string("/") + mesh_name
		+ std::string("/") + mesh_name + std::string("_0.arff");

	std::string symmetry_sample_filepath = exp_root_path + std::string("/") + mesh_name
		+ std::string("/") + mesh_name + std::string("_0_symmetry.pts");
	std::string symmetry_sample_label_filepath = exp_root_path + std::string("/") + mesh_name
		+ std::string("/") + mesh_name + std::string("_0_symmetry_label.arff");

	std::string database_sample_filepath = exp_root_path + std::string("/") + mesh_name
		+ std::string("/") + mesh_name + std::string("_0_database.pts");
	std::string database_sample_label_filepath = exp_root_path + std::string("/") + mesh_name
		+ std::string("/") + mesh_name + std::string("_0_database_label.arff");

	std::string pose_filename = exp_root_path + std::string("../pose.txt");
	std::string occlusion_pose_filename = exp_root_path + std::string("../occlusion_pose.txt");


	double snapshot_modelview_matrix[16];
	double occlusion_modelview_matrix[16];

	ret = load_result_info(mesh_, cuboid_structure_,
		mesh_filepath.c_str(), sample_filepath.c_str(),
		sample_label_filepath.c_str(), cuboid_filepath.c_str());
	updateGL();

	open_modelview_matrix_file(pose_filename.c_str());
	memcpy(snapshot_modelview_matrix, modelview_matrix(), 16 * sizeof(double));

	open_modelview_matrix_file(occlusion_pose_filename.c_str());
	memcpy(occlusion_modelview_matrix, modelview_matrix(), 16 * sizeof(double));


	MyMesh symmetry_mesh;
	MeshCuboidStructure symmetry_cuboid_structure(&symmetry_mesh);
	ret = load_result_info(symmetry_mesh, symmetry_cuboid_structure,
		mesh_filepath.c_str(), symmetry_sample_filepath.c_str(),
		symmetry_sample_label_filepath.c_str(), cuboid_filepath.c_str());

	MyMesh database_mesh;
	MeshCuboidStructure database_cuboid_structure(&database_mesh);
	ret = load_result_info(database_mesh, database_cuboid_structure,
		mesh_filepath.c_str(), database_sample_filepath.c_str(),
		database_sample_label_filepath.c_str(), cuboid_filepath.c_str());

	MeshCuboidStructure original_cuboid_structure = cuboid_structure_;

	reconstruct_fusion(mesh_filepath.c_str(),
		snapshot_modelview_matrix, occlusion_modelview_matrix,
		original_cuboid_structure, symmetry_cuboid_structure,
		database_cuboid_structure, cuboid_structure_);


	std::cout << "[Symmetry] # of points = " << symmetry_cuboid_structure.num_sample_points() << std::endl;
	std::cout << "[Database] # of points = " << database_cuboid_structure.num_sample_points() << std::endl;
	std::cout << "[Output] # of points = " << cuboid_structure_.num_sample_points() << std::endl;

	cuboid_structure_.save_sample_points_to_ply("test1");

	set_modelview_matrix(snapshot_modelview_matrix);
	setDrawMode(COLORED_POINT_SAMPLES);
	updateGL();
	snapshot("test_1");
}
