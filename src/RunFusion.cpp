#include "MeshViewerCore.h"
#include "MeshCuboidParameters.h"
#include "MeshCuboidSolver.h"

#include <sstream>
#include <Eigen/Geometry>
#include <QFileInfo>


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

	if (_sample_label_filepath || _cuboid_filepath)
	{
		std::vector<LabelIndex> sample_point_label_indices = _cuboid_structure.get_sample_point_label_indices();
		assert(sample_point_label_indices.size() == _cuboid_structure.num_sample_points());

		for (SamplePointIndex sample_point_index = 0; sample_point_index < _cuboid_structure.num_sample_points();
			++sample_point_index)
		{
			MeshSamplePoint* sample_point = _cuboid_structure.sample_points_[sample_point_index];
			assert(sample_point);
			LabelIndex label_index = sample_point_label_indices[sample_point_index];
			assert(label_index < _cuboid_structure.num_labels());

			MeshCuboid *cuboid = NULL;
			// NOTE:
			// The current implementation assumes that there is only one part for each label.
			assert(_cuboid_structure.label_cuboids_[label_index].size() <= 1);
			if (!_cuboid_structure.label_cuboids_[label_index].empty())
				cuboid = _cuboid_structure.label_cuboids_[label_index].front();

			assert(cuboid);
			cuboid->add_sample_point(sample_point);
		}
	}

	mesh_.clear_colors();

	return true;
}

void MeshViewerCore::run_test()
{
	setDrawMode(CUSTOM_VIEW);

	bool ret;
	std::string mesh_name("4665");
	std::string exp_root_path = "C:/project/app/cuboid-prediction/experiments/exp5_assembly_chairs/output/";

	std::string mesh_filepath = FLAGS_data_root_path + FLAGS_mesh_path + std::string("/") + mesh_name + std::string(".off");
	std::string sample_filepath = FLAGS_data_root_path + FLAGS_sample_path + std::string("/") + mesh_name + std::string(".pts");

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
		mesh_filepath.c_str(), sample_filepath.c_str(), NULL, NULL);
	MeshCuboidStructure input_points = cuboid_structure_;

	open_modelview_matrix_file(pose_filename.c_str());
	memcpy(snapshot_modelview_matrix, modelview_matrix(), 16 * sizeof(double));

	open_modelview_matrix_file(occlusion_pose_filename.c_str());
	memcpy(occlusion_modelview_matrix, modelview_matrix(), 16 * sizeof(double));


	ret = load_result_info(mesh_, cuboid_structure_,
		mesh_filepath.c_str(), symmetry_sample_filepath.c_str(),
		symmetry_sample_label_filepath.c_str(), cuboid_filepath.c_str());
	MeshCuboidStructure symmetry_reconstruction = cuboid_structure_;

	ret = load_result_info(mesh_, cuboid_structure_,
		mesh_filepath.c_str(), database_sample_filepath.c_str(),
		database_sample_label_filepath.c_str(), cuboid_filepath.c_str());
	MeshCuboidStructure database_reconstruction = cuboid_structure_;

	MyMesh original_mesh;
	MeshCuboidStructure original_cuboid_structure(&original_mesh);
	ret = load_object_info(original_mesh, original_cuboid_structure, mesh_filepath.c_str(), false, false, false);

	updateGL();

	reconstruct_fusion(mesh_filepath.c_str(),
		snapshot_modelview_matrix, occlusion_modelview_matrix,
		original_cuboid_structure, symmetry_reconstruction, database_reconstruction, cuboid_structure_);
	

	std::cout << "[Symmetry] # of points = " << symmetry_reconstruction.num_sample_points() << std::endl;
	std::cout << "[Database] # of points = " << database_reconstruction.num_sample_points() << std::endl;
	std::cout << "[Output] # of points = " << cuboid_structure_.num_sample_points() << std::endl;

	cuboid_structure_.save_sample_points_to_ply("test1");

	set_modelview_matrix(snapshot_modelview_matrix);
	setDrawMode(COLORED_POINT_SAMPLES);
	updateGL();
	snapshot("test_1");
}
