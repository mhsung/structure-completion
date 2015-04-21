#include "MeshViewerCore.h"
#include "MeshCuboidParameters.h"
#include "MeshCuboidSolver.h"

#include <sstream>
#include <Eigen/Geometry>
#include <QFileInfo>


typedef OpenMesh::Vec3i VoxelIndex3D;
class VoxelGrid
{
public:
	VoxelGrid(MyMesh::Point _min, MyMesh::Point _max, Real _unit_size)
		: min_(_min)
		, max_(_max)
	{
		MyMesh::Normal diff = max_ - min_;
		for (unsigned int i = 0; i < 3; ++i)
		{
			assert(diff[i] > 0);
			n_voxels_[i] = static_cast<int>(std::round(diff[i] / _unit_size));
			n_voxels_[i] = std::max(n_voxels_[i], 2);
		}
	}

	int n_voxels() const
	{
		return (n_voxels_[0] * n_voxels_[1] * n_voxels_[2]);
	}

	int n_axis_voxels(const int _axis_index) const
	{
		assert(_axis_index < 3);
		return n_voxels_[_axis_index];
	}

	int get_voxel_index(VoxelIndex3D _xyz_index) const
	{
		assert(_xyz_index[0] >= 0);
		assert(_xyz_index[1] >= 0);
		assert(_xyz_index[2] >= 0);
		assert(_xyz_index[0] < n_voxels_[0]);
		assert(_xyz_index[1] < n_voxels_[1]);
		assert(_xyz_index[2] < n_voxels_[2]);

		int voxel_index =
			_xyz_index[0] * (n_voxels_[1] * n_voxels_[2])
			+ _xyz_index[1] * n_voxels_[2]
			+ _xyz_index[2];

		assert(voxel_index >= 0);
		assert(voxel_index < n_voxels());

		// DEBUG.
		VoxelIndex3D test_xyz_index = get_voxel_index(voxel_index);
		assert(_xyz_index[0] == test_xyz_index[0]);
		assert(_xyz_index[1] == test_xyz_index[1]);
		assert(_xyz_index[2] == test_xyz_index[2]);

		return voxel_index;
	}

	VoxelIndex3D get_voxel_index(const int _voxel_index) const
	{
		assert(_voxel_index >= 0);
		assert(_voxel_index < n_voxels());
		VoxelIndex3D xyz_index;

		int temp = _voxel_index;
		xyz_index[0] = temp / (n_voxels_[1] * n_voxels_[2]);
		temp -= (xyz_index[0] * n_voxels_[1] * n_voxels_[2]);

		xyz_index[1] = temp / n_voxels_[2];
		temp -= (xyz_index[1] * n_voxels_[2]);

		xyz_index[2] = temp;

		assert(xyz_index[0] >= 0);
		assert(xyz_index[1] >= 0);
		assert(xyz_index[2] >= 0);
		assert(xyz_index[0] < n_voxels_[0]);
		assert(xyz_index[1] < n_voxels_[1]);
		assert(xyz_index[2] < n_voxels_[2]);
		return xyz_index;
	}

	int get_voxel_index(const MyMesh::Point _point) const
	{
		MyMesh::Normal diff = max_ - min_;
		VoxelIndex3D xyz_index;

		for (unsigned int i = 0; i < 3; ++i)
		{
			assert(diff[i] > 0);
			xyz_index[i] = static_cast<int>((_point[i] - min_[i]) / diff[i] * n_voxels_[i]);
			xyz_index[i] = std::max(xyz_index[i], 0);
			xyz_index[i] = std::min(xyz_index[i], n_voxels_[i] - 1);
		}

		return get_voxel_index(xyz_index);
	}

	MyMesh::Point get_center(const int _voxel_index) const
	{
		VoxelIndex3D xyz_index = get_voxel_index(_voxel_index);
		MyMesh::Normal diff = max_ - min_;
		MyMesh::Point center;

		for (unsigned int i = 0; i < 3; ++i)
		{
			assert(diff[i] > 0);
			Real ratio = static_cast<Real>(xyz_index[i] + 0.5) / n_voxels_[i];
			center[i] = min_[i] + ratio * diff[i];
		}

		return center;
	}

	void get_centers(std::vector<MyMesh::Point> &_centers) const
	{
		_centers.clear();
		_centers.resize(n_voxels());
		for (int voxel_index = 0; voxel_index < n_voxels(); ++voxel_index)
			_centers[voxel_index] = get_center(voxel_index);
	}

	void get_point_correspondences(
		const std::vector<MeshSamplePoint *> &_points,
		std::vector<int> &_points_to_voxels,
		std::vector< std::list<int> > &_voxels_to_points)
	{
		const unsigned int num_points = _points.size();
		_points_to_voxels.clear();
		_voxels_to_points.clear();
		_points_to_voxels.resize(num_points);
		_voxels_to_points.resize(n_voxels());

		for (unsigned int point_index = 0; point_index < num_points; ++point_index)
		{
			assert(_points[point_index]);
			int voxel_index = get_voxel_index(_points[point_index]->point_);;
			assert(voxel_index < n_voxels());
			_points_to_voxels[point_index] = voxel_index;
			_voxels_to_points[voxel_index].push_back(point_index);
		}
	}

private:
	MyMesh::Point min_;
	MyMesh::Point max_;
	VoxelIndex3D n_voxels_;
};

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
	std::string mesh_name("4064");
	std::string exp_root_path = "C:/project/app/cuboid-prediction/experiments/exp1_assembly_chairs/output/";

	std::string mesh_filepath = FLAGS_data_root_path + FLAGS_mesh_path + std::string("/") + mesh_name + std::string(".off");
	std::string sample_filepath = FLAGS_data_root_path + FLAGS_sample_path + std::string("/") + mesh_name + std::string(".pts");
	std::string cuboid_filepath = exp_root_path + std::string("/") + mesh_name + std::string("_0.arff");

	std::string symmetry_sample_filepath = exp_root_path + std::string("/") + mesh_name + std::string("_0_symmetry.pts");
	std::string symmetry_sample_label_filepath = exp_root_path + std::string("/") + mesh_name + std::string("_0_symmetry_label.arff");
	
	std::string database_sample_filepath = exp_root_path + std::string("/") + mesh_name + std::string("_0_database.pts");
	std::string database_sample_label_filepath = exp_root_path + std::string("/") + mesh_name + std::string("_0_database_label.arff");


	double snapshot_modelview_matrix[16];
	double occlusion_modelview_matrix[16];

	ret = load_result_info(mesh_, cuboid_structure_,
		mesh_filepath.c_str(), sample_filepath.c_str(), NULL, NULL);
	MeshCuboidStructure input_points = cuboid_structure_;

	open_modelview_matrix_file(FLAGS_pose_filename.c_str());
	memcpy(snapshot_modelview_matrix, modelview_matrix(), 16 * sizeof(double));

	if (FLAGS_occlusion_pose_filename == "")
		set_random_view_direction(true);
	else
		open_modelview_matrix_file(FLAGS_occlusion_pose_filename.c_str());
	memcpy(occlusion_modelview_matrix, modelview_matrix(), 16 * sizeof(double));

	set_modelview_matrix(snapshot_modelview_matrix);


	ret = load_result_info(mesh_, cuboid_structure_,
		mesh_filepath.c_str(), symmetry_sample_filepath.c_str(),
		symmetry_sample_label_filepath.c_str(), cuboid_filepath.c_str());
	MeshCuboidStructure symmetry_reconstruction = cuboid_structure_;

	ret = load_result_info(mesh_, cuboid_structure_,
		mesh_filepath.c_str(), database_sample_filepath.c_str(),
		database_sample_label_filepath.c_str(), cuboid_filepath.c_str());
	MeshCuboidStructure database_reconstruction = cuboid_structure_;

	cuboid_structure_.clear_sample_points();
	const Real radius = FLAGS_param_observed_point_radius
		* cuboid_structure_.mesh_->get_object_diameter();
	

	assert(symmetry_reconstruction.num_labels() == database_reconstruction.num_labels());
	unsigned int num_labels = symmetry_reconstruction.num_labels();

	for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
	{
		if (symmetry_reconstruction.label_cuboids_[label_index].empty()
			|| database_reconstruction.label_cuboids_[label_index].empty()
			|| cuboid_structure_.label_cuboids_[label_index].empty())
			continue;

		MeshCuboid *symmetry_cuboid = symmetry_reconstruction.label_cuboids_[label_index].front();
		MeshCuboid *database_cuboid = database_reconstruction.label_cuboids_[label_index].front();
		MeshCuboid *output_cuboid = cuboid_structure_.label_cuboids_[label_index].front();
		assert(symmetry_cuboid);
		assert(database_cuboid);
		assert(output_cuboid);

		if (symmetry_cuboid->num_sample_points() == 0
			|| database_cuboid->num_sample_points() == 0)
			continue;


		// 1. Per-part ICP.
		Eigen::MatrixXd symmetry_sample_points(3, symmetry_cuboid->num_sample_points());
		Eigen::MatrixXd database_sample_points(3, database_cuboid->num_sample_points());

		for (SamplePointIndex sample_point_index = 0; sample_point_index < symmetry_cuboid->num_sample_points();
			++sample_point_index)
		{
			assert(symmetry_cuboid->get_sample_point(sample_point_index));
			MyMesh::Point point = symmetry_cuboid->get_sample_point(sample_point_index)->point_;
			for (unsigned int i = 0; i < 3; ++i)
				symmetry_sample_points.col(sample_point_index)(i) = point[i];
		}

		for (SamplePointIndex sample_point_index = 0; sample_point_index < database_cuboid->num_sample_points();
			++sample_point_index)
		{
			assert(database_cuboid->get_sample_point(sample_point_index));
			MyMesh::Point point = database_cuboid->get_sample_point(sample_point_index)->point_;
			for (unsigned int i = 0; i < 3; ++i)
				database_sample_points.col(sample_point_index)(i) = point[i];
		}

		const Real neighbor_distance = FLAGS_param_sample_point_neighbor_distance * mesh_.get_object_diameter();

		Eigen::Matrix3d rotation_mat;
		Eigen::Vector3d translation_vec;
		double icp_error = ICP::run_iterative_closest_points(database_sample_points, symmetry_sample_points,
			rotation_mat, translation_vec, &neighbor_distance);
		//std::cout << "ICP Error = " << icp_error << std::endl;

		for (SamplePointIndex sample_point_index = 0; sample_point_index < database_cuboid->num_sample_points();
			++sample_point_index)
		{
			MeshSamplePoint* database_sample_point = database_cuboid->get_sample_point(sample_point_index);
			assert(database_sample_point);

			for (unsigned int i = 0; i < 3; ++i)
				database_sample_point->point_[i] = database_sample_points.col(sample_point_index)[i];
		}


		// 2. Voxel visibility test.
		Eigen::Vector3d bbox_min_vec = symmetry_sample_points.rowwise().minCoeff();
		Eigen::Vector3d bbox_max_vec = symmetry_sample_points.rowwise().maxCoeff();

		MyMesh::Point bbox_min, bbox_max;
		for (unsigned int i = 0; i < 3; ++i)
		{
			bbox_min[i] = bbox_min_vec[i];
			bbox_max[i] = bbox_max_vec[i];
		}

		VoxelGrid voxels(bbox_min, bbox_max, 0.5 * radius);

		std::vector<MyMesh::Point> voxel_centers;
		std::vector<Real> voxel_visibility_values;
		voxels.get_centers(voxel_centers);

		std::cout << "Compute visibility values... ";
		MeshCuboid::compute_cuboid_surface_point_visibility(
			occlusion_modelview_matrix, radius, input_points.sample_points_,
			voxel_centers, NULL, voxel_visibility_values);
		std::cout << "Done." << std::endl;


		// 3. Voxel filling.
		std::cout << "Filling voxels... ";
		std::vector<int> symmetry_points_to_voxels;
		std::vector< std::list<int> > symmetry_voxels_to_points;
		voxels.get_point_correspondences(symmetry_cuboid->get_sample_points(),
			symmetry_points_to_voxels, symmetry_voxels_to_points);

		std::vector<int> database_points_to_voxels;
		std::vector< std::list<int> > database_voxels_to_points;
		voxels.get_point_correspondences(database_cuboid->get_sample_points(),
			database_points_to_voxels, database_voxels_to_points);

		for (unsigned int voxel_index = 0; voxel_index < voxels.n_voxels(); ++voxel_index)
		{
			assert(voxel_index < symmetry_voxels_to_points.size());
			assert(voxel_index < database_voxels_to_points.size());
			VoxelIndex3D xyz_index = voxels.get_voxel_index(voxel_index);

			const std::list<int>& symmetry_voxel_point_indices = symmetry_voxels_to_points[voxel_index];
			const std::list<int>& database_voxel_point_indices = database_voxels_to_points[voxel_index];

			//if (voxel_visibility_values[voxel_index] > 0.5)
			{
				for (std::list<int>::const_iterator it = symmetry_voxel_point_indices.begin();
					it != symmetry_voxel_point_indices.end(); ++it)
				{
					MeshSamplePoint *sample_point = symmetry_cuboid->get_sample_point(*it);
					assert(sample_point);

					MeshSamplePoint *new_sample_point = cuboid_structure_.add_sample_point(
						sample_point->point_, sample_point->normal_);
					new_sample_point->error_ = 0.0;
					output_cuboid->add_sample_point(new_sample_point);
				}
			}
			//else
			{
				for (std::list<int>::const_iterator it = database_voxel_point_indices.begin();
					it != database_voxel_point_indices.end(); ++it)
				{
					MeshSamplePoint *sample_point = database_cuboid->get_sample_point(*it);
					assert(sample_point);

					MeshSamplePoint *new_sample_point = cuboid_structure_.add_sample_point(
						sample_point->point_, sample_point->normal_);
					new_sample_point->error_ = 1.0;
					output_cuboid->add_sample_point(new_sample_point);
				}
			}
		}
		std::cout << "Done." << std::endl;
	}

	std::cout << "[Symmetry] # of points = " << symmetry_reconstruction.num_sample_points() << std::endl;
	std::cout << "[Database] # of points = " << database_reconstruction.num_sample_points() << std::endl;
	std::cout << "[Output] # of points = " << cuboid_structure_.num_sample_points() << std::endl;

	cuboid_structure_.save_sample_points_to_ply("test1");
	database_reconstruction.save_sample_points_to_ply("test2");

	setDrawMode(COLORED_POINT_SAMPLES);
	updateGL();
	snapshot("test_1");

	//cuboid_structure_.clear_cuboids();
	//setDrawMode(CUSTOM_VIEW);
	//updateGL();
	//snapshot("test_2");
}