#include "MeshCuboidFusion.h"

#include "MeshCuboidParameters.h"
#include "ICP.h"
#include "Utilities.h"

#include <Eigen/Core>
#include <MRFEnergy.h>


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

void reconstruct_fusion_simple(
	const MeshCuboidStructure &_symmetry_reconstruction,
	const MeshCuboidStructure &_database_reconstruction,
	MeshCuboidStructure &_output_cuboid_structure)
{
	assert(_symmetry_reconstruction.num_labels() == _database_reconstruction.num_labels());

	const Real neighbor_distance = FLAGS_param_sample_point_neighbor_distance
		* _output_cuboid_structure.mesh_->get_object_diameter();

	_output_cuboid_structure = _symmetry_reconstruction;
	unsigned int num_labels = _symmetry_reconstruction.num_labels();

	for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
	{
		if (_symmetry_reconstruction.label_cuboids_[label_index].empty()
			|| _database_reconstruction.label_cuboids_[label_index].empty()
			|| _output_cuboid_structure.label_cuboids_[label_index].empty())
			continue;

		MeshCuboid *symmetry_cuboid = _symmetry_reconstruction.label_cuboids_[label_index].front();
		MeshCuboid *database_cuboid = _database_reconstruction.label_cuboids_[label_index].front();
		MeshCuboid *output_cuboid = _output_cuboid_structure.label_cuboids_[label_index].front();
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

		Eigen::Matrix3d rotation_mat;
		Eigen::Vector3d translation_vec;
		double icp_error = ICP::run_iterative_closest_points(database_sample_points, symmetry_sample_points,
			rotation_mat, translation_vec, &neighbor_distance);
		//std::cout << "ICP Error = " << icp_error << std::endl;

		for (SamplePointIndex sample_point_index = 0; sample_point_index < database_cuboid->num_sample_points();
			++sample_point_index)
		{
			MeshSamplePoint *sample_point = database_cuboid->get_sample_point(sample_point_index);
			assert(sample_point);

			MyMesh::Point new_point;
			for (unsigned int i = 0; i < 3; ++i)
				new_point[i] = database_sample_points.col(sample_point_index)[i];

			MeshSamplePoint *new_sample_point = _output_cuboid_structure.add_sample_point(
				new_point, sample_point->normal_);
			output_cuboid->add_sample_point(new_sample_point);
		}
	}

	/*
	ret = cuboid_structure_.load_dense_sample_points(dense_sample_filepath.c_str());
	assert(ret);
	set_modelview_matrix(_occlusion_modelview_matrix, false);
	remove_occluded_points();
	set_modelview_matrix(_snapshot_modelview_matrix);

	cuboid_structure_.copy_sample_points_to_symmetric_position();

	std::vector<LabelIndex> reconstructed_label_indices;
	reconstructed_label_indices.push_back(0);
	reconstructed_label_indices.push_back(1);
	reconstruct_using_database(&reconstructed_label_indices);
	*/
}

void reconstruct_fusion(const char *_mesh_filepath,
	const double *_snapshot_modelview_matrix,
	const double *_occlusion_modelview_matrix,
	const MeshCuboidStructure &_original_cuboid_structure,
	const MeshCuboidStructure &_symmetry_reconstruction,
	const MeshCuboidStructure &_database_reconstruction,
	MeshCuboidStructure &_output_cuboid_structure)
{
	assert(_symmetry_reconstruction.num_labels() == _database_reconstruction.num_labels());

	const Real neighbor_distance = FLAGS_param_sample_point_neighbor_distance
		* _output_cuboid_structure.mesh_->get_object_diameter();
	const Real occlusion_radius = FLAGS_param_occlusion_test_neighbor_distance
		* _output_cuboid_structure.mesh_->get_object_diameter();


	_output_cuboid_structure.clear_sample_points();

	// Check whether a sample point in the symmetry reconstruction is visited or not.
	bool *is_symmetry_point_visited = new bool[_symmetry_reconstruction.num_sample_points()];
	memset(is_symmetry_point_visited, false, _symmetry_reconstruction.num_sample_points() * sizeof(bool));

	unsigned int num_labels = _symmetry_reconstruction.num_labels();
	for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
	{
		if (_symmetry_reconstruction.label_cuboids_[label_index].empty()
			|| _database_reconstruction.label_cuboids_[label_index].empty()
			|| _output_cuboid_structure.label_cuboids_[label_index].empty())
			continue;

		MeshCuboid *symmetry_cuboid = _symmetry_reconstruction.label_cuboids_[label_index].front();
		MeshCuboid *database_cuboid = _database_reconstruction.label_cuboids_[label_index].front();
		MeshCuboid *output_cuboid = _output_cuboid_structure.label_cuboids_[label_index].front();
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
			MeshSamplePoint *sample_point = symmetry_cuboid->get_sample_point(sample_point_index);
			assert(sample_point);
			MyMesh::Point point = symmetry_cuboid->get_sample_point(sample_point_index)->point_;
			for (unsigned int i = 0; i < 3; ++i)
				symmetry_sample_points.col(sample_point_index)(i) = point[i];

			//
			is_symmetry_point_visited[sample_point->sample_point_index_] = true;
			//
		}

		for (SamplePointIndex sample_point_index = 0; sample_point_index < database_cuboid->num_sample_points();
			++sample_point_index)
		{
			MeshSamplePoint *sample_point = database_cuboid->get_sample_point(sample_point_index);
			assert(sample_point);
			MyMesh::Point point = database_cuboid->get_sample_point(sample_point_index)->point_;
			for (unsigned int i = 0; i < 3; ++i)
				database_sample_points.col(sample_point_index)(i) = point[i];
		}

		Eigen::Matrix3d rotation_mat;
		Eigen::Vector3d translation_vec;
		double icp_error = ICP::run_iterative_closest_points(database_sample_points, symmetry_sample_points,
			rotation_mat, translation_vec, &neighbor_distance);

		if (translation_vec.norm() > 2 * neighbor_distance)
			icp_error = -1;

		// If ICP failed, recover the original database point cloud.
		if (icp_error < 0)
		{
			for (SamplePointIndex sample_point_index = 0; sample_point_index < database_cuboid->num_sample_points();
				++sample_point_index)
			{
				MeshSamplePoint *sample_point = database_cuboid->get_sample_point(sample_point_index);
				assert(sample_point);
				MyMesh::Point point = database_cuboid->get_sample_point(sample_point_index)->point_;
				for (unsigned int i = 0; i < 3; ++i)
					database_sample_points.col(sample_point_index)(i) = point[i];
			}
		}


		for (SamplePointIndex sample_point_index = 0; sample_point_index < database_cuboid->num_sample_points();
			++sample_point_index)
		{
			MeshSamplePoint* database_sample_point = database_cuboid->get_sample_point(sample_point_index);
			assert(database_sample_point);

			for (unsigned int i = 0; i < 3; ++i)
				database_sample_point->point_[i] = database_sample_points.col(sample_point_index)[i];
		}


		// 2. Voxel visibility test.
		Eigen::Vector3d symmetry_bbox_min_vec = symmetry_sample_points.rowwise().minCoeff();
		Eigen::Vector3d symmetry_bbox_max_vec = symmetry_sample_points.rowwise().maxCoeff();
		Eigen::Vector3d database_bbox_min_vec = database_sample_points.rowwise().minCoeff();
		Eigen::Vector3d database_bbox_max_vec = database_sample_points.rowwise().maxCoeff();

		MyMesh::Point bbox_min, bbox_max;
		for (unsigned int i = 0; i < 3; ++i)
		{
			bbox_min[i] = std::min(symmetry_bbox_min_vec[i], database_bbox_min_vec[i]);
			bbox_max[i] = std::max(symmetry_bbox_max_vec[i], database_bbox_max_vec[i]);
		}

		VoxelGrid voxels(bbox_min, bbox_max, occlusion_radius);

		std::vector<MyMesh::Point> voxel_centers;
		std::vector<Real> voxel_visibility_values;
		voxels.get_centers(voxel_centers);

		std::cout << "Compute visibility values... ";
		MeshCuboid::compute_cuboid_surface_point_visibility(
			_occlusion_modelview_matrix, occlusion_radius, _original_cuboid_structure.sample_points_,
			voxel_centers, NULL, voxel_visibility_values);
		std::cout << "Done." << std::endl;


		// 3. Denoising (MRF).
		MRFEnergy<TypeBinary>* mrf;
		MRFEnergy<TypeBinary>::NodeId* nodes;
		MRFEnergy<TypeBinary>::Options options;
		TypeBinary::REAL energy, lowerBound;

		const int nodeNum = voxels.n_voxels(); // number of nodes
		mrf = new MRFEnergy<TypeBinary>(TypeBinary::GlobalSize());
		nodes = new MRFEnergy<TypeBinary>::NodeId[nodeNum];

		// construct energy
		for (unsigned int voxel_index = 0; voxel_index < voxels.n_voxels(); ++voxel_index)
		{
			nodes[voxel_index] = mrf->AddNode(TypeBinary::LocalSize(), TypeBinary::NodeData(
				voxel_visibility_values[voxel_index], 1.0 - voxel_visibility_values[voxel_index]));
		}

		for (unsigned int voxel_index = 0; voxel_index < voxels.n_voxels(); ++voxel_index)
		{
			VoxelIndex3D xyz_index = voxels.get_voxel_index(voxel_index);

			for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
			{
				VoxelIndex3D n_xyz_index = xyz_index;
				if (n_xyz_index[axis_index] + 1 == voxels.n_axis_voxels(axis_index))
					continue;

				++n_xyz_index[axis_index];
				unsigned int n_voxel_index = voxels.get_voxel_index(n_xyz_index);
				assert(voxel_index < n_voxel_index);
				assert(n_voxel_index < voxels.n_voxels());
				mrf->AddEdge(nodes[voxel_index], nodes[n_voxel_index], TypeBinary::EdgeData(0, 1.0, 1.0, 0));
			}
		}

		/////////////////////// TRW-S algorithm //////////////////////
		options.m_iterMax = 30; // maximum number of iterations
		mrf->Minimize_TRW_S(options, lowerBound, energy);

		// read solution
		std::vector<Real> smoothed_voxel_visibility_values(voxels.n_voxels(), 0.0);
		for (unsigned int voxel_index = 0; voxel_index < voxels.n_voxels(); ++voxel_index)
		{
			smoothed_voxel_visibility_values[voxel_index] =
				static_cast<Real>(mrf->GetSolution(nodes[voxel_index]));
		}

		delete nodes;
		delete mrf;


		// 4. Voxel filling.
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

			if (smoothed_voxel_visibility_values[voxel_index] > 0.5)
			{
				for (std::list<int>::const_iterator it = symmetry_voxel_point_indices.begin();
					it != symmetry_voxel_point_indices.end(); ++it)
				{
					MeshSamplePoint *sample_point = symmetry_cuboid->get_sample_point(*it);
					assert(sample_point);

					MeshSamplePoint *new_sample_point = _output_cuboid_structure.add_sample_point(
						sample_point->point_, sample_point->normal_);
					output_cuboid->add_sample_point(new_sample_point);
				}
			}
			else
			{
				for (std::list<int>::const_iterator it = database_voxel_point_indices.begin();
					it != database_voxel_point_indices.end(); ++it)
				{
					MeshSamplePoint *sample_point = database_cuboid->get_sample_point(*it);
					assert(sample_point);

					MeshSamplePoint *new_sample_point = _output_cuboid_structure.add_sample_point(
						sample_point->point_, sample_point->normal_);
					output_cuboid->add_sample_point(new_sample_point);
				}
			}
		}

		std::cout << "Done." << std::endl;
	}

	// Add unvisited (unsegmented) sample points in symmetry reconstruction.
	for (SamplePointIndex sample_point_index = 0; sample_point_index < _symmetry_reconstruction.num_sample_points();
		++sample_point_index)
	{
		if (!is_symmetry_point_visited[sample_point_index])
		{
			const MeshSamplePoint *sample_point = _symmetry_reconstruction.sample_points_[sample_point_index];
			assert(sample_point);
			MeshSamplePoint *new_sample_point = _output_cuboid_structure.add_sample_point(
				sample_point->point_, sample_point->normal_);
		}
	}

	delete[] is_symmetry_point_visited;
}