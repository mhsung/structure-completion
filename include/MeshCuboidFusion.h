#ifndef _MESH_CUBOID_FUSION_H_
#define _MESH_CUBOID_FUSION_H_

#include "MeshCuboidStructure.h"
#include "ICP.h"

#include <Eigen/Core>


typedef OpenMesh::Vec3i MeshCuboidVoxelIndex3D;
class MeshCuboidVoxelGrid
{
public:
	MeshCuboidVoxelGrid(MyMesh::Point _min, MyMesh::Point _max, Real _unit_size);
	~MeshCuboidVoxelGrid();

	int n_voxels() const;
	int n_axis_voxels(const int _axis_index) const;
	int get_voxel_index(MeshCuboidVoxelIndex3D _xyz_index) const;
	MeshCuboidVoxelIndex3D get_voxel_index(const int _voxel_index) const;
	int get_voxel_index(const MyMesh::Point _point) const;
	MyMesh::Point get_center(const int _voxel_index) const;
	void get_centers(std::vector<MyMesh::Point> &_centers) const;
	void get_point_correspondences(
		const std::vector<MyMesh::Point> &_points,
		std::vector<int> &_points_to_voxels,
		std::vector< std::list<int> > &_voxels_to_points) const;
	void get_voxel_occupancies(
		const std::vector<MyMesh::Point> &_points,
		Eigen::VectorXd &_voxel_occupancies) const;
	void get_distance_map(
		ANNpointArray &_ann_points, ANNkd_tree *_ann_kd_tree,
		Eigen::VectorXd &_voxel_to_point_distances) const;

private:
	MyMesh::Point min_;
	MyMesh::Point max_;
	MeshCuboidVoxelIndex3D n_voxels_;
};

void run_part_ICP(
	MeshCuboidStructure &_input,
	const MeshCuboidStructure &_ground_truth);

void create_voxel_grid(
	const MeshCuboid *_symmetry_cuboid,
	const MeshCuboid *_database_cuboid,
	MyMesh::Point &_local_coord_bbox_min,
	MyMesh::Point &_local_coord_bbox_max);

void get_smoothed_voxel_visibility(
	const MeshCuboidVoxelGrid &_local_coord_voxels,
	const MeshCuboid *_ground_truth_cuboid,
	const double *_occlusion_modelview_matrix,
	const MeshCuboidStructure &_original_cuboid_structure,
	const Real &_occlusion_radius,
	const Real &_smoothing_parameter,
	std::vector<Real> &_voxel_visibility);

void merge_symmetric_cuboids_visibility(
	const MeshCuboidSymmetryGroup *_symmetry_group,
	const MeshCuboidVoxelGrid &_voxels_1,
	const MeshCuboidVoxelGrid &_voxels_2,
	std::vector<Real> &_voxel_visibility_1,
	std::vector<Real> &_voxel_visibility_2);

void fill_voxels_using_visibility(
	const MeshCuboidVoxelGrid &_voxels,
	const std::vector<Real> &_voxel_visibility_values,
	const MeshCuboid *symmetry_cuboid,
	const MeshCuboid *database_cuboid,
	MeshCuboidStructure &_output_cuboid_structure,
	MeshCuboid *output_cuboid);

void reconstruct_fusion_simple(
	const MeshCuboidStructure &_symmetry_reconstruction,
	const MeshCuboidStructure &_database_reconstruction,
	MeshCuboidStructure &_output_cuboid_structure);

void reconstruct_fusion(const char *_mesh_filepath,
	const double *_snapshot_modelview_matrix,
	const double *_occlusion_modelview_matrix,
	const MeshCuboidStructure &_original_cuboid_structure,
	const MeshCuboidStructure &_symmetry_reconstruction,
	const MeshCuboidStructure &_database_reconstruction,
	MeshCuboidStructure &_output_cuboid_structure,
	bool _add_outliers = true);

#endif	// _MESH_CUBOID_FUSION_H_