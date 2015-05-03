/**
* file:	MeshCuboidSymmetryGroup.h
* description:	Define a partial symmetry group of cuboids.
*
* author:	Minhyuk Sung
* date:	March 2015
*/

#ifndef _MESH_CUBOID_SYMMETRY_GROUP_H_
#define _MESH_CUBOID_SYMMETRY_GROUP_H_

#include "MyMesh.h"
#include "MeshCuboid.h"

#include "ANN/ANN.h"
#include <array>
#include <vector>


typedef enum
{
	ReflectionSymmetryType,
	RotationSymmetryType
}
MeshCuboidSymmetryGroupType;


struct MeshCuboidSymmetryGroupInfo
{
	MeshCuboidSymmetryGroupInfo();
	MeshCuboidSymmetryGroupInfo(unsigned int _reflection_axis_index);
	MeshCuboidSymmetryGroupInfo(const MeshCuboidSymmetryGroupInfo& _other);

	MeshCuboidSymmetryGroupType symmetry_type_;
	unsigned int reflection_axis_index_;
	std::vector<LabelIndex> single_label_indices_;
	std::vector< std::pair<LabelIndex, LabelIndex> > pair_label_indices_;
};


class MeshCuboidSymmetryGroup
{
public:
	// NOTE:
	// Use factory constructor.
	static MeshCuboidSymmetryGroup* constructor(
		const MeshCuboidSymmetryGroupInfo &_info,
		const std::vector<MeshCuboid *>& _cuboids);
	~MeshCuboidSymmetryGroup();

	virtual MeshCuboidSymmetryGroup *copy_constructor() = 0;

protected:
	MeshCuboidSymmetryGroup(const MeshCuboidSymmetryGroupInfo &_info);
	MeshCuboidSymmetryGroup(const MeshCuboidSymmetryGroupInfo &_info,
		const MyMesh::Normal _n, const double _t);

public:
	struct WeightedPointPair
	{
		WeightedPointPair(const Real _weight, const MyMesh::Point _p1, const MyMesh::Point _p2)
			: weight_(_weight), p1_(_p1), p2_(_p2) {}
		Real weight_;
		MyMesh::Point p1_;
		MyMesh::Point p2_;
	};

	virtual bool compute_symmetry_axis(const std::vector<MeshCuboid *>& _cuboids) = 0;
	virtual MyMesh::Point get_symmetric_point(const MyMesh::Point& _point) const = 0;
	virtual MyMesh::Normal get_symmetric_normal(const MyMesh::Normal& _normal) const = 0;

	MeshCuboidSymmetryGroupType get_symmetry_type() { return info_.symmetry_type_; }
	unsigned int get_reflection_axis_index() const { return info_.reflection_axis_index_; }
	void get_single_cuboid_indices(const std::vector<MeshCuboid *>& _cuboids,
		std::vector<unsigned int> &_single_cuboid_indices) const;
	void get_pair_cuboid_indices(const std::vector<MeshCuboid *>& _cuboids,
		std::vector< std::pair<unsigned int, unsigned int> > &_pair_cuboid_indices) const;
	void get_axis_parameters(MyMesh::Normal &_n, double &_t) const;
	void set_axis_parameters(const MyMesh::Normal &_n, const double &_t);

	void get_symmetric_sample_point_pairs(
		const std::vector<MeshCuboid *> &_cuboids,
		const std::vector<ANNpointArray> &_cuboid_ann_points,
		const std::vector<ANNkd_tree *> &_cuboid_ann_kd_tree,
		const Real _neighbor_distance,
		std::list<WeightedPointPair> &_sample_point_pairs) const;

	virtual void get_symmetric_sample_point_pairs(
		const MeshCuboid *_cuboid_1,
		const ANNpointArray &_cuboid_ann_points_2,
		ANNkd_tree *_cuboid_ann_kd_tree_2,
		const Real _neighbor_distance,
		std::list<WeightedPointPair> &_sample_point_pairs) const = 0;

protected:
	const MeshCuboidSymmetryGroupInfo info_;

	// For any point p on the reflection plane, dot(n, p) = t.
	MyMesh::Normal n_;
	double t_;
};


class MeshCuboidReflectionSymmetryGroup : public MeshCuboidSymmetryGroup
{
public:
	MeshCuboidReflectionSymmetryGroup(const MeshCuboidSymmetryGroupInfo &_info);
	MeshCuboidReflectionSymmetryGroup(const MeshCuboidSymmetryGroupInfo &_info,
		const MyMesh::Normal _n, const double _t);
	virtual MeshCuboidSymmetryGroup *copy_constructor();

	virtual bool compute_symmetry_axis(const std::vector<MeshCuboid *>& _cuboids);
	virtual MyMesh::Point get_symmetric_point(const MyMesh::Point& _point) const;
	virtual MyMesh::Normal get_symmetric_normal(const MyMesh::Normal& _normal) const;
	virtual void get_symmetric_sample_point_pairs(
		const MeshCuboid *_cuboid_1,
		const ANNpointArray &_cuboid_ann_points_2,
		ANNkd_tree *_cuboid_ann_kd_tree_2,
		const Real _neighbor_distance,
		std::list<WeightedPointPair> &_sample_point_pairs) const;

	void get_reflection_plane_corners(const MyMesh::Point &_point, const Real _size,
		std::array<MyMesh::Point, 4>& _corners) const;

	// static functions.
	static void add_symmety_cuboid_corner_points(
		const MeshCuboid *cuboid_1, const MeshCuboid *cuboid_2,
		const std::vector<MeshCuboid *>& _cuboids,
		const unsigned int _reflection_axis_index,
		std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs);

	static void compute_reflection_plane(
		const std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs,
		MyMesh::Normal &_n, double &_t);
};


class MeshCuboidRotationSymmetryGroup : public MeshCuboidSymmetryGroup
{
public:
	MeshCuboidRotationSymmetryGroup(const MeshCuboidSymmetryGroupInfo &_info);
	MeshCuboidRotationSymmetryGroup(const MeshCuboidSymmetryGroupInfo &_info,
		const MyMesh::Normal _n, const double _t);
	virtual MeshCuboidSymmetryGroup *copy_constructor();

	virtual bool compute_symmetry_axis(const std::vector<MeshCuboid *>& _cuboids);
	virtual MyMesh::Point get_symmetric_point(const MyMesh::Point& _point) const;
	virtual MyMesh::Normal get_symmetric_normal(const MyMesh::Normal& _normal) const;

	virtual void get_symmetric_sample_point_pairs(
		const MeshCuboid *_cuboid_1,
		const ANNpointArray &_cuboid_ann_points_2,
		ANNkd_tree *_cuboid_ann_kd_tree_2,
		const Real _neighbor_distance,
		std::list<WeightedPointPair> &_sample_point_pairs) const;

	// static functions.
	static void add_symmety_cuboid_corner_points(
		const MeshCuboid *cuboid_1, const MeshCuboid *cuboid_2,
		const std::vector<MeshCuboid *>& _cuboids,
		const unsigned int _reflection_axis_index,
		std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs);

	static void compute_rotation_plane(
		const std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs,
		MyMesh::Normal &_n, double &_t);
};

#endif	// _MESH_CUBOID_SYMMETRY_GROUP_H_