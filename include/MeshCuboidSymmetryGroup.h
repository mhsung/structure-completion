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
	MeshCuboidSymmetryGroupInfo(MeshCuboidSymmetryGroupType _symmetry_type,
		unsigned int _aligned_global_axis_index);
	MeshCuboidSymmetryGroupInfo(const MeshCuboidSymmetryGroupInfo& _other);

	MeshCuboidSymmetryGroupType symmetry_type_;
	unsigned int aligned_global_axis_index_;
	std::vector<LabelIndex> single_label_indices_;
	std::vector< std::pair<LabelIndex, LabelIndex> > pair_label_indices_;
};


class MeshCuboidSymmetryGroup
{
public:
	MeshCuboidSymmetryGroup();
	MeshCuboidSymmetryGroup(const MeshCuboidSymmetryGroupInfo &_info);
	MeshCuboidSymmetryGroup(const MeshCuboidSymmetryGroup &_other);
	~MeshCuboidSymmetryGroup();

	struct WeightedPointPair
	{
		WeightedPointPair(const MyMesh::Point _p1, const MyMesh::Point _p2,
			const Real _weight = 1.0, Real _angle = 0.0)
			: p1_(_p1), p2_(_p2), weight_(_weight), angle_(_angle) {}
		MyMesh::Point p1_;
		MyMesh::Point p2_;
		Real weight_;
		Real angle_;
	};


	// Virtual functions.
	virtual bool compute_symmetry_axis(const std::vector<MeshCuboid *>& _cuboids) = 0;

	virtual MyMesh::Point get_symmetric_point(const MyMesh::Point& _point, unsigned int _symmetry_order) const = 0;
	virtual MyMesh::Normal get_symmetric_normal(const MyMesh::Normal& _normal, unsigned int _symmetry_order) const = 0;

	virtual MeshCuboidSymmetryGroupType get_symmetry_type() const = 0;
	virtual unsigned int num_symmetry_orders() const;

	Real get_rotation_angle() const;
	unsigned int get_aligned_global_axis_index() const;
	void get_single_cuboid_indices(const std::vector<MeshCuboid *>& _cuboids,
		std::vector<unsigned int> &_single_cuboid_indices) const;
	void get_pair_cuboid_indices(const std::vector<MeshCuboid *>& _cuboids,
		std::vector< std::pair<unsigned int, unsigned int> > &_pair_cuboid_indices) const;

	void get_symmetric_sample_point_pairs(
		const std::vector<MeshCuboid *> &_cuboids,
		const std::vector<ANNpointArray> &_cuboid_ann_points,
		const std::vector<ANNkd_tree *> &_cuboid_ann_kd_tree,
		const Real _squared_neighbor_distance,
		std::list<WeightedPointPair> &_sample_point_pairs) const;

	void get_symmetric_sample_point_pairs(
		const MeshCuboid *_cuboid_1,
		const ANNpointArray &_cuboid_ann_points_2,
		ANNkd_tree *_cuboid_ann_kd_tree_2,
		const Real _squared_neighbor_distance,
		std::list<WeightedPointPair> &_sample_point_pairs) const;

	void get_symmetric_sample_point_pairs(
		const std::vector<MyMesh::Point> &_cuboid_1_sample_points,
		const ANNpointArray &_cuboid_ann_points_2,
		ANNkd_tree *_cuboid_ann_kd_tree_2,
		const Real _squared_neighbor_distance,
		std::list<WeightedPointPair> &_sample_point_pairs) const;

	MeshCuboidSymmetryGroupInfo get_symmetry_group_info()const { return info_; }

protected:
	const MeshCuboidSymmetryGroupInfo info_;
	unsigned int num_symmetry_orders_;
};


class MeshCuboidReflectionSymmetryGroup : public MeshCuboidSymmetryGroup
{
public:
	static MeshCuboidReflectionSymmetryGroup* constructor(
		const MeshCuboidSymmetryGroupInfo &_info,
		const std::vector<MeshCuboid *>& _cuboids);
	virtual ~MeshCuboidReflectionSymmetryGroup();

	MeshCuboidReflectionSymmetryGroup(const MyMesh::Normal _n, const double _t);
	MeshCuboidReflectionSymmetryGroup(const MeshCuboidSymmetryGroupInfo &_info);
	MeshCuboidReflectionSymmetryGroup(const MeshCuboidReflectionSymmetryGroup &_other);


	// Virtual functions.
	virtual bool compute_symmetry_axis(const std::vector<MeshCuboid *>& _cuboids);

	virtual MyMesh::Point get_symmetric_point(const MyMesh::Point& _point, unsigned int _symmetry_order) const;
	virtual MyMesh::Normal get_symmetric_normal(const MyMesh::Normal& _normal, unsigned int _symmetry_order) const;
	virtual MyMesh::Point get_symmetric_point(const MyMesh::Point& _point) const;
	virtual MyMesh::Normal get_symmetric_normal(const MyMesh::Normal& _normal) const;

	virtual MeshCuboidSymmetryGroupType get_symmetry_type() const;
	virtual unsigned int num_symmetry_orders() const;
	virtual void set_num_symmetry_order(unsigned int _symmetry_order);

	void get_reflection_plane(MyMesh::Normal &_n, double &_t) const;
	void set_reflection_plane(const MyMesh::Normal &_n, const double &_t);
	void get_reflection_plane_corners(const MyMesh::Point &_point, const Real _size,
		std::array<MyMesh::Point, 4>& _corners) const;


	// static functions.
	static unsigned int num_axis_parameters();

	static void add_symmety_cuboid_corner_points(
		const MeshCuboid *cuboid_1, const MeshCuboid *cuboid_2,
		const std::vector<MeshCuboid *>& _cuboids,
		const unsigned int _reflection_axis_index,
		std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs);

	static void compute_reflection_plane(
		const std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs,
		MyMesh::Normal &_n, double &_t);

private:
	// For any point p on the reflection plane, dot(n, p) = t.
	MyMesh::Normal n_;
	double t_;
};


class MeshCuboidRotationSymmetryGroup : public MeshCuboidSymmetryGroup
{
public:
	static MeshCuboidRotationSymmetryGroup* constructor(
		const MeshCuboidSymmetryGroupInfo &_info,
		const std::vector<MeshCuboid *>& _cuboids);
	virtual ~MeshCuboidRotationSymmetryGroup();
	
	MeshCuboidRotationSymmetryGroup(const MeshCuboidSymmetryGroupInfo &_info);
	MeshCuboidRotationSymmetryGroup(const MeshCuboidRotationSymmetryGroup &_other);


	// Virtual functions.
	virtual bool compute_symmetry_axis(const std::vector<MeshCuboid *>& _cuboids);
	
	virtual MyMesh::Point get_symmetric_point(const MyMesh::Point& _point, unsigned int _symmetry_order) const;
	virtual MyMesh::Normal get_symmetric_normal(const MyMesh::Normal& _normal, unsigned int _symmetry_order) const;

	virtual MeshCuboidSymmetryGroupType get_symmetry_type() const;


	bool compute_rotation_angle(const std::vector<MeshCuboid *> &_cuboids);

	void get_rotation_axis(MyMesh::Normal &_n, MyMesh::Point &_t) const;
	void set_rotation_axis(const MyMesh::Normal &_n, const MyMesh::Point &_t);
	void get_rotation_axis_corners(const MyMesh::Point &_point, const Real _size,
		std::array<MyMesh::Point, 2>& _corners) const;


	// static functions.
	static unsigned int num_axis_parameters();

	static void add_symmety_cuboid_corner_points(
		const MeshCuboid *cuboid_1, const MeshCuboid *cuboid_2,
		const std::vector<MeshCuboid *>& _cuboids,
		const unsigned int _reflection_axis_index,
		std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs);

	static void compute_rotation_plane(
		const std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs,
		MyMesh::Normal &_n, double &_t);

private:
	// For any point p on the rotation axis, p = nx + t (x is scalar parameter).
	MyMesh::Normal n_;
	MyMesh::Point t_;
};

#endif	// _MESH_CUBOID_SYMMETRY_GROUP_H_