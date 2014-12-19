/**
* file:	MeshCuboidStats.h
* description:	Define a mesh proxy part relations.
*
* author:	Minhyuk Sung
* date:	October 2014
*/

#ifndef _MESH_CUBOID_STAT_H_
#define _MESH_CUBOID_STAT_H_

//#define COMPUTE_LINEAR_MAPS


#include "MyMesh.h"
#include "MeshCuboid.h"
#include "MeshCuboidStructure.h"
#include "MeshCuboidRelation.h"

#include <list>
#include <map>
#include <set>
#include <string>

#include <Eigen/Core>


class MeshCuboidStats {
public:
	MeshCuboidStats();
	~MeshCuboidStats();

	friend class MeshCuboidManualFeatures;

	typedef struct  
	{
		unsigned int count_;
		Real mean_;
		Real stdev_;
	}
	StatParams;	

	void initialize();

	unsigned int num_stats()const;

	bool load_stats(const char* _filename, bool _verbose = true);
	bool save_stats(const char* _filename, bool _verbose = true);
	void compute_stats(const std::vector<std::string>& _keys,
		const std::vector<std::vector<Real>>& _values);

private:
	std::vector<std::string> keys_;
	std::vector<StatParams> stats_;
};

class MeshCuboidManualFeatures {
public:
	MeshCuboidManualFeatures(const std::string _object_name);
	virtual ~MeshCuboidManualFeatures();

	const MyMesh::Normal k_up_direction = MeshCuboidAttributes::k_up_direction;

	void clear();

	static bool get_keys_and_values(const std::list<MeshCuboidManualFeatures *>& _stats,
		std::vector<std::string>& _keys,
		std::vector<std::vector<Real>>& _values);

	static void get_object_names(const std::list<MeshCuboidManualFeatures *>& _stats,
		std::vector<std::string>& _object_names);

	static bool save_keys_and_values(const std::list<MeshCuboidManualFeatures *>& _stats,
		const char* _filename);

	static bool save_stats(
		const std::list<MeshCuboidManualFeatures *>& _stats, const char* _filename);

	Real compute_log_probability(const MeshCuboidStats& _stats)const;


protected:
	std::string object_name_;
	std::vector<std::string> keys_;
	std::vector<Real> values_;

#ifdef COMPUTE_LINEAR_MAPS
	// x = Attributes.
	// y = Features.
	// A = Linear map.
	// y = A*x.
	std::vector<std::vector<Eigen::VectorXd>> linear_map_;
#endif
};


class MeshSingleCuboidManualFeatures : public MeshCuboidManualFeatures {
public:
	MeshSingleCuboidManualFeatures(const std::string _object_name);
	virtual ~MeshSingleCuboidManualFeatures();

	virtual void compute_values(const MeshCuboid *_cuboid);
	//virtual void compute_values(MeshCuboidStructure *_cuboid_structure);

private:
	// NOTE:
	// The prefix should be the name of probability distribution.
	//inline std::string get_key_bbox_axes(
	//	const Label _label, const unsigned int _axis_index, const unsigned int _i)const;
	//inline std::string get_key_bbox_center(
	//	const Label _label, const unsigned int _axis_index)const;
	inline std::string get_key_bbox_size(
		const unsigned int _axis_index)const;
	//inline std::string get_key_bbox_height(
	//	const unsigned int _i)const;
	inline std::string get_key_bbox_corner_height(
		const unsigned int _corner_index)const;
	inline std::string get_key_bbox_axes_up_angle(
		const unsigned int _axis_index)const;

	//Real get_value_bbox_axes(
	//	const MeshCuboid *_cuboid, const unsigned int _axis_index, const unsigned int _i)const;
	//Real get_value_bbox_center(
	//	const MeshCuboid *_cuboid, const unsigned int _axis_index)const;
	Real get_value_bbox_size(
		const MeshCuboid *_cuboid, const unsigned int _axis_index)const;
	//Real get_value_bbox_height(
	//	const MeshCuboid *_cuboid, const unsigned int _i)const;
	Real get_value_bbox_corner_height(
		const MeshCuboid *_cuboid, const unsigned int _corner_index)const;
	Real get_value_bbox_axis_up_angle(
		const MeshCuboid *_cuboid, const unsigned int _axis_index)const;
};


class MeshPairCuboidManualFeatures : public MeshCuboidManualFeatures {
public:
	MeshPairCuboidManualFeatures(const std::string _object_name);
	virtual ~MeshPairCuboidManualFeatures();

	typedef enum
	{
		CUBOID_1_2,
		CUBOID_2_1,
	} CuboidPairOrder;

	virtual void compute_values(const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2);
	//virtual void compute_values(MeshCuboidStructure *_cuboid_structure);

private:
	// NOTE:
	// The prefix should be the name of probability distribution.
	inline std::string get_key_bbox_axis_angle(
		const unsigned int _axis_index_1, const unsigned int _axis_index_2)const;
	//inline std::string get_key_bbox_size_ratio(
	//	const unsigned int _axis_index_1, const unsigned int _axis_index_2)const;
	inline std::string get_key_bbox_face_corner_distance(
		const CuboidPairOrder _pair_order,
		const unsigned int _face_index_1, const unsigned int _corner_index_2)const;
	inline std::string get_key_bbox_axis_center_distance(
		const CuboidPairOrder _pair_order,
		const unsigned int _axis_index_1)const;

	Real get_value_bbox_axis_angle(
		const MeshCuboid* _cuboid_1, const MeshCuboid* _cuboid_2,
		const unsigned int _axis_index_1, const unsigned int _axis_index_2)const;
	//Real get_value_bbox_size_ratio(
	//	const MeshCuboid* _cuboid_1, const MeshCuboid* _cuboid_2,
	//	const unsigned int _axis_index_1, const unsigned int _axis_index_2)const;
	Real get_value_bbox_face_corner_distance(
		const MeshCuboid* _cuboid_1, const MeshCuboid* _cuboid_2,
		const unsigned int _face_index_1, const unsigned int _corner_index_2)const;
	Real get_value_bbox_axis_center_distance(
		const MeshCuboid* _cuboid_1, const MeshCuboid* _cuboid_2,
		const unsigned int _axis_index_1)const;

	void compute_values_bbox_axis_angle(
		const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2);
	//void compute_values_bbox_size_ratio(
	//	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2);
	void compute_values_bbox_face_corner_distance(
		const CuboidPairOrder _pair_order,
		const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2);
	void compute_values_bbox_axis_center_distance(
		const CuboidPairOrder _pair_order,
		const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2);
};

#endif	// _MESH_CUBOID_STAT_H_