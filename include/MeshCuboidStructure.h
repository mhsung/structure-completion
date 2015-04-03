/**
* file:	MeshCuboidStructure.h
* description:	Define a mesh proxy structure.
*
* author:	Minhyuk Sung
* date:	October 2014
*/

#ifndef _MESH_CUBOID_STRUCTURE_H_
#define _MESH_CUBOID_STRUCTURE_H_

#include "MyMesh.h"
#include "MeshCuboid.h"
#include "MeshCuboidSymmetryGroup.h"

#include <vector>
#include <set>


class MeshCuboidStructure {
public:
	MeshCuboidStructure(const MyMesh *_mesh);
	MeshCuboidStructure(const MeshCuboidStructure& _other);	// Copy constructor.
	~MeshCuboidStructure();

	MeshCuboidStructure& operator=(const MeshCuboidStructure& _other);	// Copy assignment.
	void deep_copy(const MeshCuboidStructure& _other);

	void clear();
	void clear_sample_points();
	void clear_cuboids();
	void clear_labels();

	bool load_cuboids(const std::string _filename, bool _verbose = true);
	bool save_cuboids(const std::string _filename, bool _verbose = true) const;

	void apply_mesh_transformation();

	inline unsigned int num_sample_points()const {
		return static_cast<unsigned int>(sample_points_.size());
	}

	inline unsigned int num_labels()const {
		return static_cast<unsigned int>(labels_.size());
	}

	Label get_label(const LabelIndex _label_index)const;

	bool exist_label(const Label _label, LabelIndex* _label_index = NULL)const;

	LabelIndex get_label_index(const Label _label)const;

	LabelIndex get_label_index(const std::string _label_name)const;

	bool load_labels(const char *_filename, bool _verbose = true);

	bool load_label_symmetries(const char *_filename, bool _verbose = true);

	bool load_symmetry_groups(const char *_filename, bool _verbose = true);

	bool load_sample_points(const char *_filename, bool _verbose = true);
	bool save_sample_points(const char *_filename, bool _verbose = true);

	bool load_sample_point_labels(const char *_filename, bool _verbose = true);

	std::vector<MeshCuboid *> get_all_cuboids()const;

	// Get sample point labels from the label confidence values.
	std::vector<LabelIndex> get_sample_point_label_indices();

	MeshSamplePoint *add_sample_point(const MyMesh::Point& _pos);

	void add_sample_points_from_mesh_vertices();

	void compute_label_cuboids();

	// Apple mesh face labels to both sample points and parts,
	// but not modify existing parts except their labels.
	void apply_mesh_face_labels_to_cuboids();

	// Apple mesh face labels to sample points,
	// and re-create parts based on sample point labels.
	void get_mesh_face_label_cuboids(bool _add_sample_points_from_mesh_vertices = false);

	// Find the largest part for each label.
	void find_the_largest_label_cuboids();

	void split_label_cuboids();

	void print_label_cuboids(const LabelIndex _label_index)const;

	void remove_occluded_sample_points(const std::set<FaceIndex>& _visible_face_indices);

	void remove_symmetric_cuboids();

	// NOTICE:
	// Symmetry group functions are replaced with 'MeshCuboidSymmetryGroup' class.
	//bool is_label_group(LabelIndex _label_index);
	//void clear_symmetric_group_labels();
	//void add_symmetric_group_labels();
	//void create_symmetric_group_cuboids();
	//int find_parent_label_index(const LabelIndex _label_index_1, const LabelIndex _label_index_2);

	void compute_symmetry_groups();

	void copy_sample_points_to_symmetric_position();

	void copy_sample_points_to_symmetric_position(
		const MeshCuboidSymmetryGroup* _symmetry_group);

	void copy_sample_points_to_symmetric_position(
		const MeshCuboidSymmetryGroup* _symmetry_group,
		const MeshCuboid *_cuboid_1, MeshCuboid *_cuboid_2);


	// TEST.
	bool test_load_cuboids(const char *_filename, bool _verbose = true);

	// Apply labels in new indices.
	//bool set_new_label_indices(const std::vector<Label>& _labels);

	//typedef struct
	//{
	//	bool is_multiple_cuboids_;
	//	std::list<Label> mapped_cuboid_labels_;
	//}
	//PointCuboidLabelMap;

	//bool apply_point_cuboid_label_map(
	//	const std::vector<PointCuboidLabelMap>& _point_cuboid_label_maps,
	//	const std::vector<Label>& _all_cuboid_labels);


private:
	inline Label get_new_label()const;

	// Apple mesh face labels to sample points
	void apply_mesh_face_labels_to_sample_points();


public:
	const MyMesh *mesh_;

	std::vector<MeshSamplePoint *> sample_points_;
	std::vector<Label> labels_;
	std::vector<std::string> label_names_;
	std::vector< std::list<LabelIndex> > label_symmetries_;
	std::vector< std::vector<MeshCuboid *> > label_cuboids_;
	//std::vector< std::list<LabelIndex> > label_children_;

	std::vector< MeshCuboidSymmetryGroupInfo > symmetry_group_info_;
	std::vector< MeshCuboidSymmetryGroup* > symmetry_groups_;

	MyMesh::Normal translation_;
	Real scale_;
	LabelIndex query_label_index_;

	void translate(const MyMesh::Normal _translate);
	void scale(const Real _scale);
	void reset_transformation();
};

#endif	// _MESH_CUBOID_STRUCTURE_H_