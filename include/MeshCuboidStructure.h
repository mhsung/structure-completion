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

#include <vector>
#include <set>


class MeshCuboidStructure {
public:
	MeshCuboidStructure(const MyMesh *_mesh);
	~MeshCuboidStructure();


	typedef struct
	{
		bool is_multiple_cuboids_;
		std::list<Label> mapped_cuboid_labels_;
	}
	PointCuboidLabelMap;


	void clear();

	void clear_labels();

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

	bool load_sample_points(const char *_filename, bool _verbose = true);

	bool load_sample_point_labels(const char *_filename, bool _verbose = true);

	// FIXME.
	// Test code.
	bool load_cuboids(const char *_filename, bool _verbose = true);

	std::vector<MeshCuboid *> get_all_cuboids()const;

	void make_mesh_vertices_as_sample_points();

	// Get sample point labels from the label confidence values.
	std::vector<LabelIndex> get_sample_point_label_indices();

	void compute_label_cuboids();

	// Apple mesh face labels to both sample points and parts,
	// but not modify existing parts except their labels.
	void apply_mesh_face_labels_to_cuboids();

	bool apply_point_cuboid_label_map(
		const std::vector<PointCuboidLabelMap>& _point_cuboid_label_maps,
		const std::vector<Label>& _all_cuboid_labels);

	void apply_test();

	// Apple mesh face labels to sample points,
	// and re-create parts based on sample point labels.
	void get_mesh_face_label_cuboids();

	// Find the largest part for each label.
	void find_the_largest_label_cuboids();

	// Apply labels in new indices.
	bool set_new_label_indices(const std::vector<Label>& _labels);

	void split_label_cuboids();

	void print_label_cuboids(const LabelIndex _label_index)const;

	void remove_occluded_sample_points(const std::set<FaceIndex>& _visible_face_indices);


private:
	inline Label get_new_label()const;

	// Apple mesh face labels to sample points
	void apply_mesh_face_labels_to_sample_points();


public:
	const MyMesh *mesh_;

	std::vector<MeshSamplePoint *> sample_points_;
	std::vector<Label> labels_;
	std::vector< std::vector<MeshCuboid *> > label_cuboids_;

	LabelIndex query_label_index_;

	MyMesh::Normal translation_;
	Real scale_;

	void translate(const MyMesh::Normal _translate);
	void scale(const Real _scale);
	void reset_transformation();
};

#endif	// _MESH_CUBOID_STRUCTURE_H_