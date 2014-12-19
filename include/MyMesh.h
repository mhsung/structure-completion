/**
* file:	MyMesh.h
* description:	Define a mesh.
*
* author:	Minhyuk Sung
* date:	August 2009
*		(Revised) October 2014
*/

#ifndef _MY_MESH_H_
#define _MY_MESH_H_

#define CURVATURE_SMOOTHING			8
#define LABEL_COLORING_RANDOM_SEED	20130923

#define SWAP(a, b, t)	((t) = (a), (a) = (b), (b) = (t))
#define NEXT(i) ((i)<2 ? (i)+1 : (i)-2)		// Trimesh
#define PREV(i) ((i)>0 ? (i)-1 : (i)+2)

#include <list>
#include <vector>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

typedef double Real;
typedef std::vector<Real> RealArray;

// NOTE:
// Negative index indicates undefined.

typedef int VertexIndex;
typedef std::vector<VertexIndex> VertexIndexArray;

typedef int EdgeIndex;
typedef std::vector<EdgeIndex> EdgeIndexArray;

typedef int FaceIndex;
typedef std::vector<FaceIndex> FaceIndexArray;

typedef int Label;
typedef std::vector<Label> LabelArray;


struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef Real				Scalar;
	typedef OpenMesh::Vec3d		Point;
	typedef OpenMesh::Vec3d		Normal;
	typedef OpenMesh::Vec3uc	Color;

	VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
	FaceAttributes( OpenMesh::Attributes::Color );

	/*
	VertexTraits
	{
	public:
		Scalar sample;

	public:
		VertexT(): sample(0.0f) { }
	};

	MyMesh::VertexData &vd = this->data(v_it.handle());
	vd.sample = value;
	*/
};

//typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  MyMesh;


struct MyMesh : OpenMesh::TriMesh_ArrayKernelT<MyTraits>
{
public:
	MyMesh();
	virtual ~MyMesh();

	typedef enum
	{
		VERTEX_COLOR,
		FACE_COLOR,
	}
	MeshColorMap;

	virtual void clear();
	void initialize();

	void translate(const Normal _translate);
	void scale(const Real _scale);
	void reset_transformation();

	// (vertex area) = (sum of adjacent face areas) / (num of adjacent faces)
	void request_vertex_areas();
	void release_vertex_areas();

	void request_face_areas();
	void release_face_areas();

	void request_curvatures();
	void release_curvatures();

	MyMesh::Point get_bbox_center()const { return bbox_center_; }
	MyMesh::Normal get_bbox_size()const { return bbox_size_; }
	Real get_object_diameter()const { return object_diameter_; }

	MyMesh::Point get_translation()const { return translation_; }
	Real get_scale()const { return scale_; }

	bool load_feature_vertices(const char *_filename, bool _verbose = true);
	bool save_feature_vertices(const char *_filename, bool _verbose = true) const;

	bool load_vertex_color_map(const char *_filename, bool _verbose = true);
	bool save_vertex_color_map(const char *_filename, bool _verbose = true) const;

	bool load_face_color_map(const char *_filename, bool _verbose = true);
	bool save_face_color_map(const char *_filename, bool _verbose = true) const;

	bool load_face_label_simple(const char *_filename, bool _verbose = true);
	bool save_face_label_simple(const char *_filename, bool _verbose = true) const;

	bool load_face_label(const char *_filename, bool _verbose = true);
	bool save_face_label(const char *_filename, bool _verbose = true) const;

	void get_all_label_faces(std::list<std::vector<FaceIndex>> &_all_label_faces) const;
	void set_vertex_label_from_face_label();

	void clear_colors();
	
	void set_vertex_color_map();
	void set_vertex_color_map(const RealArray &_values);

	void set_face_color_map();
	void set_face_color_map(const RealArray &_values);

	void set_face_label_colors(bool _verbose = true);

	void extract_local_min_max_feature_vertices(const RealArray &_values,
		const unsigned int _local_search_depth,
		const bool _extract_min, const bool _extract_max);

	void extract_zero_value_feature_vertices(const RealArray &_values,
		const Real _zero_threshold = 0.01);

	void make_face_normal_consistent();

	void print_vertex_information(unsigned int _vid) const;
	void print_user_defined_verices_information() const;

	inline void get_edge_neighbor_vertex_indices(const unsigned int _eid,
		VertexIndex &_v1, VertexIndex &_v2) const;

	static Color get_label_color(const Label _label);

	inline static void gray_to_rgb_color(const Real _gray, Real &_r, Real &_g, Real &_b);


private:
	bool load_color_map(const char *_filename, bool _verbose = true);
	bool save_color_map(const char *_filename, bool _verbose = true) const;

	void is_local_min_max(const RealArray &_values,
		const VertexIndex _vid, const VertexIndex _n_vid,
		const unsigned int _local_search_depth,
		bool &is_local_min, bool &_is_local_max) const;


public:
	// Areas.
	OpenMesh::FPropHandleT<MyMesh::Scalar> face_area_;
	OpenMesh::VPropHandleT<MyMesh::Scalar> vertex_area_;

	// Curvatures.
	OpenMesh::VPropHandleT<MyMesh::Scalar> principal_curvature_1_, principal_curvature_2_;
	OpenMesh::VPropHandleT<MyMesh::Normal> principal_direction_1_, principal_direction_2_;
	OpenMesh::VPropHandleT<MyMesh::Scalar> mean_curvature_, gaussian_curvature_;

	// Color map values.
	OpenMesh::VPropHandleT<MyMesh::Scalar> vertex_color_map_value_;
	OpenMesh::FPropHandleT<MyMesh::Scalar> face_color_map_value_;

	// Labels.
	OpenMesh::VPropHandleT<Label> vertex_label_;
	OpenMesh::FPropHandleT<Label> face_label_;

	// User-defined vertex indices.
	VertexIndex seed_vertex_index_;
	VertexIndex query_vertex_index_;
	VertexIndexArray feature_vertex_indices_;

	MeshColorMap mesh_coloring_option_;


private:
	// Check whether custom properties are initialized.
	bool prop_face_areas_;
	bool prop_vertex_areas_;
	bool prop_curvatures_;

	// Axis-aligned bounding box.
	MyMesh::Point bbox_center_;
	MyMesh::Normal bbox_size_;
	Real object_diameter_;

	MyMesh::Normal translation_;
	Real scale_;
};

template<typename T> bool print_vector(const std::string &_filename, std::vector<T> &_ar)
{
	std::ofstream file(_filename.c_str());
	if(!file)
	{
		std::cerr << "Can't open loading file: \"" << _filename << "\"" << std::endl;
		return false;
	}
	std::cout << "Saving " << _filename << "..." << std::endl;

	for(std::vector<T>::iterator it = _ar.begin(); it != _ar.end(); ++it)
		file << *it << std::endl;

	file.close();
	std::cout << "Done." << std::endl;
	return true;
}

#endif	// _MY_MESH_H_