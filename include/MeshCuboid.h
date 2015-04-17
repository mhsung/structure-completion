#ifndef _MESH_CUBOID_H_
#define _MESH_CUBOID_H_

//#define AXIS_ALIGNED_INITIAL_CUBOID
#define MIN_CUBOID_SIZE						1.0E-12
#define CUBOID_SURFACE_SAMPLING_RANDOM_SEED	20130923

#include "ICP.h"
#include "MyMesh.h"

#include <array>
#include <map>
#include <ANN/ANN.h>
#include <Eigen/Core>


// NOTE:
// Each label has a label index.
// In 'MeshCuboidStructure' class, the label can be obtained with a label index
// with the function 'get_label()'. 
typedef unsigned int LabelIndex;
typedef unsigned int SamplePointIndex;

class MeshSamplePoint
{
public:
	MeshSamplePoint(
		SamplePointIndex _sample_point_index,
		FaceIndex _corr_fid,
		MyMesh::Point _bary_coord,
		MyMesh::Point _point,
		MyMesh::Normal _normal)
		: sample_point_index_(_sample_point_index)
		, corr_fid_(_corr_fid)
		, bary_coord_(_bary_coord)
		, point_(_point)
		, normal_(_normal)
		, error_(0)
	{}

	MeshSamplePoint(
		const MeshSamplePoint& _other)
		: sample_point_index_(_other.sample_point_index_)
		, corr_fid_(_other.corr_fid_)
		, bary_coord_(_other.bary_coord_)
		, point_(_other.point_)
		, normal_(_other.normal_)
		, label_index_confidence_(_other.label_index_confidence_)
		, error_(_other.error_)
	{}

	SamplePointIndex sample_point_index_;
	FaceIndex corr_fid_;
	MyMesh::Point bary_coord_;
	MyMesh::Point point_;
	MyMesh::Normal normal_;
	std::vector<Real> label_index_confidence_;
	Real error_;
};

class MeshCuboidSurfacePoint
{
public:
	MeshCuboidSurfacePoint(
		const MyMesh::Point _point,
		const MyMesh::Normal _normal,
		unsigned int _cuboid_face_index,
		const std::array<Real, 8>& _corner_weights)
		: point_(_point)
		, normal_(_normal)
		, cuboid_face_index_(_cuboid_face_index)
		, corner_weights_(_corner_weights)
		, visibility_(0) 
	{}

	MeshCuboidSurfacePoint(
		const MeshCuboidSurfacePoint& _other)
		: point_(_other.point_)
		, normal_(_other.normal_)
		, cuboid_face_index_(_other.cuboid_face_index_)
		, corner_weights_(_other.corner_weights_)
		, visibility_(_other.visibility_)
	{}

	MyMesh::Point point_;
	MyMesh::Normal normal_;
	unsigned int cuboid_face_index_;
	std::array<Real, 8> corner_weights_;
	Real visibility_;
};

class MeshCuboid
{
public:
	MeshCuboid(const LabelIndex _label_index);
	MeshCuboid(const MeshCuboid& _other);	// Copy constructor.
	virtual ~MeshCuboid();

	friend class MeshCuboidStructure;

	MeshCuboid& operator=(const MeshCuboid& _other);	// Copy assignment.
	void deep_copy(const MeshCuboid& _other);

	static const unsigned int k_num_corners = 8;
	static const unsigned int k_num_faces = 6;
	static const unsigned int k_num_edges = 12;
	static const unsigned int k_num_face_corners = 4;
	static const unsigned int k_face_corner_indices[k_num_faces][k_num_face_corners];
	static const unsigned int k_face_edges[k_num_edges][2];

	typedef enum {
		POSITIVE_X_AXIS,
		NEGATIVE_X_AXIS,
		POSITIVE_Y_AXIS,
		NEGATIVE_Y_AXIS,
		POSITIVE_Z_AXIS,
		NEGATIVE_Z_AXIS
	} AXIS_DIRECTION;

	// List of all axis configurations.
	// first - x-axis, second - y-axis (z-axis is determined by x- and y-axis).
	static const std::vector < std::pair<AXIS_DIRECTION, AXIS_DIRECTION> >
		k_all_axis_configuration;

	static unsigned int num_axis_configurations() {
		return static_cast<unsigned int>(k_all_axis_configuration.size());
	};

	static std::array<MyMesh::Normal, 3>
		get_transformed_axes(const unsigned int _axis_configuration_index,
		const std::array<MyMesh::Normal, 3> &_axes);

	void set_axis_configuration(const unsigned int _axis_configuration_index);


	void clear_sample_points();

	void add_sample_point(MeshSamplePoint *_point);

	void add_sample_points(const std::vector<MeshSamplePoint *> _points);

	void remove_sample_points(const bool *is_sample_point_removed);

	bool compute_bbox();


	unsigned int num_sample_points()const {
		return static_cast<unsigned int>(sample_points_.size());
	}
	unsigned int num_cuboid_surface_points()const {
		return static_cast<unsigned int>(cuboid_surface_points_.size());
	}

	bool is_point_inside_cuboid(const MyMesh::Point& _point)const;

	bool is_group_cuboid()const { return is_group_cuboid_;  }

	void set_group_cuboid(bool _is_group_cuboid) { is_group_cuboid_ = _is_group_cuboid; }

	LabelIndex get_label_index()const { return label_index_; }

	const std::vector<MeshSamplePoint *> &get_sample_points()const;

	const MeshSamplePoint *get_sample_point(const unsigned int _point_index)const;

	const std::vector<MeshCuboidSurfacePoint *> &get_cuboid_surface_points()const;

	MeshCuboidSurfacePoint *get_cuboid_surface_point(const unsigned int _point_index)const;

	void set_label_index(LabelIndex _label_index) { label_index_ = _label_index; }

	MyMesh::Point get_bbox_min()const;

	MyMesh::Point get_bbox_max()const;

	MyMesh::Normal get_bbox_axis(const unsigned int _axis_index)const;

	std::array<MyMesh::Normal, 3> get_bbox_axes()const;

	MyMesh::Point get_bbox_center()const;

	MyMesh::Normal get_bbox_size()const;

	MyMesh::Point get_bbox_corner(const unsigned int _corner_index,
		MyMesh::Normal *_axis_direction = NULL)const;

	std::array<MyMesh::Point, k_num_corners> get_bbox_corners()const;

	std::pair<MyMesh::Normal, MyMesh::Point> get_bbox_faces(const unsigned int _face_index,
		MyMesh::Normal *_axis_direction = NULL)const;

	std::vector< std::pair<MyMesh::Normal, MyMesh::Point> > get_bbox_faces()const;

	Real get_bbox_volume()const;

	Real get_bbox_diag_length()const;

	Real get_bbox_face_area(const unsigned int _face_index)const;

	const std::vector<int>& get_sample_to_cuboid_surface_correspondences() const;

	int get_sample_to_cuboid_surface_correspondences(const unsigned int _point_index) const;

	const std::vector<int>& get_cuboid_surface_to_sample_correspondence() const;

	int get_cuboid_surface_to_sample_correspondence(const unsigned int _point_index) const;


	void set_bbox_center(const MyMesh::Point &_bbox_center);

	void set_bbox_size(const MyMesh::Normal &_bbox_size, bool _update_corners = true);

	void set_bbox_axes(const std::array<MyMesh::Normal, 3> &_bbox_axes, bool _update_corners = true);

	void set_bbox_corners(const std::array<MyMesh::Point, k_num_corners> &_bbox_corners);

	void translate(const Eigen::Vector3d _translation_vec);

	void rotate(const Eigen::Matrix3d _rotation_mat, bool _update_center_size = true);

	void flip_axis(const unsigned int _axis_index);
	
	void update_corner_points();

	void update_center_size_corner_points();

	void update_axes_center_size_corner_points();

	// Update the label based on the label confidence values of sample points.
	void update_label_using_sample_points();

	static MeshCuboid *merge_cuboids(const LabelIndex _label_index,
		const std::vector<MeshCuboid *> _cuboids);

	std::vector<MeshCuboid *> split_cuboid(const Real _object_diameter);

	void clear_cuboid_surface_points();

	void create_random_points_on_cuboid_surface(
		const unsigned int _num_cuboid_surface_points);

	void create_grid_points_on_cuboid_surface(
		const unsigned int _num_cuboid_surface_points);

	void compute_cuboid_surface_point_visibility(
		const Real _modelview_matrix[16],
		const Real _radius,
		const std::vector<MeshSamplePoint *>& _all_sample_points);

	Real get_cuboid_overvall_visibility();

	void update_point_correspondences();

	// NOTE:
	// Do not use corner points, but use only center, size, and axes.
	// The distance becomes less than zero when the point is inside the cuboid.
	void points_to_cuboid_distances(const Eigen::MatrixXd& _points,
		Eigen::VectorXd &_distances);

	static Real distance_between_cuboids(
		const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2);

	void print_cuboid()const;

	void draw_cuboid()const;


protected:
	bool is_group_cuboid_;

	LabelIndex label_index_;
	std::vector<MeshSamplePoint *> sample_points_;
	std::vector<MeshCuboidSurfacePoint *> cuboid_surface_points_;

	std::vector<int> sample_to_cuboid_surface_correspondence_;
	std::vector<int> cuboid_surface_to_sample_corresopndence_;

	std::array<MyMesh::Normal, 3> bbox_axes_;
	MyMesh::Point bbox_center_;
	MyMesh::Normal bbox_size_;
	std::array<MyMesh::Point, k_num_corners> bbox_corners_;

	void compute_axis_aligned_bbox();

	void compute_oriented_bbox();

	void create_sub_cuboids(const Real _object_diameter, ANNkd_tree* _kd_tree, std::vector<MeshCuboid *> &_sub_cuboids);

	void remove_small_sub_cuboids(std::vector<MeshCuboid *> &_sub_cuboids);

	//void align_sub_cuboids(const Real _object_diameter, std::vector<MeshCuboid *> &_sub_cuboids);
	//void split_cuboid_recursive(ANNkd_tree* _kd_tree, std::vector<MeshCuboid *> &_sub_cuboids);
	//void split_cuboid_once(ANNkd_tree* _kd_tree, std::vector<MeshCuboid *> &_sub_cuboids);
};

#endif // _MESH_CUBOID_H_