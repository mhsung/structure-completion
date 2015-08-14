#include "MeshCuboid.h"

#include "MeshCuboidParameters.h"
#include "Utilities.h"
#include "simplerandom.h"

#include <bitset>
#include <deque>
#include <functional>
#include <limits>
#include <random>
#include <queue>
#include <time.h>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <Eigen/Geometry>

#if __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

// Corners.
// [0]: - - -
// [1]: - - +
// [2]: - + -
// [3]: - + +
// [4]: + - -
// [5]: + - +
// [6]: + + -
// [7]: + + +

const unsigned int MeshCuboid::k_face_corner_indices[k_num_faces][k_num_face_corners] = {
	{ 1u, 3u, 7u, 5u }, // POSITIVE_X_AXIS.
	{ 0u, 4u, 6u, 2u }, // NEGATIVE_X_AXIS.
	{ 2u, 6u, 7u, 3u }, // POSITIVE_Y_AXIS.
	{ 0u, 1u, 5u, 4u }, // NEGATIVE_Y_AXIS.
	{ 4u, 5u, 7u, 6u }, // POSITIVE_Z_AXIS.
	{ 0u, 2u, 3u, 1u }  // NEGATIVE_Z_AXIS.
};

const unsigned int MeshCuboid::k_face_edges[k_num_edges][2] = {
	{ 0u, 4u }, { 1u, 5u }, { 2u, 6u }, { 3u, 7u },
	{ 0u, 2u }, { 1u, 3u }, { 4u, 6u }, { 5u, 7u },
	{ 0u, 1u }, { 2u, 3u }, { 4u, 5u }, { 6u, 7u }
};

const std::vector < std::pair<MeshCuboid::AXIS_DIRECTION, MeshCuboid::AXIS_DIRECTION> >
MeshCuboid::k_all_axis_configuration = {
	std::make_pair(POSITIVE_X_AXIS, POSITIVE_Y_AXIS),	// +X, +Y, +Z.
	std::make_pair(POSITIVE_X_AXIS, NEGATIVE_Y_AXIS),	// +X, -Y, -Z.
	std::make_pair(POSITIVE_X_AXIS, POSITIVE_Z_AXIS),	// +X, +Z, -Y.
	std::make_pair(POSITIVE_X_AXIS, NEGATIVE_Z_AXIS),	// +X, -Z, +Y.

	std::make_pair(NEGATIVE_X_AXIS, POSITIVE_Y_AXIS),	// -X, +Y, -Z.
	std::make_pair(NEGATIVE_X_AXIS, NEGATIVE_Y_AXIS),	// -X, -Y, +Z.
	std::make_pair(NEGATIVE_X_AXIS, POSITIVE_Z_AXIS),	// -X, +Z, +Y.
	std::make_pair(NEGATIVE_X_AXIS, NEGATIVE_Z_AXIS),	// -X, -Z, -Y.

	std::make_pair(POSITIVE_Y_AXIS, POSITIVE_Z_AXIS),	// +Y, +Z, +Y.
	std::make_pair(POSITIVE_Y_AXIS, NEGATIVE_Z_AXIS),	// +Y, -Z, -Y.
	std::make_pair(POSITIVE_Y_AXIS, POSITIVE_X_AXIS),	// +Y, +X, -Z.
	std::make_pair(POSITIVE_Y_AXIS, NEGATIVE_X_AXIS),	// +Y, -X, +Z.

	std::make_pair(NEGATIVE_Y_AXIS, POSITIVE_Z_AXIS),	// -Y, +Z, -Y.
	std::make_pair(NEGATIVE_Y_AXIS, NEGATIVE_Z_AXIS),	// -Y, -Z, +Y.
	std::make_pair(NEGATIVE_Y_AXIS, POSITIVE_X_AXIS),	// -Y, +X, +Z.
	std::make_pair(NEGATIVE_Y_AXIS, NEGATIVE_X_AXIS),	// -Y, -X, -Z.

	std::make_pair(POSITIVE_Z_AXIS, POSITIVE_X_AXIS),	// +Z, +X, +Y.
	std::make_pair(POSITIVE_Z_AXIS, NEGATIVE_X_AXIS),	// +Z, -X, -Y.
	std::make_pair(POSITIVE_Z_AXIS, POSITIVE_Y_AXIS),	// +Z, +Y, -X.
	std::make_pair(POSITIVE_Z_AXIS, NEGATIVE_Y_AXIS),	// +Z, -Y, +X.

	std::make_pair(NEGATIVE_Z_AXIS, POSITIVE_X_AXIS),	// -Z, +X, -Y.
	std::make_pair(NEGATIVE_Z_AXIS, NEGATIVE_X_AXIS),	// -Z, -X, +Y.
	std::make_pair(NEGATIVE_Z_AXIS, POSITIVE_Y_AXIS),	// -Z, +Y, +X.
	std::make_pair(NEGATIVE_Z_AXIS, NEGATIVE_Y_AXIS)	// -Z, -Y, -X.
};

MyMesh::Normal get_axis(const MeshCuboid::AXIS_DIRECTION _axis_direction)
{
	MyMesh::Normal axis;
	switch (_axis_direction)
	{
	case MeshCuboid::POSITIVE_X_AXIS:
		axis = MyMesh::Normal(1.0, 0.0, 0.0);
		break;
	case MeshCuboid::NEGATIVE_X_AXIS:
		axis = MyMesh::Normal(-1.0, 0.0, 0.0);
		break;
	case MeshCuboid::POSITIVE_Y_AXIS:
		axis = MyMesh::Normal(0.0, 1.0, 0.0);
		break;
	case MeshCuboid::NEGATIVE_Y_AXIS:
		axis = MyMesh::Normal(0.0, -1.0, 0.0);
		break;
	case MeshCuboid::POSITIVE_Z_AXIS:
		axis = MyMesh::Normal(0.0, 0.0, 1.0);
		break;
	case MeshCuboid::NEGATIVE_Z_AXIS:
		axis = MyMesh::Normal(0.0, 0.0, -1.0);
		break;
	default:
		assert(false);
		break;
	}
	return axis;
}

std::array<MyMesh::Normal, 3>
MeshCuboid::get_transformed_axes(const unsigned int _axis_configuration_index,
	const std::array<MyMesh::Normal, 3> &_axes)
{
	assert(_axis_configuration_index < num_axis_configurations());
	std::array<MyMesh::Normal, 3> rotation_axes;
	rotation_axes[0] = get_axis(k_all_axis_configuration[_axis_configuration_index].first);
	rotation_axes[1] = get_axis(k_all_axis_configuration[_axis_configuration_index].second);
	rotation_axes[2] = cross(rotation_axes[0], rotation_axes[1]);

	Eigen::Matrix3d rotation_mat;
	Eigen::Matrix3d axes_mat;

	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		for (unsigned int i = 0; i < 3; ++i)
		{
			rotation_mat.col(axis_index)[i] = rotation_axes[axis_index][i];
			axes_mat.row(axis_index)[i] = _axes[axis_index][i];
		}
	}

	Eigen::Matrix3d new_axes_mat = rotation_mat * axes_mat;
	std::array<MyMesh::Normal, 3> new_axes;

	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		for (unsigned int i = 0; i < 3; ++i)
		{
			new_axes[axis_index][i] = new_axes_mat.row(axis_index)[i];
		}
	}

	return new_axes;
}

void MeshCuboid::set_axis_configuration(const unsigned int _axis_configuration_index)
{
	assert(_axis_configuration_index < num_axis_configurations());
	std::array<MyMesh::Normal, 3> rotation_axes;
	rotation_axes[0] = get_axis(k_all_axis_configuration[_axis_configuration_index].first);
	rotation_axes[1] = get_axis(k_all_axis_configuration[_axis_configuration_index].second);
	rotation_axes[2] = cross(rotation_axes[0], rotation_axes[1]);

	Eigen::Matrix3d rotation_mat;
	Eigen::Matrix3d bbox_axes_mat;
	Eigen::Vector3d bbox_size_vec;

	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		for (unsigned int i = 0; i < 3; ++i)
		{
			rotation_mat.col(axis_index)[i] = rotation_axes[axis_index][i];
			bbox_axes_mat.row(axis_index)[i] = bbox_axes_[axis_index][i];
		}
		bbox_size_vec(axis_index) = bbox_size_[axis_index];
	}

	Eigen::Matrix3d new_bbox_axes_mat = rotation_mat * bbox_axes_mat;

	// NOTE:
	// For size, the +/- sign of axis should NOT be considered,
	// but just the order of axes should be changed.
	// The center position is NOT changed.
	Eigen::Matrix3d abs_rotation_mat = rotation_mat.array().cwiseAbs();
	Eigen::Vector3d new_bbox_size_vec = abs_rotation_mat * bbox_size_vec;

	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		for (unsigned int i = 0; i < 3; ++i)
		{
			bbox_axes_[axis_index][i] = new_bbox_axes_mat.row(axis_index)[i];
		}
		bbox_size_[axis_index] = new_bbox_size_vec(axis_index);
	}

	update_corner_points();
}

MeshCuboid::MeshCuboid(const LabelIndex _label_index)
	: label_index_(_label_index)
	, bbox_center_(0.0)
	, bbox_size_(0.0)
	//, is_group_cuboid_(false)
{
	bbox_axes_[0] = MyMesh::Normal(1.0, 0.0, 0.0);
	bbox_axes_[1] = MyMesh::Normal(0.0, 1.0, 0.0);
	bbox_axes_[2] = MyMesh::Normal(0.0, 0.0, 1.0);

	for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
		bbox_corners_[corner_index] = MyMesh::Point(0.0);
}

MeshCuboid::MeshCuboid(const MeshCuboid& _other)
{
	deep_copy(_other);
}

MeshCuboid& MeshCuboid::operator=(const MeshCuboid& _other)
{
	deep_copy(_other);
	return (*this);
}

void MeshCuboid::deep_copy(const MeshCuboid& _other)
{
	//this->is_group_cuboid_ = _other.is_group_cuboid_;
	this->label_index_ = _other.label_index_;

	// NOTE:
	// Do not deep copy sample points.
	this->sample_points_ = _other.sample_points_;

	this->sample_to_cuboid_surface_correspondence_ = _other.sample_to_cuboid_surface_correspondence_;
	this->cuboid_surface_to_sample_corresopndence_ = _other.cuboid_surface_to_sample_corresopndence_;

	this->bbox_axes_ = _other.bbox_axes_;
	this->bbox_center_ = _other.bbox_center_;
	this->bbox_size_ = _other.bbox_size_;
	this->bbox_corners_ = _other.bbox_corners_;

	// Deep copy cuboid surface points.
	assert(_other.cuboid_surface_points_.size() == _other.num_cuboid_surface_points());
	this->cuboid_surface_points_.clear();
	this->cuboid_surface_points_.reserve(_other.num_cuboid_surface_points());

	for (std::vector<MeshCuboidSurfacePoint *>::const_iterator it =
		_other.cuboid_surface_points_.begin();
		it != _other.cuboid_surface_points_.end(); ++it)
	{
		MeshCuboidSurfacePoint *cuboid_surface_point = new MeshCuboidSurfacePoint(**it);
		this->cuboid_surface_points_.push_back(cuboid_surface_point);
	}
}

MeshCuboid::~MeshCuboid()
{
	clear_cuboid_surface_points();
}

unsigned int MeshCuboid::num_sample_points()const {
	return static_cast<unsigned int>(sample_points_.size());
}

unsigned int MeshCuboid::num_cuboid_surface_points()const {
	return static_cast<unsigned int>(cuboid_surface_points_.size());
}

const std::vector<MeshSamplePoint *> &MeshCuboid::get_sample_points() const
{
	return sample_points_;
}

void MeshCuboid::get_sample_points(std::vector<MyMesh::Point> &_sample_points) const
{
	_sample_points.clear();
	_sample_points.reserve(num_sample_points());

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points();
		++sample_point_index)
	{
		MeshSamplePoint *sample_point = get_sample_point(sample_point_index);
		assert(sample_point);
		MyMesh::Point point = sample_point->point_;
		_sample_points.push_back(point);
	}
}

void MeshCuboid::get_sample_points(Eigen::MatrixXd &_sample_points) const
{
	_sample_points = Eigen::MatrixXd(3, num_sample_points());

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points();
		++sample_point_index)
	{
		MeshSamplePoint *sample_point = get_sample_point(sample_point_index);
		assert(sample_point);
		MyMesh::Point point = sample_point->point_;

		for (unsigned int i = 0; i < 3; ++i)
			_sample_points.col(sample_point_index)(i) = point[i];
	}
}

MeshSamplePoint *MeshCuboid::get_sample_point(
	const unsigned int _point_index) const
{
	assert(_point_index < sample_points_.size());
	return sample_points_[_point_index];
}

const std::vector<MeshCuboidSurfacePoint *>&
MeshCuboid::get_cuboid_surface_points() const
{
	return cuboid_surface_points_;
}

MeshCuboidSurfacePoint *MeshCuboid::get_cuboid_surface_point(
	const unsigned int _point_index) const
{
	assert(_point_index < cuboid_surface_points_.size());
	return cuboid_surface_points_[_point_index];
}

MyMesh::Point MeshCuboid::get_bbox_min() const
{
	MyMesh::Point bbox_min = bbox_center_;
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		bbox_min -= (0.5 * bbox_size_[axis_index]) * bbox_axes_[axis_index];
	}
	return bbox_min;
}

MyMesh::Point MeshCuboid::get_bbox_max() const
{
	MyMesh::Point bbox_max = bbox_center_;
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		bbox_max += (0.5 * bbox_size_[axis_index]) * bbox_axes_[axis_index];
	}
	return bbox_max;
}

MyMesh::Normal MeshCuboid::get_bbox_axis(const unsigned int _axis_index) const
{
	assert(_axis_index < 3);
	return bbox_axes_[_axis_index];
}

std::array<MyMesh::Normal, 3> MeshCuboid::get_bbox_axes() const
{
	return bbox_axes_;
}

MyMesh::Point MeshCuboid::get_bbox_center() const
{
	return bbox_center_;
}

MyMesh::Normal MeshCuboid::get_bbox_size() const
{
	return bbox_size_;
}

MyMesh::Point MeshCuboid::get_bbox_corner(const unsigned int _corner_index,
	MyMesh::Normal *_axis_direction) const
{
	assert(_corner_index < k_num_corners);

	// Optional.
	if (_axis_direction)
	{
		std::bitset<3> bits(_corner_index);
		Eigen::Vector3d axis_directions_vec = Eigen::Vector3d::Zero();

		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		{
			if (bits[axis_index])	axis_directions_vec[axis_index] = 1.0;
			else					axis_directions_vec[axis_index] = -1.0;
		}

		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
			(*_axis_direction)[axis_index] = axis_directions_vec[axis_index];
	}

	return bbox_corners_[_corner_index];
}

std::array<MyMesh::Point, MeshCuboid::k_num_corners> MeshCuboid::get_bbox_corners() const
{
	return bbox_corners_;
}

std::pair<MyMesh::Normal, MyMesh::Point>
MeshCuboid::get_bbox_faces(const unsigned int _face_index,
	MyMesh::Normal *_axis_direction) const
{
	const unsigned int num_faces = 6;
	assert(_face_index < num_faces);

	unsigned int face_axis_index = _face_index / 2;

	MyMesh::Normal face_normal;
	MyMesh::Point face_point;

	//if ((_face_index % 2) == 0)	face_normal = bbox_axes_[face_axis_index];
	//else  face_normal = -bbox_axes_[face_axis_index];

	//face_point = bbox_center_ + (0.5 * bbox_size_[face_axis_index]) * face_normal;

	Eigen::Vector3d axis_directions_vec = Eigen::Vector3d::Zero();
	Eigen::Vector3d bbox_center_vec;
	Eigen::Vector3d bbox_size_vec;
	Eigen::Matrix3d bbox_axes_mat;

	if ((face_axis_index % 2) == 0)
		axis_directions_vec[face_axis_index] = 1.0;
	else
		axis_directions_vec[face_axis_index] = -1.0;

	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		bbox_center_vec[axis_index] = bbox_center_[axis_index];
		bbox_size_vec[axis_index] = bbox_size_[axis_index];
		for (unsigned int i = 0; i < 3; i++)
			bbox_axes_mat.col(axis_index)[i] = bbox_axes_[axis_index][i];
	}

	Eigen::Vector3d face_normal_vec = bbox_axes_mat * axis_directions_vec;

	Eigen::Vector3d face_point_vec = bbox_center_vec +
		(bbox_axes_mat * (0.5 * axis_directions_vec).asDiagonal() * bbox_size_vec);

	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		face_normal[axis_index] = face_normal_vec[axis_index];
		face_point[axis_index] = face_point_vec[axis_index];

		// Optional.
		if (_axis_direction) (*_axis_direction)[axis_index] = axis_directions_vec[axis_index];
	}

	return std::make_pair(face_normal, face_point);
}

std::vector< std::pair<MyMesh::Normal, MyMesh::Point> > MeshCuboid::get_bbox_faces() const
{
	const unsigned int num_faces = 6;
	std::vector< std::pair<MyMesh::Normal, MyMesh::Point> > bbox_faces(num_faces);

	for (unsigned int face_index = 0; face_index < num_faces; ++face_index)
		bbox_faces[face_index] = get_bbox_faces(face_index);

	return bbox_faces;
}

Real MeshCuboid::get_bbox_volume() const
{
	return bbox_size_[0] * bbox_size_[1] * bbox_size_[2];
}

Real MeshCuboid::get_bbox_diag_length() const
{
	return bbox_size_.length();
}

Real MeshCuboid::get_bbox_face_area(const unsigned int _face_index) const
{
	// http://mathworld.wolfram.com/Quadrilateral.html.

	const unsigned int *corner_indices = k_face_corner_indices[_face_index];
	std::array<MyMesh::Point, k_num_face_corners> corner_point;
	for (unsigned int i = 0; i < k_num_face_corners; ++i)
		corner_point[i] = bbox_corners_[corner_indices[i]];

	MyMesh::Normal p = corner_point[2] - corner_point[0];
	MyMesh::Normal q = corner_point[3] - corner_point[1];
	Real area = 0.5 * (cross(p, q)).length();

	return area;
}

MyMesh::Point MeshCuboid::get_local_coord(const MyMesh::Point _pos) const
{
	// Global coordinates -> Local coordinates (Bounding box center and axes).
	MyMesh::Point pos = _pos - bbox_center_;
	MyMesh::Point local_pos;
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		local_pos[axis_index] = dot(bbox_axes_[axis_index], pos);

	MyMesh::Point global_pos = get_global_coord(local_pos);
	//CHECK_NUMERICAL_ERROR(__FUNCTION__, (global_pos - _pos).norm());

	return local_pos;
}

MyMesh::Point MeshCuboid::get_global_coord(const MyMesh::Point _pos) const
{
	// Local coordinates (Bounding box center and axes) -> Global coordinates.
	MyMesh::Point global_pos(0.0);
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		global_pos += (_pos[axis_index] * bbox_axes_[axis_index]);
	global_pos += bbox_center_;
	return global_pos;
}

void MeshCuboid::get_local_coord_sample_points(std::vector<MyMesh::Point> &_sample_points) const
{
	_sample_points.clear();
	_sample_points.resize(num_sample_points());

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points();
		++sample_point_index)
	{
		MeshSamplePoint *sample_point = get_sample_point(sample_point_index);
		assert(sample_point);
		MyMesh::Point point = sample_point->point_;
		MyMesh::Point local_coord = get_local_coord(point);
		_sample_points[sample_point_index] = local_coord;
	}
}

void MeshCuboid::get_local_coord_sample_points(Eigen::MatrixXd &_sample_points) const
{
	_sample_points = Eigen::MatrixXd(3, num_sample_points());

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points();
		++sample_point_index)
	{
		MeshSamplePoint *sample_point = get_sample_point(sample_point_index);
		assert(sample_point);
		MyMesh::Point point = sample_point->point_;
		MyMesh::Point local_coord = get_local_coord(point);

		for (unsigned int i = 0; i < 3; ++i)
			_sample_points.col(sample_point_index)(i) = local_coord[i];
	}
}

const std::vector<int>&
MeshCuboid::get_sample_to_cuboid_surface_correspondences() const
{
	return sample_to_cuboid_surface_correspondence_;
}

int MeshCuboid::get_sample_to_cuboid_surface_correspondences(
	const unsigned int _point_index) const
{
	if (num_cuboid_surface_points() == 0)
		return -1;

	assert(_point_index < num_sample_points());
	return sample_to_cuboid_surface_correspondence_[_point_index];
}

const std::vector<int>&
MeshCuboid::get_cuboid_surface_to_sample_correspondence() const
{
	return cuboid_surface_to_sample_corresopndence_;
}

int MeshCuboid::get_cuboid_surface_to_sample_correspondence(
	const unsigned int _point_index) const
{
	assert(_point_index < num_cuboid_surface_points());
	return cuboid_surface_to_sample_corresopndence_[_point_index];
}

void MeshCuboid::set_bbox_center(const MyMesh::Point &_bbox_center)
{
	bbox_center_ = _bbox_center;
}

void MeshCuboid::set_bbox_size(const MyMesh::Normal &_bbox_size,
	bool _update_corners)
{
	bbox_size_ = _bbox_size;

	if (_update_corners)
	{
		update_corner_points();
	}
}

void MeshCuboid::set_bbox_axes(const std::array<MyMesh::Normal, 3> &_bbox_axes,
	bool _update_corners)
{
	bbox_axes_ = _bbox_axes;

	if (_update_corners)
	{
		update_corner_points();
	}
}

void MeshCuboid::set_bbox_corners(const std::array<MyMesh::Point, k_num_corners> &_bbox_corners)
{
	bbox_corners_ = _bbox_corners;
}

void MeshCuboid::translate(const Eigen::Vector3d _translation_vec)
{
	for (unsigned int i = 0; i < 3; ++i)
	{
		bbox_center_[i] = bbox_center_[i] + _translation_vec[i];

		for (unsigned int corner_index = 0; corner_index < k_num_corners; ++corner_index)
			bbox_corners_[corner_index][i] = bbox_corners_[corner_index][i] + _translation_vec[i];
	}
}

void MeshCuboid::rotate(const Eigen::Matrix3d _rotation_mat, bool _update_center_size)
{
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		Eigen::Vector3d bbox_axis_vec;
		for (unsigned int i = 0; i < 3; ++i)
			bbox_axis_vec[i] = bbox_axes_[axis_index][i];

		bbox_axis_vec = _rotation_mat.transpose() * bbox_axis_vec;

		for (unsigned int i = 0; i < 3; ++i)
			bbox_axes_[axis_index][i] = bbox_axis_vec[i];

		bbox_axes_[axis_index].normalize();
	}

	if (_update_center_size)
	{
		assert(num_sample_points() >= FLAGS_param_min_num_cuboid_sample_points);

		Eigen::Matrix3d local_coord_rotation_mat;
		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
			for (unsigned int i = 0; i < 3; ++i)
				local_coord_rotation_mat.row(axis_index)(i) = bbox_axes_[axis_index][i];

		Eigen::MatrixXd sample_points_mat = Eigen::MatrixXd::Zero(3, num_sample_points());
		for (SamplePointIndex sapmle_point_index = 0; sapmle_point_index < num_sample_points();
			++sapmle_point_index)
		{
			const MyMesh::Point sample_point = sample_points_[sapmle_point_index]->point_;
			for (unsigned int i = 0; i < 3; ++i)
				sample_points_mat(i, sapmle_point_index) = sample_point[i];
		}

		Eigen::MatrixXd rotated_sample_points_mat = local_coord_rotation_mat * sample_points_mat;
		Eigen::Vector3d bbox_size_vec = rotated_sample_points_mat.rowwise().maxCoeff()
			- rotated_sample_points_mat.rowwise().minCoeff();
		Eigen::Vector3d bbox_center_vec = 0.5 * (rotated_sample_points_mat.rowwise().minCoeff()
			+ rotated_sample_points_mat.rowwise().maxCoeff());
		bbox_center_vec = local_coord_rotation_mat.transpose() * bbox_center_vec;
		

		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		{
			bbox_center_[axis_index] = bbox_center_vec[axis_index];
			bbox_size_[axis_index] = bbox_size_vec[axis_index];
		}
	}
}

void MeshCuboid::flip_axis(const unsigned int _axis_index)
{
	// Bounding box center and sizes are not changed.
	assert(_axis_index < 3);

	bbox_axes_[_axis_index] = -bbox_axes_[_axis_index];

	// Flip corner points.
	std::array<MyMesh::Point, k_num_corners> new_bbox_corners;
	for (unsigned int corner_index = 0; corner_index < k_num_corners; ++corner_index)
	{
		std::bitset<3> bits(corner_index);
		bits.flip(_axis_index);
		unsigned int new_corner_index = static_cast<unsigned int>(bits.to_ulong());
		new_bbox_corners[new_corner_index] = bbox_corners_[corner_index];
	}

	bbox_corners_.swap(new_bbox_corners);
}

void MeshCuboid::clear_sample_points()
{
	sample_points_.clear();
	sample_to_cuboid_surface_correspondence_.clear();
	cuboid_surface_to_sample_corresopndence_.clear();

	cuboid_surface_to_sample_corresopndence_.resize(num_cuboid_surface_points(), -1);
}

void MeshCuboid::add_sample_point(MeshSamplePoint *_point)
{
	if (sample_points_.size() == sample_to_cuboid_surface_correspondence_.size())
		sample_to_cuboid_surface_correspondence_.push_back(-1);
	sample_points_.push_back(_point);
}

void MeshCuboid::add_sample_points(const std::vector<MeshSamplePoint *> _points)
{
	if (sample_points_.size() == sample_to_cuboid_surface_correspondence_.size())
		sample_to_cuboid_surface_correspondence_.resize(
		sample_points_.size() + _points.size(), - 1);
	sample_points_.insert(sample_points_.end(), _points.begin(), _points.end());
}

void MeshCuboid::remove_sample_points(const bool *is_sample_point_removed)
{
	assert(is_sample_point_removed);

	for (std::vector<MeshSamplePoint *>::iterator it = sample_points_.begin();
		it != sample_points_.end(); )	// No increment.
	{
		assert(*it);
		if (is_sample_point_removed[(*it)->sample_point_index_])
			it = sample_points_.erase(it);
		else
			++it;
	}
}

void MeshCuboid::clear_cuboid_surface_points()
{
	for (std::vector<MeshCuboidSurfacePoint *>::iterator it = cuboid_surface_points_.begin();
		it != cuboid_surface_points_.end(); ++it)
		delete (*it);

	cuboid_surface_points_.clear();
}

void MeshCuboid::create_random_points_on_cuboid_surface(
	const unsigned int _num_cuboid_surface_points)
{
	clear_cuboid_surface_points();
	cuboid_surface_points_.reserve(_num_cuboid_surface_points);

	Real all_faces_area = 0;
	for (unsigned int face_index = 0; face_index < k_num_faces; ++face_index)
		all_faces_area += get_bbox_face_area(face_index);
	assert(all_faces_area > 0);

	// Sample points on each face.
	static SimpleRandomCong_t rng_cong;
	simplerandom_cong_seed(&rng_cong, CUBOID_SURFACE_SAMPLING_RANDOM_SEED);
	

	for (unsigned int face_index = 0; face_index < k_num_faces; ++face_index)
	{
		const unsigned int *corner_indices = k_face_corner_indices[face_index];
		std::array<MyMesh::Point, k_num_face_corners> corner_point;
		for (unsigned int i = 0; i < k_num_face_corners; ++i)
			corner_point[i] = bbox_corners_[corner_indices[i]];

		MyMesh::Normal normal = get_bbox_axis(face_index / 2);
		if ((face_index % 2) != 0)
			normal = -normal;

		Real face_area = get_bbox_face_area(face_index);
		int num_face_points = std::round(face_area / all_faces_area * _num_cuboid_surface_points);

		for (int point_index = 0; point_index < num_face_points
			&& cuboid_surface_points_.size() < _num_cuboid_surface_points; ++point_index)
		{
			Real w1 = static_cast<Real>(simplerandom_cong_next(&rng_cong))
				/ std::numeric_limits<uint32_t>::max();
			Real w2 = static_cast<Real>(simplerandom_cong_next(&rng_cong))
				/ std::numeric_limits<uint32_t>::max();

			MyMesh::Point p1 = w1 * (corner_point[1] - corner_point[0]) + corner_point[0];
			MyMesh::Point p2 = w1 * (corner_point[2] - corner_point[3]) + corner_point[3];
			MyMesh::Point point = w2 * (p2 - p1) + p1;

			Real sum_corner_weights = 0;
			std::array<Real, k_num_corners> corner_weights;
			corner_weights.fill(0.0);

			for (unsigned int i = 0; i < k_num_face_corners; ++i)
			{
				corner_weights[ corner_indices[i] ] = (point - corner_point[i]).length();
				sum_corner_weights += corner_weights[i];
			}

			assert(sum_corner_weights > 0);
			for (unsigned int i = 0; i < k_num_face_corners; ++i)
			{
				corner_weights[ corner_indices[i] ] =
					(sum_corner_weights - corner_weights[ corner_indices[i] ]) / sum_corner_weights;
			}

#ifdef DEBUG_TEST
			MyMesh::Point same_point(0.0);
			for (unsigned int i = 0; i < k_num_corners; ++i)
				same_point += corner_weights[i] * bbox_corners_[i];
			Real error = (same_point - point).length();
			CHECK_NUMERICAL_ERROR(__FUNCTION__, error);
#endif

			MeshCuboidSurfacePoint *cuboid_surface_point = new MeshCuboidSurfacePoint(
				point, normal, face_index, corner_weights);
			cuboid_surface_points_.push_back(cuboid_surface_point);
		}
	}


	// Update point correspondences.
	update_point_correspondences();
}

void MeshCuboid::create_grid_points_on_cuboid_surface(
	const unsigned int _num_cuboid_surface_points)
{
	clear_cuboid_surface_points();
	cuboid_surface_points_.reserve(_num_cuboid_surface_points);

	Real all_faces_area = 0;
	for (unsigned int face_index = 0; face_index < k_num_faces; ++face_index)
		all_faces_area += get_bbox_face_area(face_index);


	// NOTE:
	// The number of points may not be exactly the same with the given '_num_cuboid_surface_points'.
	for (unsigned int face_index = 0; face_index < k_num_faces; ++face_index)
	{
		const unsigned int *corner_indices = k_face_corner_indices[face_index];
		std::array<MyMesh::Point, k_num_face_corners> corner_point;
		for (unsigned int i = 0; i < k_num_face_corners; ++i)
			corner_point[i] = bbox_corners_[corner_indices[i]];

		MyMesh::Normal normal = get_bbox_axis(face_index / 2);
		if ((face_index % 2) != 0)
			normal = -normal;

		//
		Real face_area = get_bbox_face_area(face_index);
		Real num_face_points = (face_area / all_faces_area * _num_cuboid_surface_points);

		Real axis_length[2];
		axis_length[0] = (corner_point[1] - corner_point[0]).length()
			+ (corner_point[2] - corner_point[3]).length();
		axis_length[1] = (corner_point[3] - corner_point[0]).length()
			+ (corner_point[2] - corner_point[1]).length();

		if (axis_length[0] <= MIN_CUBOID_SIZE || axis_length[1] <= MIN_CUBOID_SIZE)
			continue;

		Real ratio = std::sqrt(num_face_points / (axis_length[0] * axis_length[1]));

		unsigned int num_axis_points[2];
		for (unsigned int i = 0; i < 2; ++i)
		{
			num_axis_points[i] = static_cast<unsigned int>(std::round(ratio * axis_length[i]));
			num_axis_points[i] = std::min(num_axis_points[i], static_cast<unsigned int>(num_face_points / 2));
			num_axis_points[i] = std::max(num_axis_points[i], 2u);
		}
		//

		for (unsigned int point_index_1 = 0; point_index_1 < num_axis_points[0]; ++point_index_1)
		{
			Real w1 = static_cast<Real>(point_index_1) / (num_axis_points[0] - 1);

			for (unsigned int point_index_2 = 0; point_index_2 < num_axis_points[1]; ++point_index_2)
			{
				Real w2 = static_cast<Real>(point_index_2) / (num_axis_points[1] - 1);

				MyMesh::Point p1 = w1 * (corner_point[1] - corner_point[0]) + corner_point[0];
				MyMesh::Point p2 = w1 * (corner_point[2] - corner_point[3]) + corner_point[3];
				MyMesh::Point point = w2 * (p2 - p1) + p1;

				Real sum_corner_weights = 0;
				std::array<Real, k_num_corners> corner_weights;
				corner_weights.fill(0.0);

				//corner_weights[0] = (1 - w1)*(1 - w2);
				//corner_weights[1] = (w1)*(1 - w2);
				//corner_weights[2] = (w1)*(w2);
				//corner_weights[3] = (1 - w1)*(w2);

				corner_weights[ corner_indices[0] ] = (1 - w1)*(1 - w2);
				corner_weights[ corner_indices[1] ] = (w1)*(1 - w2);
				corner_weights[ corner_indices[2] ] = (w1)*(w2);
				corner_weights[ corner_indices[3] ] = (1 - w1)*(w2);

#ifdef DEBUG_TEST
				MyMesh::Point same_point(0.0);
				for (unsigned int i = 0; i < k_num_corners; ++i)
					same_point += corner_weights[i] * bbox_corners_[i];
				Real error = (same_point - point).length();
				CHECK_NUMERICAL_ERROR(__FUNCTION__, error);
#endif

				MeshCuboidSurfacePoint *cuboid_surface_point = new MeshCuboidSurfacePoint(
					point, normal, face_index, corner_weights);
				cuboid_surface_points_.push_back(cuboid_surface_point);
			}
		}
	}

	// Update point correspondences.
	update_point_correspondences();
}

void MeshCuboid::update_point_correspondences()
{
	sample_to_cuboid_surface_correspondence_.clear();
	cuboid_surface_to_sample_corresopndence_.clear();

	sample_to_cuboid_surface_correspondence_.resize(num_sample_points(), -1);
	cuboid_surface_to_sample_corresopndence_.resize(num_cuboid_surface_points(), -1);


	// NOTE:
	// X: sample points, Y: cuboid surface_points. 
	unsigned int num_X_points = num_sample_points();
	unsigned int num_Y_points = num_cuboid_surface_points();

	if (num_X_points == 0 || num_Y_points == 0)
		return;

	Eigen::MatrixXd X_points(3, num_X_points);
	Eigen::MatrixXd Y_points(3, num_Y_points);

	// FIXME:
	// The type of indices should integer.
	// But, it causes compile errors in the 'ICP::get_closest_points' function.
	Eigen::MatrixXd X_indices(1, num_X_points);
	Eigen::MatrixXd Y_indices(1, num_Y_points);

	for (unsigned int X_point_index = 0; X_point_index < num_X_points; ++X_point_index)
	{
		for (unsigned int i = 0; i < 3; ++i)
			X_points.col(X_point_index)(i) = get_sample_points()[X_point_index]->point_[i];
		X_indices.col(X_point_index)(0) = X_point_index;
	}

	for (unsigned int Y_point_index = 0; Y_point_index < num_Y_points; ++Y_point_index)
	{
		for (unsigned int i = 0; i < 3; ++i)
			Y_points.col(Y_point_index)(i) = get_cuboid_surface_points()[Y_point_index]->point_[i];
		Y_indices.col(Y_point_index)(0) = Y_point_index;
	}

	ANNpointArray X_ann_points = NULL;
	ANNkd_tree* X_ann_kd_tree = ICP::create_kd_tree(X_points, X_ann_points);
	assert(X_ann_kd_tree);

	ANNpointArray Y_ann_points = NULL;
	ANNkd_tree* Y_ann_kd_tree = ICP::create_kd_tree(Y_points, Y_ann_points);
	assert(Y_ann_kd_tree);

	// X -> Y.
	Eigen::MatrixXd closest_Y_indices;
	ICP::get_closest_points(Y_ann_kd_tree, X_points, Y_indices, closest_Y_indices);
	assert(closest_Y_indices.cols() == num_X_points);

	// Y -> X.
	Eigen::MatrixXd closest_X_indices;
	ICP::get_closest_points(X_ann_kd_tree, Y_points, X_indices, closest_X_indices);
	assert(closest_X_indices.cols() == num_Y_points);


	// NOTE:
	// X: sample points, Y: cuboid surface_points.
	for (unsigned int X_point_index = 0; X_point_index < num_X_points; ++X_point_index)
	{
		assert(closest_Y_indices.col(X_point_index)(0) < num_Y_points);
		sample_to_cuboid_surface_correspondence_[X_point_index] =
			static_cast<int>(closest_Y_indices.col(X_point_index)(0));
	}

	for (unsigned int Y_point_index = 0; Y_point_index < num_Y_points; ++Y_point_index)
	{
		assert(closest_X_indices.col(Y_point_index)(0) < num_X_points);
		cuboid_surface_to_sample_corresopndence_[Y_point_index] =
			static_cast<int>(closest_X_indices.col(Y_point_index)(0));
	}

	if (X_ann_points) annDeallocPts(X_ann_points);
	if (Y_ann_points) annDeallocPts(Y_ann_points);
	delete X_ann_kd_tree;
	delete Y_ann_kd_tree;
}

void MeshCuboid::compute_cuboid_surface_point_visibility(
	const Real _modelview_matrix[16],
	const Real _radius,
	const std::vector<MeshSamplePoint *> &_given_sample_points,
	const std::vector<MyMesh::Point> &_test_points,
	const std::vector<MyMesh::Normal> *_test_normals,
	std::vector<Real> &_visibility_values)
{
	assert(_modelview_matrix);
	assert(_radius > 0);

	Eigen::Matrix4d modelview_matrix;
	for (unsigned int col = 0; col < 4; ++col)
		for (unsigned int row = 0; row < 4; ++row)
			modelview_matrix(row, col) = _modelview_matrix[4 * col + row];

	MyMesh::Normal view_direction(
		-_modelview_matrix[2], -_modelview_matrix[6], -_modelview_matrix[10]);

	unsigned int num_test_points = _test_points.size();
	assert(!_test_normals || (*_test_normals).size() == num_test_points);
	_visibility_values.resize(num_test_points);

	for (unsigned int test_point_index = 0; test_point_index < num_test_points; ++test_point_index)
	{
		Real &visibility = _visibility_values[test_point_index];
		visibility = 1.0;

		// Ignore a surface point if its normal is not heading to the viewing direction.
		if (_test_normals && dot((*_test_normals)[test_point_index], view_direction) >= 0)
		{
			visibility = 0.0;
			continue;
		}

		Eigen::Vector3d surface_point;
		for (unsigned int i = 0; i < 3; ++i)
			surface_point[i] = _test_points[test_point_index][i];

		Eigen::Vector4d surface_point_4;
		surface_point_4 << surface_point, 1.0;

		for (std::vector<MeshSamplePoint *>::const_iterator jt = _given_sample_points.begin();
			jt != _given_sample_points.end(); jt++)
		{
			Eigen::Vector3d observed_point;
			for (unsigned int i = 0; i < 3; ++i)
				observed_point[i] = (*jt)->point_[i];

			Eigen::Vector4d observed_point_4;
			observed_point_4 << observed_point, 1.0;

			// Positions in the model view coordinates.
			Eigen::Vector4d lc_observed_point_4 = modelview_matrix * observed_point_4;
			Eigen::Vector3d lc_observed_point = lc_observed_point_4.topRows(3) / lc_observed_point_4[3];
			lc_observed_point[2] -= _radius;

			Eigen::Vector4d lc_surface_point_4 = modelview_matrix * surface_point_4;
			Eigen::Vector3d lc_surface_point = lc_surface_point_4.topRows(3) / lc_surface_point_4[3];

			Real lc_surface_point_len = lc_surface_point.norm();
			Real lc_observed_point_len = lc_observed_point.norm();

			Real dot_prod = lc_surface_point.dot(lc_observed_point);

			// Ignore a sample point if either it or the voxel point is not visible
			// from this model view.
			if (lc_observed_point[2] >= 0 || lc_surface_point[2] >= 0
				|| lc_observed_point[2] <= lc_surface_point[2]
				|| lc_surface_point_len == 0 || lc_observed_point_len == 0)
				continue;

			Real cos_angle = dot_prod / (lc_surface_point_len * lc_observed_point_len);
			if (std::acos(cos_angle) <= std::atan(_radius / lc_observed_point_len))
			{
				visibility = 0;
			}

			//Real cos_angle = dot_prod / (lc_surface_point_len * lc_observed_point_len);
			//Real sin_angle = std::sqrt(1 - cos_angle * cos_angle);
			//Real tan_angle = sin_angle / cos_angle;

			//if (dot_prod < lc_observed_point_len * lc_observed_point_len)
			//{
			//	// NOTE:
			//	// If the following code is uncommented, a point is considered as a sphere.
			//	//Real distance = (surface_point - observed_point).norm() / _radius;
			//	//visibility = min(frag_color.a, distance);
			//}
			//else
			//{
			//	Real distance = (lc_observed_point_len * tan_angle) / _radius;
			//	//visibility = min(visibility, distance);
			//	// NOTE:
			//	//
			//	if (distance < 1) visibility = 0;
			//}
		}

		assert(visibility >= 0.0);
		assert(visibility <= 1.0);
	}


	// Test 2D view plane mask for occlusion.
	//
	if (FLAGS_use_view_plane_mask)
	{
		std::list<SamplePointIndex> occluded_test_point_indices;
		MeshCuboid::compute_view_plane_mask_visibility(_modelview_matrix,
			_test_points, occluded_test_point_indices);

		for (std::list<SamplePointIndex>::iterator it = occluded_test_point_indices.begin();
			it != occluded_test_point_indices.end(); ++it)
		{
			SamplePointIndex test_point_index = *it;
			assert(test_point_index < _visibility_values.size());
			_visibility_values[test_point_index] = 0.0;
		}
	}
	//
}

void MeshCuboid::compute_view_plane_mask_visibility(const Real _modelview_matrix[16],
	const std::vector<MyMesh::Point>& _points,
	std::list<SamplePointIndex> &_masked_point_indices)
{
	Eigen::Matrix3d rotation_mat;
	Eigen::Vector3d translation_vec;

	for (int row = 0; row < 3; ++row)
		for (int col = 0; col < 3; ++col)
			rotation_mat(row, col) = _modelview_matrix[4 * col + row];

	for (int row = 0; row < 3; ++row)
		translation_vec(row) = _modelview_matrix[4 * 3 + row];

	unsigned int num_points = _points.size();
	if (num_points == 0)
		return;

	Eigen::MatrixXd view_plane_points_mat = Eigen::MatrixXd(3, num_points);
	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_points;
		++sample_point_index)
	{
		Eigen::Vector3d point_vec;
		for (int i = 0; i < 3; ++i)
			point_vec(i) = _points[sample_point_index][i];

		// Transformation.
		point_vec = rotation_mat * point_vec + translation_vec;
		view_plane_points_mat.col(sample_point_index) = point_vec;
	}

	Eigen::Vector2d view_plane_mask_range_min;
	Eigen::Vector2d view_plane_mask_range_max;
	view_plane_mask_range_min << FLAGS_param_view_plane_mask_min_x, FLAGS_param_view_plane_mask_min_y;
	view_plane_mask_range_max << FLAGS_param_view_plane_mask_max_x, FLAGS_param_view_plane_mask_max_y;

	
	_masked_point_indices.clear();
	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_points;
		++sample_point_index)
	{
		Eigen::Vector3d sample_point_vec = view_plane_points_mat.col(sample_point_index);
		if (sample_point_vec[0] >= view_plane_mask_range_min[0]
			&& sample_point_vec[1] >= view_plane_mask_range_min[1]
			&& sample_point_vec[0] <= view_plane_mask_range_max[0]
			&& sample_point_vec[1] <= view_plane_mask_range_max[1])
			_masked_point_indices.push_back(sample_point_index);
	}
}

void MeshCuboid::compute_cuboid_surface_point_visibility(
	const Real _modelview_matrix[16],
	const Real _radius,
	const std::vector<MeshSamplePoint *>& _given_sample_points,
	bool _use_cuboid_normal)
{
	std::vector<MyMesh::Point> test_points(num_cuboid_surface_points());
	std::vector<MyMesh::Normal> test_normals(num_cuboid_surface_points());
	
	for (unsigned int point_index = 0; point_index < num_cuboid_surface_points(); ++point_index)
	{
		MeshCuboidSurfacePoint *cuboid_surface_point = cuboid_surface_points_[point_index];
		assert(cuboid_surface_point);
		test_points[point_index] = cuboid_surface_point->point_;
		test_normals[point_index] = cuboid_surface_point->normal_;
	}

	std::vector<Real> visibility_values;

	if (_use_cuboid_normal)
		compute_cuboid_surface_point_visibility(_modelview_matrix, _radius,
		_given_sample_points, test_points, &test_normals, visibility_values);
	else
		compute_cuboid_surface_point_visibility(_modelview_matrix, _radius,
		_given_sample_points, test_points, NULL, visibility_values);

	assert(visibility_values.size() == num_cuboid_surface_points());

	for (unsigned int point_index = 0; point_index < num_cuboid_surface_points(); ++point_index)
	{
		MeshCuboidSurfacePoint *cuboid_surface_point = cuboid_surface_points_[point_index];
		assert(cuboid_surface_point);
		cuboid_surface_point->visibility_ = visibility_values[point_index];
	}
}

Real MeshCuboid::get_cuboid_overvall_visibility() const
{
	// Assume that visibility of each surface point is already computed.
	Real sum_visibility = 0.0;
	if (cuboid_surface_points_.empty())
		return sum_visibility;

	for (std::vector<MeshCuboidSurfacePoint *>::const_iterator it = cuboid_surface_points_.begin();
		it != cuboid_surface_points_.end(); it++)
	{
		sum_visibility += (*it)->visibility_;
	}

	sum_visibility /= static_cast<Real>(cuboid_surface_points_.size());
	return sum_visibility;
}

bool MeshCuboid::is_point_inside_cuboid(const MyMesh::Point& _point) const
{
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		MyMesh::Normal axis = bbox_axes_[axis_index].normalized();
		if (std::abs(dot(axis, _point - bbox_center_)) > 0.5 * bbox_size_[axis_index])
			return false;
	}

	return true;
}

void MeshCuboid::update_corner_points()
{
	for (unsigned int corner_index = 0; corner_index < k_num_corners; ++corner_index)
	{
		std::bitset<3> bits(corner_index);

		//bbox_corners_[corner_index] = bbox_center_;
		//for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		//{
		//	if (bits[axis_index])
		//		bbox_corners_[corner_index] += (0.5 * bbox_size_[axis_index]) * bbox_axes_[axis_index];
		//	else
		//		bbox_corners_[corner_index] -= (0.5 * bbox_size_[axis_index]) * bbox_axes_[axis_index];
		//}

		Eigen::Vector3d axis_directions_vec = Eigen::Vector3d::Zero();
		Eigen::Vector3d bbox_center_vec;
		Eigen::Vector3d bbox_size_vec;
		Eigen::Matrix3d bbox_axes_mat;

		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		{
			if (bits[axis_index])	axis_directions_vec[axis_index] = 1.0;
			else					axis_directions_vec[axis_index] = -1.0;

			bbox_center_vec[axis_index] = bbox_center_[axis_index];
			bbox_size_vec[axis_index] = bbox_size_[axis_index];
			for (unsigned int i = 0; i < 3; i++)
				bbox_axes_mat.col(axis_index)[i] = bbox_axes_[axis_index][i];
		}

		Eigen::Vector3d bbox_corner_vec = bbox_center_vec +
			(bbox_axes_mat * (0.5 * axis_directions_vec).asDiagonal() * bbox_size_vec);

		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
			bbox_corners_[corner_index][axis_index] = bbox_corner_vec[axis_index];
	}
}

void MeshCuboid::update_center_size_corner_points()
{
	// error = (center + sum_i( +/-0.5 * scale_i * axis_i)) - corner.
	// Find center and size minimizing the error.

	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3 * MeshCuboid::k_num_corners, 6);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(3 * MeshCuboid::k_num_corners);

	for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
	{
		std::bitset<3> bits(corner_index);

		A.block<3, 3>(3 * corner_index, 0).setIdentity();

		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		{
			for (unsigned int i = 0; i < 3; ++i)
			{
				if (!bits[axis_index])
					A(3 * corner_index + i, 3 + axis_index) = -0.5 * bbox_axes_[axis_index][i];
				else
					A(3 * corner_index + i, 3 + axis_index) = +0.5 * bbox_axes_[axis_index][i];

				b(3 * corner_index + i) = bbox_corners_[corner_index][i];
			}
		}
	}

	// Least square.
	Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
	assert(x.rows() == 6);
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		bbox_center_[axis_index] = x(axis_index);
		bbox_size_[axis_index] = x(3 + axis_index);
	}

	update_corner_points();
}

void MeshCuboid::update_axes_center_size_corner_points()
{
	const unsigned int num_face_corners = MeshCuboid::k_num_face_corners;
	Eigen::Matrix3d axes;

	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		Eigen::MatrixXd A(num_face_corners, 3);
		double sum_length = 0.0;

		for (unsigned int face_corner_index = 0; face_corner_index < k_num_face_corners; ++face_corner_index)
		{
			std::bitset<3> bits;
			bits[(axis_index + 1) % 3] = ((face_corner_index / 2) == 0);
			bits[(axis_index + 2) % 3] = ((face_corner_index % 2) == 0);

			bits[axis_index] = true;
			unsigned int pos_corner_index = static_cast<unsigned int>(bits.to_ulong());
			MyMesh::Point pos_corner_point = bbox_corners_[pos_corner_index];

			bits[axis_index] = false;
			unsigned int neg_corner_index = static_cast<unsigned int>(bits.to_ulong());
			MyMesh::Point neg_corner_point = bbox_corners_[neg_corner_index];

			MyMesh::Normal direction = pos_corner_point - neg_corner_point;
			sum_length += direction.norm();
			direction.normalize();

			for (unsigned int i = 0; i < 3; ++i)
				A.row(face_corner_index)[i] = direction[i];
		}

		if (sum_length == 0)
			sum_length += MIN_CUBOID_SIZE;

		// The first column of matrix V in SVD is the vector maximizing
		// cos(angle) with all columns in the given matrix A.
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
		axes.col(axis_index) = svd.matrixV().col(0);

		// Weight each axis based on the length of cuboid edges.
		axes.col(axis_index) = sum_length * axes.col(axis_index);

		// FIXME:
		// Hack. Check whether each axis is flipped...
		Eigen::Vector3d axis;
		for (unsigned int i = 0; i < 3; ++i)
			axis[i] = bbox_axes_[axis_index][i];

		if (axes.col(axis_index).dot(axis) < 0)
			axes.col(axis_index) = -axes.col(axis_index);
	}


	// Find the nearest orthogonal matrix.
	Eigen::JacobiSVD<Eigen::Matrix3d> svd(axes, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::Matrix3d new_axes = svd.matrixU() * svd.matrixV().transpose();
	new_axes.colwise().normalize();

	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		for (unsigned int i = 0; i < 3; ++i)
			bbox_axes_[axis_index][i] = new_axes.col(axis_index)(i);


	update_center_size_corner_points();
}

bool MeshCuboid::compute_bbox()
{
	if (num_sample_points() < FLAGS_param_min_num_cuboid_sample_points)
		return false;

#ifdef AXIS_ALIGNED_INITIAL_CUBOID
	compute_axis_aligned_bbox();
#else
	compute_oriented_bbox();
#endif

	// Compute box corners.
	update_corner_points();

	return true;
}

void MeshCuboid::compute_oriented_bbox()
{
	assert(num_sample_points() >= FLAGS_param_min_num_cuboid_sample_points);

	Eigen::MatrixXd sample_points_mat(3, num_sample_points());
	for (SamplePointIndex sapmle_point_index = 0; sapmle_point_index < num_sample_points();
		++sapmle_point_index)
	{
		const MyMesh::Point sample_point = sample_points_[sapmle_point_index]->point_;
		for (unsigned int i = 0; i < 3; ++i)
			sample_points_mat(i, sapmle_point_index) = sample_point[i];
	}

	Eigen::Vector3d sample_point_mean = sample_points_mat.rowwise().mean();
	Eigen::MatrixXd centered_sample_points_mat = sample_points_mat.colwise()
		- sample_point_mean;

	Eigen::Matrix3d rotation_mat;
	Eigen::Vector3d translation_vec;

	// Initialize (PCA).
	//Eigen::MatrixXd cov = centered_sample_points_mat * centered_sample_points_mat.adjoint();
	//Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(centered_sample_points_mat, Eigen::ComputeFullU);


	// Note:
	// The eigenvalues are sorted in increasing order.
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		//rotation_mat.row(axis_index) = eig.eigenvectors().col(3 - axis_index - 1);
		rotation_mat.row(axis_index) = svd.matrixU().col(axis_index);

	// Fix z-axis direction.
	if ((rotation_mat.row(0).cross(rotation_mat.row(1))).dot(rotation_mat.row(2)) < 0)
		rotation_mat.row(2) = -rotation_mat.row(2);

	// Set identity first, and rotate.
	bbox_axes_[0] = MyMesh::Normal(1.0, 0.0, 0.0);
	bbox_axes_[1] = MyMesh::Normal(0.0, 1.0, 0.0);
	bbox_axes_[2] = MyMesh::Normal(0.0, 0.0, 1.0);

	rotate(rotation_mat, true);

	// Rotate by each axis and find the minimum size of bounding box.
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		Eigen::Vector3d axis_vec;
		for (unsigned int i = 0; i < 3; ++i)
			axis_vec[i] = bbox_axes_[axis_index][i];

		const double angle_interval = 1.0 / 180.0 * M_PI;
		double min_angle = 0.0;
		double min_bbox_size = std::numeric_limits<double>::max();

		for (double angle = -M_PI_4; angle <= M_PI_4; angle += angle_interval)
		{
			MeshCuboid copy_cuboid(*this);
			Eigen::AngleAxisd axis_rotation(angle, axis_vec);
			copy_cuboid.rotate(axis_rotation.toRotationMatrix(), true);

			Real bbox_size = 1.0;
			for (unsigned int i = 0; i < 3; ++i)
				if (i != axis_index)
					bbox_size *= copy_cuboid.get_bbox_size()[i];
			
			if (bbox_size < min_bbox_size)
			{
				min_angle = angle;
				min_bbox_size = bbox_size;
			}
		}

		Eigen::AngleAxisd axis_rotation(min_angle, axis_vec);
		rotate(axis_rotation.toRotationMatrix(), true);
	}

	//// ICP-style iterative optimization.
	//double prev_error = std::numeric_limits<double>::max();

	//const unsigned int max_num_iterations = 1;
	//const double min_angle_difference = static_cast<double>(MIN_BBOX_ORIENTATION_ANGLE_DIFFERENCE) / 180.0 * M_PI;
	//const double min_translation = MIN_BBOX_ORIENTATION_TRANSLATION;

	//std::vector<MyMesh::Point> cuboid_surface_points;
	//Eigen::MatrixXd random_surface_points_mat(3, 1000);

	//for (unsigned int iteration = 0; iteration < max_num_iterations; ++iteration)
	//{
	//	cuboid_surface_points.clear();
	//	create_grid_points_on_cuboid_surface();
	//	for (SamplePointIndex point_index = 0; point_index < random_surface_points_mat.cols(); ++point_index)
	//		for (unsigned int i = 0; i < 3; ++i)
	//			random_surface_points_mat(i, point_index) = cuboid_surface_points[point_index][i];

	//	double error = ICP::run_iterative_closest_points(sample_points_mat, random_surface_points_mat,
	//		rotation_mat, translation_vec);

	//	if (error > prev_error)	break;
	//	prev_error = error;

	//	Eigen::AngleAxisd rotation_angle;
	//	rotation_angle.fromRotationMatrix(rotation_mat);
	//	if (rotation_angle.angle() < min_angle_difference && translation_vec.norm() < min_translation)
	//		break;

	//	rotate_cuboid(rotation_mat, true);
	//}


	// Find the best order of bounding box axes so that each axis in order 0, 1, and 2
	// is close to x, y, and z axis, respectively.
	Real min_angle = std::numeric_limits<Real>::max();
	unsigned int best_axis_configuration_index = 0;

	for (unsigned int axis_configuration_axis = 0; axis_configuration_axis < num_axis_configurations();
		++axis_configuration_axis)
	{
		std::array<MyMesh::Normal, 3> curr_bbox_axes =
			get_transformed_axes(axis_configuration_axis, bbox_axes_);

		// Make a rotation matrix using current bounding box axes,
		// and compute the angle of rotation (rotation from identity matrix).
		Eigen::Matrix3d new_rotation_mat;
		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
			for (unsigned int i = 0; i < 3; ++i)
				new_rotation_mat.col(axis_index)(i) = curr_bbox_axes[axis_index][i];

		Eigen::AngleAxisd rotation_difference;
		rotation_difference.fromRotationMatrix(new_rotation_mat);

		if (rotation_difference.angle() < min_angle)
		{
			min_angle = rotation_difference.angle();
			best_axis_configuration_index = axis_configuration_axis;
		}
	}

	set_axis_configuration(best_axis_configuration_index);
}

void MeshCuboid::compute_axis_aligned_bbox()
{
	assert(num_sample_points() >= FLAGS_param_min_num_cuboid_sample_points);

	bbox_axes_[0] = MyMesh::Normal(1.0, 0.0, 0.0);
	bbox_axes_[1] = MyMesh::Normal(0.0, 1.0, 0.0);
	bbox_axes_[2] = MyMesh::Normal(0.0, 0.0, 1.0);

	MyMesh::Point bbox_min(+std::numeric_limits<Real>::max());
	MyMesh::Point bbox_max(-std::numeric_limits<Real>::max());

	for (SamplePointIndex sapmle_point_index = 0; sapmle_point_index < num_sample_points();
		++sapmle_point_index)
	{
		const MyMesh::Point sample_point = sample_points_[sapmle_point_index]->point_;
		for (unsigned int i = 0; i < 3; ++i)
		{
			bbox_min[i] = std::min(bbox_min[i], sample_point[i]);
			bbox_max[i] = std::max(bbox_max[i], sample_point[i]);
		}
	}

	bbox_center_ = 0.5 * (bbox_min + bbox_max);
	bbox_size_ = bbox_max - bbox_min;
}

void MeshCuboid::update_label_using_sample_points()
{
	// Update the label based on the label confidence values of sample points.
	if (num_sample_points() == 0)
		return;
	assert(!sample_points_.empty());

	unsigned int num_labels = sample_points_.front()->label_index_confidence_.size();
	assert(num_labels > 0);

	std::vector<Real> accumulated_label_confidence(num_labels, 0.0);

	for (std::vector<MeshSamplePoint *>::iterator it = sample_points_.begin();
		it != sample_points_.end(); ++it)
	{
		const std::vector<Real>& label_confidence = (*it)->label_index_confidence_;
		assert(label_confidence.size() == num_labels);

		for (LabelIndex label_index = 0; label_index < num_labels; ++label_index)
		{
			accumulated_label_confidence[label_index] += label_confidence[label_index];
		}
	}

	LabelIndex new_label_index = 0;
	for (LabelIndex label_index = 1; label_index < num_labels; ++label_index)
	{
		if (accumulated_label_confidence[label_index]
			> accumulated_label_confidence[new_label_index])
		{
			new_label_index = label_index;
		}
	}

	label_index_ = new_label_index;
}

std::vector<MeshCuboid *> MeshCuboid::split_cuboid(const Real _object_diameter)
{
	std::vector<MeshCuboid *> sub_cuboids;

	// Construct a KD-tree.
	const int dim = 3;
	ANNpointArray data_pts = annAllocPts(num_sample_points(), dim);	// allocate data points

	for (SamplePointIndex sample_point_index = 0; sample_point_index < num_sample_points();
		++sample_point_index)
	{
		for (int i = 0; i < dim; i++)
			data_pts[sample_point_index][i] = sample_points_[sample_point_index]->point_[i];
	}

	ANNkd_tree *kd_tree = new ANNkd_tree(data_pts, num_sample_points(), dim);
	
	create_sub_cuboids(_object_diameter, kd_tree, sub_cuboids);
	remove_small_sub_cuboids(sub_cuboids);
	//align_sub_cuboids(_object_diameter, sub_cuboids);

	if (data_pts) annDeallocPts(data_pts);
	delete kd_tree;
	annClose();

	return sub_cuboids;
}

void MeshCuboid::create_sub_cuboids(const Real _object_diameter,
	ANNkd_tree* _kd_tree, std::vector<MeshCuboid *> &_sub_cuboids)
{
	const int dim = 3;
	const Real squared_neighbor_distance = FLAGS_param_cuboid_split_neighbor_distance *
		FLAGS_param_cuboid_split_neighbor_distance * _object_diameter;
	const int num_neighbors = std::min(FLAGS_param_num_sample_point_neighbors,
		static_cast<int>(num_sample_points()));
	
	ANNpoint q = annAllocPt(dim);
	ANNidxArray nn_idx = new ANNidx[num_neighbors];
	ANNdistArray dd = new ANNdist[num_neighbors];

	bool *is_sample_visited = new bool[num_sample_points()];
	memset(is_sample_visited, false, num_sample_points() * sizeof(bool));


	for (std::vector<MeshCuboid *>::iterator it = _sub_cuboids.begin(); it != _sub_cuboids.end(); ++it)
		delete (*it);
	_sub_cuboids.clear();


	while (true)
	{
		std::vector<MeshSamplePoint *> sub_cuboid_sample_points;
		sub_cuboid_sample_points.reserve(sample_points_.size());

		// 1. Find the seed point that has the highest confidence.
		int seed_sample_point_index = -1;

		for (SamplePointIndex sample_point_index = 0;
			sample_point_index < num_sample_points(); sample_point_index++)
		{
			if (is_sample_visited[sample_point_index])
			{
				continue;
			}
			else if (seed_sample_point_index < 0)
			{
				seed_sample_point_index = sample_point_index;
				continue;
			}
			assert(seed_sample_point_index >= 0);

			MeshSamplePoint* sample_point = sample_points_[sample_point_index];
			MeshSamplePoint* seed_sample_point = sample_points_[seed_sample_point_index];

			assert(sample_point);
			assert(seed_sample_point);
			assert(label_index_ < sample_point->label_index_confidence_.size());
			assert(label_index_ < seed_sample_point->label_index_confidence_.size());

			if (sample_point->label_index_confidence_[label_index_] >
				seed_sample_point->label_index_confidence_[label_index_])
			{
				seed_sample_point_index = sample_point_index;
			}
		}

		if (seed_sample_point_index < 0) break;
		assert(!is_sample_visited[seed_sample_point_index]);
		assert(sample_points_[seed_sample_point_index]->label_index_confidence_[label_index_]
			>= FLAGS_param_min_sample_point_confidence);


		// 2. Propagate to close points.
		std::deque<SamplePointIndex> queue;
		queue.push_back((SamplePointIndex)seed_sample_point_index);
		is_sample_visited[seed_sample_point_index] = true;

		while (!queue.empty())
		{
			int curr_sample_point_index = queue.front();
			queue.pop_front();
			assert(curr_sample_point_index >= 0 && curr_sample_point_index < num_sample_points());

			sub_cuboid_sample_points.push_back(sample_points_[curr_sample_point_index]);

			MyMesh::Point curr_pos = sample_points_[curr_sample_point_index]->point_;
			q[0] = curr_pos[0]; q[1] = curr_pos[1]; q[2] = curr_pos[2];

			int num_searched_neighbors = _kd_tree->annkFRSearch(q,
				squared_neighbor_distance, num_neighbors, nn_idx);

			for (int i = 0; i < std::min(num_neighbors, num_searched_neighbors); i++)
			{
				unsigned int next_sample_point_index = (int)nn_idx[i];
				if (!is_sample_visited[next_sample_point_index])
				{
					queue.push_back(next_sample_point_index);
					is_sample_visited[next_sample_point_index] = true;
				}
			}
		}

		if (sub_cuboid_sample_points.size() < FLAGS_param_min_num_cuboid_sample_points)
			continue;

		MeshCuboid *cuboid = new MeshCuboid(label_index_);
		cuboid->add_sample_points(sub_cuboid_sample_points);
		cuboid->compute_bbox();


		// Delete too small parts.
		MyMesh::Normal bb_size = cuboid->get_bbox_size();
		if (bb_size[0] < FLAGS_param_min_cuboid_bbox_size * _object_diameter
			&& bb_size[1] < FLAGS_param_min_cuboid_bbox_size * _object_diameter
			&& bb_size[2] < FLAGS_param_min_cuboid_bbox_size * _object_diameter)
		{
			delete cuboid;
			continue;
		}
		else
		{
			_sub_cuboids.push_back(cuboid);
		}
	}

	delete[] is_sample_visited;
	delete[] nn_idx;
	delete[] dd;
	annDeallocPt(q);
}

void MeshCuboid::remove_small_sub_cuboids(std::vector<MeshCuboid *> &_sub_cuboids)
{
	if (_sub_cuboids.empty())
		return;

	// 1. Find the largest part.
	MeshCuboid *seed_cuboid = NULL;
	Real max_bb_volume = 0;

	for (std::vector<MeshCuboid *>::iterator it = _sub_cuboids.begin();
		it != _sub_cuboids.end(); ++it)
	{
		Real bb_volume = (*it)->get_bbox_volume();
		assert(bb_volume >= 0);
		if (bb_volume > max_bb_volume)
		{
			max_bb_volume = bb_volume;
			seed_cuboid = (*it);
		}
	}

	if (!seed_cuboid)
		return;


	// 2. Remove parts that have too small volume.
	for (std::vector<MeshCuboid *>::iterator it = _sub_cuboids.begin();
		it != _sub_cuboids.end();)	// No increment
	{
		Real bb_diag_size = (*it)->get_bbox_diag_length();
		//assert(bb_diag_size <= seed_cuboid->get_bb_diag_size());
		if (bb_diag_size < (1 - FLAGS_param_min_cuboid_bbox_diag_length) * seed_cuboid->get_bbox_diag_length())
		{
			it = _sub_cuboids.erase(it);
		}
		else
		{
			++it;
		}
	}
}

MeshCuboid *MeshCuboid::merge_cuboids(const LabelIndex _label_index,
	const std::vector<MeshCuboid *> _cuboids)
{
	MeshCuboid *merged_cuboid = NULL;

	if (_cuboids.size() == 1)
	{
		merged_cuboid = new MeshCuboid(*(_cuboids.front()));
		merged_cuboid->label_index_ = _label_index;
	}
	else if (_cuboids.size() > 1)
	{
		merged_cuboid = new MeshCuboid(_label_index);

		for (std::vector<MeshCuboid *>::const_iterator it = _cuboids.begin();
			it != _cuboids.end(); ++it)
			merged_cuboid->add_sample_points((*it)->sample_points_);

		merged_cuboid->compute_bbox();
	}

	return merged_cuboid;
}

void MeshCuboid::print_cuboid()const
{
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		std::cout << " - axis (" << axis_index << "): "
			<< bbox_axes_[axis_index][0] << ", "
			<< bbox_axes_[axis_index][1] << ", "
			<< bbox_axes_[axis_index][2] << std::endl;
	}

	std::cout << " - center: "
		<< bbox_center_[0] << ", " << bbox_center_[1] << ", " << bbox_center_[2] << std::endl;

	std::cout << " - size: "
		<< bbox_size_[0] << ", " << bbox_size_[1] << ", " << bbox_size_[2] << std::endl;

	std::cout << " - volume: " << get_bbox_volume() << std::endl;
}

void MeshCuboid::draw_cuboid()const
{
	for (unsigned int face_index = 0; face_index < k_num_faces; ++face_index)
	{
		glBegin(GL_QUADS);
		for (unsigned int i = 0; i < k_num_face_corners; ++i)
		{
			unsigned int corner_index = k_face_corner_indices[face_index][i];
			glVertex3f(
				bbox_corners_[corner_index][0],
				bbox_corners_[corner_index][1],
				bbox_corners_[corner_index][2]);
		}
		glEnd();
	}
}

void MeshCuboid::points_to_cuboid_distances(const Eigen::MatrixXd& _points,
	Eigen::VectorXd &_distances)
{
	// NOTE:
	// Do not use corner points, but use only center, size, and axes.
	// The distance becomes less than zero when the point is inside the cuboid.
	unsigned int num_points = _points.cols();

	Eigen::Vector3d bbox_center_vec;
	for (unsigned int i = 0; i < 3; ++i)
		bbox_center_vec[i] = bbox_center_[i];

	Eigen::MatrixXd centered_points = (_points.colwise() - bbox_center_vec);
	Eigen::MatrixXd axis_distances(3, num_points);

	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		Eigen::Vector3d bbox_axis_vec;
		for (unsigned int i = 0; i < 3; ++i)
			bbox_axis_vec[i] = bbox_axes_[axis_index][i];
		bbox_axis_vec.normalize();

		// Distance from center along each direction.
		axis_distances.row(axis_index) =
			(bbox_axis_vec.transpose() * centered_points).array();

		axis_distances.row(axis_index) = axis_distances.row(axis_index).cwiseAbs();

		// Compute the distance from cuboid surface.
		axis_distances.row(axis_index) = axis_distances.row(axis_index).array()
			- (0.5 * bbox_size_[axis_index]);
	}

	_distances = Eigen::VectorXd(num_points);
	for (unsigned int point_index = 0; point_index < num_points; ++point_index)
	{
		Eigen::Vector3d axis_distance = axis_distances.col(point_index);
		if ((axis_distance.array() <= 0).all())
		{
			_distances[point_index] = axis_distance.maxCoeff();
		}
		else
		{
			axis_distance = (axis_distance.array() < 0).select(0, axis_distance);
			_distances[point_index] = axis_distance.norm();
		}
	}
}

Real MeshCuboid::distance_between_cuboids(
	const MeshCuboid *_cuboid_1, const MeshCuboid *_cuboid_2)
{
	// NOTE: Assume that cuboid surface points exist.
	assert(_cuboid_1);
	assert(_cuboid_2);

	unsigned int num_cuboid_surface_points_1 = _cuboid_1->cuboid_surface_points_.size();
	unsigned int num_cuboid_surface_points_2 = _cuboid_2->cuboid_surface_points_.size();

	// TODO: Pre-compute KD-trees.
	Eigen::MatrixXd cuboid_surface_points_1(3, num_cuboid_surface_points_1);
	Eigen::MatrixXd cuboid_surface_points_2(3, num_cuboid_surface_points_2);

	for (unsigned int point_index = 0; point_index < num_cuboid_surface_points_1; ++point_index)
	{
		for (unsigned int i = 0; i < 3; ++i)
			cuboid_surface_points_1.col(point_index)(i) =
			_cuboid_1->cuboid_surface_points_[point_index]->point_[i];
	}

	for (unsigned int point_index = 0; point_index < num_cuboid_surface_points_2; ++point_index)
	{
		for (unsigned int i = 0; i < 3; ++i)
			cuboid_surface_points_2.col(point_index)(i) =
			_cuboid_2->cuboid_surface_points_[point_index]->point_[i];
	}

	ANNpointArray ann_points_1 = NULL;
	ANNkd_tree* ann_kd_tree_1 = ICP::create_kd_tree(cuboid_surface_points_1, ann_points_1);
	assert(ann_kd_tree_1);

	ANNpointArray ann_points_2 = NULL;
	ANNkd_tree* ann_kd_tree_2 = ICP::create_kd_tree(cuboid_surface_points_2, ann_points_2);
	assert(ann_kd_tree_2);

	// 1 -> 2.
	Eigen::VectorXd distances_12;
	ICP::get_closest_points(ann_kd_tree_2, cuboid_surface_points_1, distances_12);
	assert(distances_12.rows() == num_cuboid_surface_points_1);

	// 2 -> 1.
	Eigen::VectorXd distances_21;
	ICP::get_closest_points(ann_kd_tree_1, cuboid_surface_points_2, distances_21);
	assert(distances_21.rows() == num_cuboid_surface_points_2);

	if (ann_points_1) annDeallocPts(ann_points_1);
	if (ann_points_2) annDeallocPts(ann_points_2);
	delete ann_kd_tree_1;
	delete ann_kd_tree_2;

	Real max_distance = (distances_12.maxCoeff(), distances_21.maxCoeff());
	return max_distance;
}


/*
void MeshCuboid::split_cuboid_recursive(ANNkd_tree* _kd_tree, std::vector<MeshCuboid *> &_sub_cuboids)
{
	for (std::vector<MeshCuboid *>::iterator it = _sub_cuboids.begin();
		it != _sub_cuboids.end(); ++it)
		delete (*it);
	_sub_cuboids.clear();


	std::vector<MeshCuboid *> child_cuboids;
	split_cuboid_once(_kd_tree, child_cuboids);

	_sub_cuboids.insert(_sub_cuboids.end(), child_cuboids.begin(), child_cuboids.end());
	
	//Real total_child_cuboids_volume = 0.0;
	//for (std::vector<MeshCuboid *>::iterator it = child_cuboids.begin(); it != child_cuboids.end(); ++it)
	//	total_child_cuboids_volume += (*it)->get_bbox_volume();

	//if (child_cuboids.empty() || total_child_cuboids_volume == 0.0
	//	|| total_child_cuboids_volume > 0.75 * get_bbox_volume())
	//{
	//	for (std::vector<MeshCuboid *>::iterator it = child_cuboids.begin(); it != child_cuboids.end(); ++it)
	//		delete (*it);

	//	MeshCuboid *cuboid = new MeshCuboid(*this);
	//	_sub_cuboids.push_back(cuboid);
	//	return;
	//}
	//else
	//{
	//	std::vector<MeshCuboid *> sub_tree_cuboids;
	//	for (std::vector<MeshCuboid *>::iterator it = _sub_cuboids.begin(); it != _sub_cuboids.end(); ++it)
	//		(*it)->split_cuboid_recursive(_kd_tree, sub_tree_cuboids);

	//	_sub_cuboids.insert(_sub_cuboids.end(), sub_tree_cuboids.begin(), sub_tree_cuboids.end());
	//}
}

void MeshCuboid::split_cuboid_once(ANNkd_tree* _kd_tree, std::vector<MeshCuboid *> &_sub_cuboids)
{
	for (std::vector<MeshCuboid *>::iterator it = _sub_cuboids.begin();
		it != _sub_cuboids.end(); ++it)
		delete (*it);
	_sub_cuboids.clear();


	if (sample_points_.empty())
		return;

	const unsigned int num_samples = num_sample_points();
	const int dim = 3;
	const int num_neighbors = 16;

	ANNpoint q = annAllocPt(dim);
	ANNidxArray nn_idx = new ANNidx[num_neighbors];
	ANNdistArray dd = new ANNdist[num_neighbors];

	class DistanceComparison
	{
	public:
		bool operator() (const std::pair<int, Real>& lhs, const std::pair<int, Real>& rhs) const
		{
			// Greater-than comparison (find the minimum).
			return (lhs.second > rhs.second);
		}
	};

	typedef std::priority_queue < std::pair<int, Real>,
		std::vector< std::pair<int, Real> >, DistanceComparison > DistanceQueue;
	int seed_sample_point_index[2];


	// The first sample point becomes a random seed point.
	MyMesh::Point random_start_point = sample_points_[0]->point_;
	seed_sample_point_index[0] = 0;
	Real max_distance = 0;

	for (int sample_point_index = 0; sample_point_index < num_samples; ++sample_point_index)
	{
		MyMesh::Point point = sample_points_[sample_point_index]->point_;
		Real distance = (point - random_start_point).norm();
		if (distance > max_distance)
		{
			max_distance = distance;
			seed_sample_point_index[0] = sample_point_index;
		}
	}


	// Find the farthest point from the first seed point.
	seed_sample_point_index[1] = seed_sample_point_index[0];
	max_distance = 0;

	for (int sample_point_index = 0; sample_point_index < num_samples; ++sample_point_index)
	{
		MyMesh::Point point = sample_points_[sample_point_index]->point_;
		Real distance = (point - sample_points_[seed_sample_point_index[0]]->point_).norm();
		if (distance > max_distance)
		{
			max_distance = distance;
			seed_sample_point_index[0] = sample_point_index;
		}
	}

	if (seed_sample_point_index[0] == seed_sample_point_index[1])
		return;


	// Divide into two groups.
	std::vector<MeshSamplePoint *> sub_sample_points[2];
	std::vector<Real> sub_distances[2];
	DistanceQueue sub_queue[2];

	for (unsigned int i = 0; i < 2; ++i)
	{
		MyMesh::Point curr_point = sample_points_[seed_sample_point_index[i]]->point_;
		q[0] = curr_point[0]; q[1] = curr_point[1]; q[2] = curr_point[2];
		_kd_tree->annkSearch(q, num_neighbors, nn_idx, dd);

		for (unsigned int j = 0; j < num_neighbors; j++)
		{
			unsigned int next_sample_point_index = (int)nn_idx[j];
			sub_sample_points[i].push_back(sample_points_[next_sample_point_index]);
		}
	}
	
	for (unsigned int i = 0; i < 2; ++i)
	{
		sub_sample_points[i].reserve(sample_points_.size());
		sub_distances[i].resize(num_samples, std::numeric_limits<Real>::max());
		sub_queue[i].push(std::make_pair(seed_sample_point_index[i], 0));
	}

	while (true)
	{
		for (unsigned int i = 0; i < 2; ++i)
			while (!sub_queue[i].empty() &&
				(sub_queue[i].top().second > sub_distances[0][sub_queue[i].top().first]
				|| sub_queue[i].top().second > sub_distances[1][sub_queue[i].top().first]))
				sub_queue[i].pop();

		unsigned int sub_cuboid_index = 0;

		if (sub_queue[0].empty() && sub_queue[1].empty()) break;
		else if (sub_queue[0].empty())	sub_cuboid_index = 1;
		else if (sub_queue[1].empty())	sub_cuboid_index = 0;
		else if (sub_queue[0].top().second > sub_queue[1].top().second)
			sub_cuboid_index = 1;
		else sub_cuboid_index = 0;


		int curr_sample_point_index = sub_queue[sub_cuboid_index].top().first;
		assert(curr_sample_point_index >= 0 && curr_sample_point_index < num_samples);

		//
		sub_sample_points[sub_cuboid_index].push_back(sample_points_[curr_sample_point_index]);
		//

		sub_distances[sub_cuboid_index][curr_sample_point_index] = -1;
		sub_queue[sub_cuboid_index].pop();

		MyMesh::Point curr_point = sample_points_[curr_sample_point_index]->point_;
		q[0] = curr_point[0]; q[1] = curr_point[1]; q[2] = curr_point[2];

		_kd_tree->annkSearch(q, num_neighbors, nn_idx, dd);

		for (unsigned int i = 0; i < num_neighbors; i++)
		{
			unsigned int next_sample_point_index = (int)nn_idx[i];
			Real distance = (sample_points_[next_sample_point_index]->point_ -
				sample_points_[seed_sample_point_index[sub_cuboid_index]]->point_).norm();
			if (distance < sub_distances[sub_cuboid_index][next_sample_point_index])
			{
				sub_queue[sub_cuboid_index].push(std::make_pair(next_sample_point_index, distance));
				sub_distances[sub_cuboid_index][next_sample_point_index] = distance;
			}
		}
	}

	for (unsigned int i = 0; i < 2; ++i)
	{
		MeshCuboid *cuboid = new MeshCuboid(label_index_);
		cuboid->add_points(sub_sample_points[i]);
		cuboid->compute_bbox();
		_sub_cuboids.push_back(cuboid);
	}


	delete[] nn_idx;
	delete[] dd;
	annDeallocPt(q);
}

void MeshCuboid::align_sub_cuboids(const Real _object_diameter, std::vector<MeshCuboid *> &_sub_cuboids)
{
	if (_sub_cuboids.size() <= 1)
		return;

	// 1. Check whether the center and size of bounding boxes are similar each other.
	unsigned num_cuboids_in_group = _sub_cuboids.size();

	bool similar_bb_center[3];
	bool similar_bb_size[3];
	memset(similar_bb_center, true, 3 * sizeof(bool));
	memset(similar_bb_size, true, 3 * sizeof(bool));


	MyMesh::Normal seed_bb_center = seed_cuboid->get_bbox_center();
	MyMesh::Normal seed_bb_size = seed_cuboid->get_bbox_size();

	MyMesh::Normal avg_bb_center(0.0);
	MyMesh::Normal avg_bb_size(0.0);


	for (std::vector<MeshCuboid *>::iterator it = _sub_cuboids.begin();
		it != _sub_cuboids.end(); ++it)
	{
		MyMesh::Normal bb_center = (*it)->get_bbox_center();
		MyMesh::Normal bb_size = (*it)->get_bbox_size();

		for (unsigned int i = 0; i < 3; i++)
		{
			Real diff_center = abs(bb_center[i] - seed_bb_center[i]);
			if (diff_center > param_sim_abs_attr_tol * _object_diameter)
				similar_bb_center[i] = false;

			Real diff_size = abs(bb_size[i] - seed_bb_size[i]);
			if (diff_size > param_sim_abs_attr_tol * _object_diameter)
				similar_bb_size[i] = false;
		}

		avg_bb_center += bb_center;
		avg_bb_size += bb_size;
	}

	avg_bb_center /= num_cuboids_in_group;
	avg_bb_size /= num_cuboids_in_group;


	// 2. Adjust center and size of bounding boxes.
	for (std::vector<MeshCuboid *>::iterator it = _sub_cuboids.begin();
		it != _sub_cuboids.end(); ++it)
	{
		MyMesh::Normal bb_center = (*it)->get_bbox_center();
		MyMesh::Normal bb_size = (*it)->get_bbox_size();

		for (unsigned int i = 0; i < 3; i++)
		{
			if (similar_bb_center[i]) bb_center[i] = avg_bb_center[i];
			if (similar_bb_size[i]) bb_size[i] = avg_bb_size[i];
		}

		(*it)->set_bbox(bb_center, bb_size);
	}
}
*/
