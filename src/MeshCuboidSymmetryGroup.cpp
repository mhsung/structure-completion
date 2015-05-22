#include "MeshCuboidSymmetryGroup.h"
#include "MeshCuboidParameters.h"
#include "Utilities.h"

#include <bitset>
#include <Eigen/Eigenvalues> 


MeshCuboidSymmetryGroupInfo::MeshCuboidSymmetryGroupInfo()
	: symmetry_type_(ReflectionSymmetryType)
	, aligned_global_axis_index_(0)
{

}

MeshCuboidSymmetryGroupInfo::MeshCuboidSymmetryGroupInfo(MeshCuboidSymmetryGroupType _symmetry_type,
	unsigned int _aligned_global_axis_index)
	: symmetry_type_(_symmetry_type)
	, aligned_global_axis_index_(_aligned_global_axis_index)
{

}

MeshCuboidSymmetryGroupInfo::MeshCuboidSymmetryGroupInfo(const MeshCuboidSymmetryGroupInfo& _other)
{
	symmetry_type_ = _other.symmetry_type_;
	aligned_global_axis_index_ = _other.aligned_global_axis_index_;
	single_label_indices_ = _other.single_label_indices_;
	pair_label_indices_ = _other.pair_label_indices_;
}

MeshCuboidSymmetryGroup::MeshCuboidSymmetryGroup()
	: num_symmetry_orders_(2)
{

}

MeshCuboidSymmetryGroup::MeshCuboidSymmetryGroup(const MeshCuboidSymmetryGroupInfo &_info)
	: info_(_info)
	, num_symmetry_orders_(2)
{

}

MeshCuboidSymmetryGroup::MeshCuboidSymmetryGroup(const MeshCuboidSymmetryGroup &_other)
	:info_(_other.info_)
	, num_symmetry_orders_(_other.num_symmetry_orders_)
{

}

MeshCuboidSymmetryGroup::~MeshCuboidSymmetryGroup()
{

}

unsigned int MeshCuboidSymmetryGroup::num_symmetry_orders() const
{
	assert(num_symmetry_orders_ >= 2);
	return num_symmetry_orders_;
}

Real MeshCuboidSymmetryGroup::get_rotation_angle() const
{
	return M_PI / 180.0 * (360.0 / num_symmetry_orders_);
}

unsigned int MeshCuboidSymmetryGroup::get_aligned_global_axis_index() const
{
	assert(info_.aligned_global_axis_index_ < 3);
	return info_.aligned_global_axis_index_;
}

void MeshCuboidSymmetryGroup::get_single_cuboid_indices(const std::vector<MeshCuboid *>& _cuboids,
	std::vector<unsigned int> &_single_cuboid_indices) const
{
	const unsigned int num_cuboids = _cuboids.size();
	_single_cuboid_indices.clear();

	unsigned int num_single_labels = info_.single_label_indices_.size();
	_single_cuboid_indices.reserve(num_single_labels);

	for (unsigned int i = 0; i < num_single_labels; ++i)
	{
		Label label_index = info_.single_label_indices_[i];
		int cuboid_index = 0;
		for (; cuboid_index < num_cuboids; ++cuboid_index)
		{
			if (_cuboids[cuboid_index]->get_label_index() == label_index) break;
		}

		if (cuboid_index < num_cuboids)
			_single_cuboid_indices.push_back(cuboid_index);
	}
}

void MeshCuboidSymmetryGroup::get_pair_cuboid_indices(const std::vector<MeshCuboid *>& _cuboids,
	std::vector< std::pair<unsigned int, unsigned int> > &_pair_cuboid_indices) const
{
	const unsigned int num_cuboids = _cuboids.size();
	_pair_cuboid_indices.clear();

	unsigned int num_pair_labels = info_.pair_label_indices_.size();
	_pair_cuboid_indices.reserve(num_pair_labels);

	for (unsigned int i = 0; i < num_pair_labels; ++i)
	{
		Label label_index_1 = info_.pair_label_indices_[i].first;
		Label label_index_2 = info_.pair_label_indices_[i].second;

		unsigned int cuboid_index_1 = 0;
		for (; cuboid_index_1 < num_cuboids; ++cuboid_index_1)
			if (_cuboids[cuboid_index_1]->get_label_index() == label_index_1) break;

		unsigned int cuboid_index_2 = 0;
		for (; cuboid_index_2 < num_cuboids; ++cuboid_index_2)
			if (_cuboids[cuboid_index_2]->get_label_index() == label_index_2) break;

		if (cuboid_index_1 < num_cuboids && cuboid_index_2 < num_cuboids)
			_pair_cuboid_indices.push_back(std::make_pair(cuboid_index_1, cuboid_index_2));
	}
}

void MeshCuboidSymmetryGroup::get_symmetric_sample_point_pairs(
	const std::vector<MeshCuboid *> &_cuboids,
	const std::vector<ANNpointArray> &_cuboid_ann_points,
	const std::vector<ANNkd_tree *> &_cuboid_ann_kd_tree,
	const Real _squared_neighbor_distance,
	std::list<WeightedPointPair> &_sample_point_pairs) const
{
	const unsigned int num_cuboids = _cuboids.size();
	assert(_cuboid_ann_points.size() == num_cuboids);
	assert(_cuboid_ann_kd_tree.size() == num_cuboids);

	_sample_point_pairs.clear();


	std::vector<unsigned int> single_cuboid_indices;
	get_single_cuboid_indices(_cuboids, single_cuboid_indices);

	for (std::vector<unsigned int>::const_iterator it = single_cuboid_indices.begin();
		it != single_cuboid_indices.end(); ++it)
	{
		const unsigned int cuboid_index = (*it);

		get_symmetric_sample_point_pairs(
			_cuboids[cuboid_index],
			_cuboid_ann_points[cuboid_index],
			_cuboid_ann_kd_tree[cuboid_index],
			_squared_neighbor_distance, _sample_point_pairs);
	}

	// NOTE:
	// Get inter-cuboid symmetric points only for reflection symmetry.
	if (get_symmetry_type() == ReflectionSymmetryType)
	{
		std::vector< std::pair<unsigned int, unsigned int> > pair_cuboid_indices;
		get_pair_cuboid_indices(_cuboids, pair_cuboid_indices);

		for (std::vector< std::pair<unsigned int, unsigned int> >::const_iterator it = pair_cuboid_indices.begin();
			it != pair_cuboid_indices.end(); ++it)
		{
			const unsigned int cuboid_index_1 = (*it).first;
			const unsigned int cuboid_index_2 = (*it).second;

			// 1 -> 2.
			get_symmetric_sample_point_pairs(
				_cuboids[cuboid_index_1],
				_cuboid_ann_points[cuboid_index_2],
				_cuboid_ann_kd_tree[cuboid_index_2],
				_squared_neighbor_distance, _sample_point_pairs);

			// 2 -> 1.
			get_symmetric_sample_point_pairs(
				_cuboids[cuboid_index_2],
				_cuboid_ann_points[cuboid_index_1],
				_cuboid_ann_kd_tree[cuboid_index_1],
				_squared_neighbor_distance, _sample_point_pairs);
		}
	}
}

void MeshCuboidSymmetryGroup::get_symmetric_sample_point_pairs(
	const MeshCuboid *_cuboid_1,
	const ANNpointArray &_cuboid_ann_points_2,
	ANNkd_tree *_cuboid_ann_kd_tree_2,
	const Real _squared_neighbor_distance,
	std::list<WeightedPointPair> &_sample_point_pairs) const
{
	if (!_cuboid_1) return;
	std::vector<MyMesh::Point> cuboid_1_sample_points;
	_cuboid_1->get_sample_points(cuboid_1_sample_points);

	get_symmetric_sample_point_pairs(
		cuboid_1_sample_points,
		_cuboid_ann_points_2,
		_cuboid_ann_kd_tree_2,
		_squared_neighbor_distance,
		_sample_point_pairs);
}

void MeshCuboidSymmetryGroup::get_symmetric_sample_point_pairs(
	const std::vector<MyMesh::Point> &_cuboid_1_sample_points,
	const ANNpointArray &_cuboid_ann_points_2,
	ANNkd_tree *_cuboid_ann_kd_tree_2,
	const Real _squared_neighbor_distance,
	std::list<WeightedPointPair> &_sample_point_pairs) const
{
	if (_cuboid_1_sample_points.size() == 0 || !_cuboid_ann_points_2 || !_cuboid_ann_kd_tree_2)
		return;

	const int num_points_1 = _cuboid_1_sample_points.size();
	const int num_points_2 = _cuboid_ann_kd_tree_2->nPoints();

	//
	ANNpoint q = annAllocPt(3);
	ANNidxArray nn_idx = new ANNidx[1];
	ANNdistArray dd = new ANNdist[1];
	//

	for (int point_index_1 = 0; point_index_1 < num_points_1; ++point_index_1)
	{
		MyMesh::Point point_1 = _cuboid_1_sample_points[point_index_1];

		for (int symmetry_order = 1; symmetry_order < num_symmetry_orders_; ++symmetry_order)
		{
			MyMesh::Point symmetric_point_1 = get_symmetric_point(point_1, symmetry_order);

			for (unsigned int i = 0; i < 3; ++i)
				q[i] = symmetric_point_1[i];

			int num_searched_neighbors = _cuboid_ann_kd_tree_2->annkFRSearch(
				q, _squared_neighbor_distance, 1, nn_idx, dd);
			if (num_searched_neighbors > 0)
			{
				int point_index_2 = (int)nn_idx[0];
				assert(point_index_2 < num_points_2);

				MyMesh::Point point_2;
				for (unsigned int i = 0; i < 3; ++i)
					point_2[i] = _cuboid_ann_points_2[point_index_2][i];

				//
				// Debug.
				Real test_distance = (symmetric_point_1 - point_2).norm();
				CHECK_NUMERICAL_ERROR(__FUNCTION__, test_distance, std::sqrt(dd[0]));
				double distance = (std::sqrt(_squared_neighbor_distance) - std::sqrt(dd[0]));
				assert(distance >= 0);
				//

				Real weight = distance * distance;
				Real angle = symmetry_order * get_rotation_angle();
				WeightedPointPair point_pair(point_1, point_2, weight, angle);
				_sample_point_pairs.push_back(point_pair);
			}
		}
	}

	//
	annDeallocPt(q);
	delete[] nn_idx;
	delete[] dd;
	//
}

MeshCuboidReflectionSymmetryGroup* MeshCuboidReflectionSymmetryGroup::constructor(
	const MeshCuboidSymmetryGroupInfo &_info,
	const std::vector<MeshCuboid *>& _cuboids)
{
	MeshCuboidReflectionSymmetryGroup *group = new MeshCuboidReflectionSymmetryGroup(_info);
	bool ret = group->compute_symmetry_axis(_cuboids);
	if (!ret)
	{
		delete group;
		group = NULL;
	}
	return group;
}

MeshCuboidReflectionSymmetryGroup::MeshCuboidReflectionSymmetryGroup(
	const MyMesh::Normal _n, const double _t)
	: MeshCuboidSymmetryGroup()
	, n_(_n)
	, t_(_t)
{

}

MeshCuboidReflectionSymmetryGroup::MeshCuboidReflectionSymmetryGroup(
	const MeshCuboidSymmetryGroupInfo &_info)
	: MeshCuboidSymmetryGroup(_info)
	, n_(MyMesh::Normal(1, 0, 0))
	, t_(0)
{
	assert(_info.symmetry_type_ == ReflectionSymmetryType);
}

MeshCuboidReflectionSymmetryGroup::MeshCuboidReflectionSymmetryGroup(
	const MeshCuboidReflectionSymmetryGroup &_other)
	: MeshCuboidSymmetryGroup(_other)
	, n_(_other.n_)
	, t_(_other.t_)
{
	n_.normalize();
	assert(info_.symmetry_type_ == ReflectionSymmetryType);
}

MeshCuboidReflectionSymmetryGroup::~MeshCuboidReflectionSymmetryGroup()
{

}

unsigned int MeshCuboidReflectionSymmetryGroup::num_axis_parameters()
{
	return (3 + 1);
}

MeshCuboidSymmetryGroupType MeshCuboidReflectionSymmetryGroup::get_symmetry_type() const
{
	assert(info_.symmetry_type_ == ReflectionSymmetryType);
	return info_.symmetry_type_;
}

unsigned int MeshCuboidReflectionSymmetryGroup::num_symmetry_orders() const
{
	// NOTE:
	// Reflection symmetry is order of 2.
	assert(num_symmetry_orders_ == 2);
	return num_symmetry_orders_;
}

void MeshCuboidReflectionSymmetryGroup::set_num_symmetry_order(unsigned int _symmetry_order)
{
	// NOTE:
	// Reflection symmetry is order of 2.
	assert(false);
}

void MeshCuboidReflectionSymmetryGroup::get_reflection_plane(MyMesh::Normal &_n, double &_t) const
{
	_n = n_;
	_t = t_;
}

void MeshCuboidReflectionSymmetryGroup::set_reflection_plane(const MyMesh::Normal &_n, const double &_t)
{
	n_ = _n;
	t_ = _t;
	n_.normalize();
}

bool MeshCuboidReflectionSymmetryGroup::compute_symmetry_axis(const std::vector<MeshCuboid *>& _cuboids)
{
	std::vector< std::pair<unsigned int, unsigned int> > pair_cuboid_indices;
	get_pair_cuboid_indices(_cuboids, pair_cuboid_indices);

	if (pair_cuboid_indices.empty())
	{
		std::vector<unsigned int> single_cuboid_indices;
		get_single_cuboid_indices(_cuboids, single_cuboid_indices);
		if (single_cuboid_indices.empty())
			return false;

		const unsigned int num_cuboids = _cuboids.size();

		std::list<MyMesh::Normal> aligned_normals;
		std::list<MyMesh::Point> aligned_centers;

		for (std::vector<unsigned int>::const_iterator it = single_cuboid_indices.begin();
			it != single_cuboid_indices.end(); ++it)
		{
			int cuboid_index = (*it);
			assert(cuboid_index < num_cuboids);
			const MeshCuboid *cuboid = _cuboids[cuboid_index];

			aligned_normals.push_back(cuboid->get_bbox_axis(get_aligned_global_axis_index()));
			aligned_centers.push_back(cuboid->get_bbox_center());
		}

		assert(!aligned_normals.empty());
		assert(!aligned_centers.empty());

		Eigen::Matrix3d A = Eigen::Matrix3d::Zero();

		for (std::list<MyMesh::Normal>::const_iterator it = aligned_normals.begin();
			it != aligned_normals.end(); ++it)
		{
			Eigen::Vector3d n_vec;
			for (int i = 0; i < 3; ++i) n_vec[i] = (*it)[i];
			A += (n_vec * n_vec.transpose());
		}

		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(A);
		Eigen::Vector3d n = es.eigenvectors().col(3 - 1);
		for (unsigned int i = 0; i < 3; ++i) n_[i] = n[i];
		n_.normalize();

		t_ = 0;
		for (std::list<MyMesh::Point>::const_iterator it = aligned_centers.begin();
			it != aligned_centers.end(); ++it)
			t_ += dot(n_, (*it));
		t_ /= aligned_centers.size();
	}
	else
	{
		const unsigned int num_cuboids = _cuboids.size();

		std::list< std::pair < MyMesh::Point, MyMesh::Point > > point_pairs;
		for (std::vector < std::pair <unsigned int, unsigned int> >::const_iterator it = pair_cuboid_indices.begin();
			it != pair_cuboid_indices.end(); ++it)
		{
			int cuboid_index_1 = (*it).first;
			int cuboid_index_2 = (*it).second;
			assert(cuboid_index_1 < num_cuboids);
			assert(cuboid_index_2 < num_cuboids);

			const MeshCuboid *cuboid_1 = _cuboids[cuboid_index_1];
			const MeshCuboid *cuboid_2 = _cuboids[cuboid_index_2];

			add_symmety_cuboid_corner_points(cuboid_1, cuboid_2, _cuboids,
				info_.aligned_global_axis_index_, point_pairs);
		}

		assert(!point_pairs.empty());
		compute_reflection_plane(point_pairs, n_, t_);
	}

	return true;
}

void MeshCuboidReflectionSymmetryGroup::add_symmety_cuboid_corner_points(
	const MeshCuboid *cuboid_1, const MeshCuboid *cuboid_2,
	const std::vector<MeshCuboid *>& _cuboids,
	const unsigned int _reflection_axis_index,
	std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs)
{
	assert(cuboid_1);
	assert(cuboid_2);

	const unsigned int dimension = 3;
	assert(_reflection_axis_index < dimension);

	// NOTE:
	// Implemented only for 3 dimension.
	assert(dimension == 3);

	for (unsigned int corner_index_1 = 0; corner_index_1 < MeshCuboid::k_num_corners; ++corner_index_1)
	{
		std::bitset<3> bits(corner_index_1);
		bits[_reflection_axis_index].flip();
		unsigned int corner_index_2 = bits.to_ulong();

		MyMesh::Point corner_point_1 = cuboid_1->get_bbox_corner(corner_index_1);
		MyMesh::Point corner_point_2 = cuboid_2->get_bbox_corner(corner_index_2);
		
		_point_pairs.push_back(std::make_pair(corner_point_1, corner_point_2));
	}
}

void MeshCuboidReflectionSymmetryGroup::compute_reflection_plane(
	const std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs,
	MyMesh::Normal &_n, double &_t)
{
	// Eq (1):
	// min {(I - nn^T)(x - y)}^2. Let d = (x - y).
	// Since (I - nn^T)^2 = (I - nn^T),
	// => min d^T(I - nn^T)d = d^Td - d^Tnn^Td = d^Td - n^Tdd^Tn.
	// => max n^Tdd^Tn.

	// Eq (2):
	// min {n^T(x + y) - 2t}^2.
	// If n is estimated, t = 0.5*n^T(x + y).

	assert(!_point_pairs.empty());

	Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
	Eigen::Vector3d b = Eigen::Vector3d::Zero();

	for (std::list< std::pair < MyMesh::Point, MyMesh::Point > >::const_iterator it = _point_pairs.begin();
		it != _point_pairs.end(); ++it)
	{
		Eigen::Vector3d sum_p, diff_p;
		for (int i = 0; i < 3; ++i)
		{
			sum_p[i] = (*it).first[i] + (*it).second[i];
			diff_p[i] = (*it).first[i] - (*it).second[i];
		}

		A += (diff_p * diff_p.transpose());
		b += (0.5 * sum_p);
	}

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(A);
	Eigen::Vector3d n = es.eigenvectors().col(3 - 1);
	for (unsigned int i = 0; i < 3; ++i) _n[i] = n[i];

	_n.normalize();
	_t = (n.transpose() * b);
	_t /= _point_pairs.size();
}

void MeshCuboidReflectionSymmetryGroup::get_reflection_plane_corners(
	const MyMesh::Point &_point, const Real _size,
	std::array<MyMesh::Point, 4>& _corners) const
{
	// Find the closest point on the plane to the given '_point'.
	// p_proj = (p - a) - (n'(p-a))n + a  ('a' is a point on the plane, n'a = t)
	// p_proj = p - (n'p)n + tn

	MyMesh::Point center = _point - dot(n_, _point) * n_ + t_ * n_;

	// Find random two axes of plane.
	MyMesh::Normal axes[2];

	axes[0] = MyMesh::Normal(0.0);
	for (unsigned int i = 0; i < 3; ++i)
	{
		MyMesh::Normal temp(0.0);
		temp[i] = 1.0;
		axes[0] = cross(n_, temp);
		if (axes[0].norm() > 0.1)
		{
			axes[0].normalize();
			break;
		}
	}
	CHECK_NUMERICAL_ERROR(__FUNCTION__, axes[0].norm(), 1.0);
	axes[1] = cross(n_, axes[0]);
	axes[1].normalized();

	_corners[0] = center - 0.5 * _size * axes[0] - 0.5 * _size * axes[1];
	_corners[1] = center + 0.5 * _size * axes[0] - 0.5 * _size * axes[1];
	_corners[2] = center + 0.5 * _size * axes[0] + 0.5 * _size * axes[1];
	_corners[3] = center - 0.5 * _size * axes[0] + 0.5 * _size * axes[1];
}

MyMesh::Point MeshCuboidReflectionSymmetryGroup::get_symmetric_point(
	const MyMesh::Point& _point, unsigned int _symmetry_order) const
{
	// NOTE:
	// Symmetry order of 1 means reflection.
	assert(_symmetry_order == 1);
	return get_symmetric_point(_point);
}

MyMesh::Normal MeshCuboidReflectionSymmetryGroup::get_symmetric_normal(
	const MyMesh::Normal& _normal, unsigned int _symmetry_order) const
{
	// NOTE:
	// Symmetry order of 1 means reflection.
	assert(_symmetry_order == 1);
	return get_symmetric_normal(_normal);
}

MyMesh::Point MeshCuboidReflectionSymmetryGroup::get_symmetric_point(const MyMesh::Point& _point) const
{
	// Assume that 'n_' is normalized.
	MyMesh::Normal plane_to_point = n_ * (dot(n_, _point) - t_);
	MyMesh::Point symmetric_point = _point - (plane_to_point * 2);
	CHECK_NUMERICAL_ERROR(__FUNCTION__, t_ - dot(n_, 0.5 * (_point + symmetric_point)));
	return symmetric_point;
}

MyMesh::Normal MeshCuboidReflectionSymmetryGroup::get_symmetric_normal(const MyMesh::Normal& _normal) const
{
	// Assume that 'n_' and 'normal_' are normalized.
	MyMesh::Normal n_direction_component = n_ * (dot(n_, _normal));
	MyMesh::Normal symmetric_normal = _normal - (n_direction_component * 2);
	CHECK_NUMERICAL_ERROR(__FUNCTION__, dot(n_, 0.5 * (_normal + symmetric_normal)));
	return symmetric_normal;
}

MeshCuboidRotationSymmetryGroup* MeshCuboidRotationSymmetryGroup::constructor(
	const MeshCuboidSymmetryGroupInfo &_info,
	const std::vector<MeshCuboid *>& _cuboids)
{
	MeshCuboidRotationSymmetryGroup *group = new MeshCuboidRotationSymmetryGroup(_info);
	bool ret = (group->compute_symmetry_axis(_cuboids) && group->compute_rotation_angle(_cuboids));
	if (!ret)
	{
		delete group;
		group = NULL;
	}
	return group;
}

MeshCuboidRotationSymmetryGroup::MeshCuboidRotationSymmetryGroup(
	const MeshCuboidSymmetryGroupInfo &_info)
	: MeshCuboidSymmetryGroup(_info)
	, n_(MyMesh::Normal(1, 0, 0))
	, t_(MyMesh::Point(0, 0, 0))
{
	assert(_info.symmetry_type_ == RotationSymmetryType);
}

MeshCuboidRotationSymmetryGroup::MeshCuboidRotationSymmetryGroup(
	const MeshCuboidRotationSymmetryGroup &_other)
	: MeshCuboidSymmetryGroup(_other)
	, n_(_other.n_)
	, t_(_other.t_)
{
	n_.normalize();
	assert(info_.symmetry_type_ == RotationSymmetryType);
}

MeshCuboidRotationSymmetryGroup::~MeshCuboidRotationSymmetryGroup()
{

}

unsigned int MeshCuboidRotationSymmetryGroup::num_axis_parameters()
{
	return (3 + 3);
}

MeshCuboidSymmetryGroupType MeshCuboidRotationSymmetryGroup::get_symmetry_type() const
{
	assert(info_.symmetry_type_ == RotationSymmetryType);
	return info_.symmetry_type_;
}

bool MeshCuboidRotationSymmetryGroup::compute_symmetry_axis(const std::vector<MeshCuboid *>& _cuboids)
{
	std::vector<unsigned int> single_cuboid_indices;
	get_single_cuboid_indices(_cuboids, single_cuboid_indices);
	if (single_cuboid_indices.empty())
		return false;

	const unsigned int num_cuboids = _cuboids.size();

	std::list<MyMesh::Normal> aligned_normals;
	std::list<MyMesh::Point> aligned_centers;

	for (std::vector<unsigned int>::const_iterator it = single_cuboid_indices.begin();
		it != single_cuboid_indices.end(); ++it)
	{
		int cuboid_index = (*it);
		assert(cuboid_index < num_cuboids);
		const MeshCuboid *cuboid = _cuboids[cuboid_index];

		aligned_normals.push_back(cuboid->get_bbox_axis(get_aligned_global_axis_index()));
		aligned_centers.push_back(cuboid->get_bbox_center());
	}

	assert(!aligned_normals.empty());
	assert(!aligned_centers.empty());

	Eigen::Matrix3d A = Eigen::Matrix3d::Zero();

	for (std::list<MyMesh::Normal>::const_iterator it = aligned_normals.begin();
		it != aligned_normals.end(); ++it)
	{
		Eigen::Vector3d n_vec;
		for (int i = 0; i < 3; ++i) n_vec[i] = (*it)[i];
		A += (n_vec * n_vec.transpose());
	}

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(A);
	Eigen::Vector3d n = es.eigenvectors().col(3 - 1);
	for (unsigned int i = 0; i < 3; ++i) n_[i] = n[i];
	n_.normalize();

	t_ = MyMesh::Point(0.0);
	for (std::list<MyMesh::Point>::const_iterator it = aligned_centers.begin();
		it != aligned_centers.end(); ++it)
		t_ += (*it);
	t_ /= aligned_centers.size();

	return true;
}

MyMesh::Point MeshCuboidRotationSymmetryGroup::get_symmetric_point(
	const MyMesh::Point& _point, unsigned int _symmetry_order) const
{
	assert(_symmetry_order > 0);
	assert(_symmetry_order < num_symmetry_orders_);

	Eigen::Vector3d axis_vec, point_vec;
	for (unsigned int i = 0; i < 3; ++i)
	{
		axis_vec[i] = n_[i];
		point_vec[i] = _point[i] - t_[i];
	}

	double angle = get_rotation_angle() * _symmetry_order;

	Eigen::AngleAxisd axis_rotation(angle, axis_vec);
	Eigen::Vector3d rot_point_vec = axis_rotation.toRotationMatrix() * point_vec;

	MyMesh::Point rot_point;
	for (unsigned int i = 0; i < 3; ++i)
		rot_point[i] = rot_point_vec[i] + t_[i];

	return rot_point;
}

MyMesh::Normal MeshCuboidRotationSymmetryGroup::get_symmetric_normal(
	const MyMesh::Normal& _normal, unsigned int _symmetry_order) const
{
	assert(_symmetry_order > 0);
	assert(_symmetry_order < num_symmetry_orders_);

	Eigen::Vector3d axis_vec, normal_vec;
	for (unsigned int i = 0; i < 3; ++i)
	{
		axis_vec[i] = n_[i];
		normal_vec[i] = _normal[i];
	}

	double angle = get_rotation_angle() * (_symmetry_order + 1);

	Eigen::AngleAxisd axis_rotation(angle, axis_vec);
	Eigen::Vector3d rot_normal_vec = axis_rotation.toRotationMatrix() * normal_vec;

	MyMesh::Point ret;
	for (unsigned int i = 0; i < 3; ++i)
	{
		ret[i] = rot_normal_vec[i];
	}

	return ret;
}

void MeshCuboidRotationSymmetryGroup::get_rotation_axis(MyMesh::Normal &_n, MyMesh::Point &_t) const
{
	_n = n_;
	_t = t_;
}

void MeshCuboidRotationSymmetryGroup::set_rotation_axis(const MyMesh::Normal &_n, const MyMesh::Point &_t)
{
	n_ = _n;
	t_ = _t;
	n_.normalize();
}

void MeshCuboidRotationSymmetryGroup::get_rotation_axis_corners(
	const MyMesh::Point &_point, const Real _size, std::array<MyMesh::Point, 2>& _corners) const
{
	// Find the closest point on the axis to the given '_point'.
	// min || n*x + t - p ||.
	// x = -(n^T (t - p)).

	double x = -dot(n_, t_ - _point);
	MyMesh::Point center = n_ * x + t_;

	_corners[0] = center - n_ * 0.5 * _size;
	_corners[1] = center + n_ * 0.5 * _size;
}

void MeshCuboidRotationSymmetryGroup::add_symmety_cuboid_corner_points(
	const MeshCuboid *cuboid_1, const MeshCuboid *cuboid_2,
	const std::vector<MeshCuboid *>& _cuboids,
	const unsigned int _reflection_axis_index,
	std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs)
{

}

void MeshCuboidRotationSymmetryGroup::compute_rotation_plane(
	const std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs,
	MyMesh::Normal &_n, double &_t)
{

}

bool MeshCuboidRotationSymmetryGroup::compute_rotation_angle(
	const std::vector<MeshCuboid *> &_cuboids)
{
	const Real squared_neighbor_distance = FLAGS_param_sparse_neighbor_distance
		* FLAGS_param_sparse_neighbor_distance;
	const unsigned int num_cuboids = _cuboids.size();

	//
	std::vector<ANNpointArray> cuboid_ann_points(num_cuboids);
	std::vector<ANNkd_tree *> cuboid_ann_kd_tree(num_cuboids);

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		cuboid_ann_kd_tree[cuboid_index] = NULL;
		cuboid_ann_points[cuboid_index] = NULL;

		MeshCuboid *cuboid = _cuboids[cuboid_index];
		unsigned int num_cuboid_sample_points = cuboid->num_sample_points();
		if (num_cuboid_sample_points == 0)
			continue;

		Eigen::MatrixXd cuboid_sample_points(3, num_cuboid_sample_points);

		for (unsigned int point_index = 0; point_index < num_cuboid_sample_points; ++point_index)
		{
			for (unsigned int i = 0; i < 3; ++i)
				cuboid_sample_points.col(point_index)(i) =
				cuboid->get_sample_point(point_index)->point_[i];
		}

		cuboid_ann_kd_tree[cuboid_index] = ICP::create_kd_tree(cuboid_sample_points,
			cuboid_ann_points[cuboid_index]);
		assert(cuboid_ann_points[cuboid_index]);
		assert(cuboid_ann_kd_tree[cuboid_index]);

		assert(cuboid_ann_points.size() == num_cuboids);
		assert(cuboid_ann_kd_tree.size() == num_cuboids);
	}
	//

	std::vector<unsigned int> single_cuboid_indices;
	get_single_cuboid_indices(_cuboids, single_cuboid_indices);
	if (single_cuboid_indices.empty())
		return false;

	// FIXME.
	const unsigned int min_num_symmetry_orders = 3;
	const unsigned int max_num_symmetry_orders = 6;
	const unsigned int default_num_symmetry_orders = 5;

	Real max_num_point_pairs = 0.0;
	unsigned int best_num_symmetry_orders = default_num_symmetry_orders;

	for (unsigned int num_symmetry_orders = min_num_symmetry_orders; num_symmetry_orders <= max_num_symmetry_orders;
		++num_symmetry_orders)
	{
		assert(num_symmetry_orders > 2);

		// Temporary set symmetry order.
		num_symmetry_orders_ = num_symmetry_orders;
		std::list<WeightedPointPair> sample_point_pairs;
		get_symmetric_sample_point_pairs(_cuboids, cuboid_ann_points, cuboid_ann_kd_tree,
			squared_neighbor_distance, sample_point_pairs);

		Real num_point_pairs = sample_point_pairs.size();

		// Divide by the number of rotations.
		num_point_pairs /= (num_symmetry_orders - 1);

		std::cout << "[" << num_symmetry_orders << "]: " << num_point_pairs << std::endl;

		if (num_point_pairs > max_num_point_pairs)
		{
			max_num_point_pairs = num_point_pairs;
			best_num_symmetry_orders = num_symmetry_orders;
		}
	}

	//
	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		if (cuboid_ann_points[cuboid_index]) annDeallocPts(cuboid_ann_points[cuboid_index]);
		if (cuboid_ann_kd_tree[cuboid_index]) delete cuboid_ann_kd_tree[cuboid_index];
	}
	cuboid_ann_points.clear();
	cuboid_ann_kd_tree.clear();
	//

	num_symmetry_orders_ = best_num_symmetry_orders;
	std::cout << "num_symmetry_orders = " << num_symmetry_orders_ << std::endl;

	return true;
}
