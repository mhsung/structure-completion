#include "MeshCuboidSymmetryGroup.h"
#include "Utilities.h"

#include <bitset>
#include <Eigen/Eigenvalues> 


MeshCuboidSymmetryGroupInfo::MeshCuboidSymmetryGroupInfo(unsigned int _reflection_axis_index)
	: reflection_axis_index_(_reflection_axis_index)
{

}

MeshCuboidSymmetryGroupInfo::MeshCuboidSymmetryGroupInfo(const MeshCuboidSymmetryGroupInfo& _other)
{
	reflection_axis_index_ = _other.reflection_axis_index_;
	single_label_indices_ = _other.single_label_indices_;
	pair_label_indices_ = _other.pair_label_indices_;
}

MeshCuboidSymmetryGroup* MeshCuboidSymmetryGroup::constructor(
	const MeshCuboidSymmetryGroupInfo &_info,
	const std::vector<MeshCuboid *>& _cuboids)
{
	MeshCuboidSymmetryGroup *group = new MeshCuboidSymmetryGroup(_info);
	bool ret = group->compute_symmetry_axis(_cuboids);
	if (!ret)
	{
		delete group;
		group = NULL;
	}
	return group;
}

MeshCuboidSymmetryGroup::MeshCuboidSymmetryGroup(const MeshCuboidSymmetryGroupInfo &_info)
	: info_(_info)
	, n_(MyMesh::Normal(1, 0, 0))
	, t_(0)
{
	
}

MeshCuboidSymmetryGroup::~MeshCuboidSymmetryGroup()
{

}

bool MeshCuboidSymmetryGroup::compute_symmetry_axis(const std::vector<MeshCuboid *>& _cuboids)
{
	std::vector< std::pair<unsigned int, unsigned int> > pair_cuboid_indices;
	get_pair_cuboid_indices(_cuboids, pair_cuboid_indices);
	if (pair_cuboid_indices.empty())
		return false;

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
			info_.reflection_axis_index_, point_pairs);
	}

	assert(!point_pairs.empty());
	compute_reflection_plane(point_pairs, n_, t_);
	return true;
}

void MeshCuboidSymmetryGroup::get_reflection_plane(MyMesh::Normal &_n, double &_t) const
{
	_n = n_;
	_t = t_;
}

void MeshCuboidSymmetryGroup::set_reflection_plane(const MyMesh::Normal &_n, const double &_t)
{
	n_ = _n;
	t_ = _t;
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
			if (_cuboids[cuboid_index]->get_label_index() == label_index) break;

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

void MeshCuboidSymmetryGroup::add_symmety_cuboid_corner_points(
	const MeshCuboid *cuboid_1, const MeshCuboid *cuboid_2,
	const std::vector<MeshCuboid *>& _cuboids,
	const unsigned int _reflection_axis_index,
	std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs)
{
	assert(cuboid_1);
	assert(cuboid_2);

	const unsigned int dimension = 3;
	assert(_reflection_axis_index < dimension);

	// NOTICE:
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

void MeshCuboidSymmetryGroup::compute_reflection_plane(
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

	_t = (n.transpose() * b);
	_t /= _point_pairs.size();
}

void MeshCuboidSymmetryGroup::get_reflection_plane_corners(MyMesh::Point &_point, double _size,
	std::array<MyMesh::Point, 4>& _corners)
{
	// Find the closest point on the plane to the given '_point'.
	// p_proj = (p - a) - (n'(p-a))n + a  ('a' is a point on the plane, n'a = t)
	// p_proj = p - (n'p)n + tn

	n_.normalize();

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

MyMesh::Point MeshCuboidSymmetryGroup::get_symmetric_point(const MyMesh::Point& _point) const
{
	MyMesh::Normal plane_to_point = n_ * (dot(n_, _point) - t_);
	MyMesh::Point symmetric_point = _point - (plane_to_point * 2);
	CHECK_NUMERICAL_ERROR(__FUNCTION__, t_ - dot(n_, 0.5 * (_point + symmetric_point)));
	return symmetric_point;
}

void MeshCuboidSymmetryGroup::get_symmetric_sample_point_pairs(
	const std::vector<MeshCuboid *> &_cuboids,
	const std::vector<ANNpointArray> &_cuboid_ann_points,
	const std::vector<ANNkd_tree *> &_cuboid_ann_kd_tree,
	const Real _neighbor_distance,
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
			_neighbor_distance, _sample_point_pairs);
	}


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
			_neighbor_distance, _sample_point_pairs);

		// 2 -> 1.
		get_symmetric_sample_point_pairs(
			_cuboids[cuboid_index_2],
			_cuboid_ann_points[cuboid_index_1],
			_cuboid_ann_kd_tree[cuboid_index_1],
			_neighbor_distance, _sample_point_pairs);
	}
}

void MeshCuboidSymmetryGroup::get_symmetric_sample_point_pairs(
	const MeshCuboid *_cuboid_1,
	const ANNpointArray &_cuboid_ann_points_2,
	ANNkd_tree *_cuboid_ann_kd_tree_2,
	const Real neighbor_distance,
	std::list<WeightedPointPair> &_sample_point_pairs) const
{
	if (!_cuboid_1 || !_cuboid_ann_points_2 || !_cuboid_ann_kd_tree_2)
		return;

	const int num_points_1 = _cuboid_1->num_sample_points();
	const int num_points_2 = _cuboid_ann_kd_tree_2->nPoints();

	//
	ANNpoint q = annAllocPt(3);
	ANNidxArray nn_idx = new ANNidx[1];
	ANNdistArray dd = new ANNdist[1];
	//

	for (int point_index_1 = 0; point_index_1 < num_points_1; ++point_index_1)
	{
		MyMesh::Point point_1 = _cuboid_1->get_sample_point(point_index_1)->point_;
		MyMesh::Point symmetric_point_1 = get_symmetric_point(point_1);

		for (unsigned int i = 0; i < 3; ++i)
			q[i] = symmetric_point_1[i];

		int num_searched_neighbors = _cuboid_ann_kd_tree_2->annkFRSearch(q, neighbor_distance, 1, nn_idx, dd);
		if (num_searched_neighbors > 0)
		{
			int point_index_2 = (int)nn_idx[0];
			assert(point_index_2 < num_points_2);

			MyMesh::Point point_2;
			for (unsigned int i = 0; i < 3; ++i)
				point_2[i] = _cuboid_ann_points_2[point_index_2][i];

			//
			double distance = (neighbor_distance - dd[0]);
			assert(distance >= 0);
			//

			double weight = distance * distance;

			WeightedPointPair point_pair(weight, point_1, point_2);
			_sample_point_pairs.push_back(point_pair);
		}
	}

	//
	annDeallocPt(q);
	delete[] nn_idx;
	delete[] dd;
	//
}