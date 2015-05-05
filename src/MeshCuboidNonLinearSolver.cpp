#include "MeshCuboidNonLinearSolver.h"

#include <bitset>

#include <Eigen/Eigenvalues> 


MeshCuboidNonLinearSolver::MeshCuboidNonLinearSolver(
	const std::vector<MeshCuboid *>& _cuboids,
	const std::vector<MeshCuboidReflectionSymmetryGroup *>& _reflection_symmetry_groups,
	const std::vector<MeshCuboidRotationSymmetryGroup *>& _rotation_symmetry_groups,
	const Real _neighbor_distance,
	const Real _symmetry_energy_term_weight)
	: cuboids_(_cuboids)
	, reflection_symmetry_groups_(_reflection_symmetry_groups)
	, rotation_symmetry_groups_(_rotation_symmetry_groups)
	, neighbor_distance_(_neighbor_distance)
	, symmetry_energy_term_weight_(_symmetry_energy_term_weight)
	, num_cuboid_corner_variables_(MeshCuboidAttributes::k_num_attributes)
	, num_cuboid_axis_variables_(3 * 3)
	, num_reflection_symmetry_group_variables_(MeshCuboidReflectionSymmetryGroup::num_axis_parameters())
	, num_rotation_symmetry_group_variables_(MeshCuboidRotationSymmetryGroup::num_axis_parameters())
	, constant_term_(0)
{
	num_cuboids_ = _cuboids.size();
	//assert(num_cuboids_ > 0);

	num_reflection_symmetry_groups_ = _reflection_symmetry_groups.size();
	//assert(num_reflection_symmetry_groups_ > 0);
	num_rotation_symmetry_groups_ = _rotation_symmetry_groups.size();
	//assert(num_rotation_symmetry_groups_ > 0);

	cuboid_corner_variable_start_index_ = 0;
	cuboid_axis_variable_start_index_ = cuboid_corner_variable_start_index_
		+ num_total_cuboid_corner_variables();
	reflection_symmetry_group_variable_start_index_ = cuboid_axis_variable_start_index_
		+ num_total_cuboid_axis_variables();
	rotation_symmetry_group_variable_start_index_ = reflection_symmetry_group_variable_start_index_
		+ num_total_reflection_symmetry_group_variables();
}

MeshCuboidNonLinearSolver::~MeshCuboidNonLinearSolver()
{

}

unsigned int MeshCuboidNonLinearSolver::num_total_cuboid_corner_variables() const
{
	return num_cuboid_corner_variables_ * num_cuboids_;
}

unsigned int MeshCuboidNonLinearSolver::num_total_cuboid_axis_variables() const
{
	return num_cuboid_axis_variables_ * num_cuboids_;
}

unsigned int MeshCuboidNonLinearSolver::num_total_reflection_symmetry_group_variables() const
{
	return num_reflection_symmetry_groups_ * num_reflection_symmetry_group_variables_;
}

unsigned int MeshCuboidNonLinearSolver::num_total_rotation_symmetry_group_variables() const
{
	return num_rotation_symmetry_groups_ * num_rotation_symmetry_group_variables_;
}

unsigned int MeshCuboidNonLinearSolver::num_total_variables() const
{
	return num_total_cuboid_corner_variables()
		+ num_total_cuboid_axis_variables()
		+ num_total_reflection_symmetry_group_variables()
		+ num_total_rotation_symmetry_group_variables();
}

NLPVectorExpression MeshCuboidNonLinearSolver::create_vector_variable(
	const std::pair<Index, Index>& _index_size_pair)
{
	assert(_index_size_pair.second > 0);
	NLPVectorExpression variable(_index_size_pair.second);
	for (unsigned int i = 0; i < _index_size_pair.second; ++i)
		variable[i] = NLPExpression(1.0, _index_size_pair.first + i);
	return variable;
}

std::pair<Index, Index> MeshCuboidNonLinearSolver::get_cuboid_corner_variable_index_size(
	unsigned int _cuboid_index, unsigned int _corner_index) const
{
	const unsigned int dimension = 3;
	assert(_corner_index < MeshCuboid::k_num_corners);
	
	unsigned int index = cuboid_corner_variable_start_index_
		+ num_cuboid_corner_variables_ * _cuboid_index
		+ dimension * _corner_index;
	return std::make_pair(index, dimension);
}

std::pair<Index, Index> MeshCuboidNonLinearSolver::get_cuboid_axis_variable_index_size(
	unsigned int _cuboid_index, unsigned int _axis_index) const
{
	const unsigned int dimension = 3;
	assert(_axis_index < dimension);
	
	unsigned int index = cuboid_axis_variable_start_index_
		+ num_cuboid_axis_variables_ * _cuboid_index
		+ dimension * _axis_index;
	return std::make_pair(index, dimension);
}

std::pair<Index, Index> MeshCuboidNonLinearSolver::get_reflection_symmetry_group_variable_n_index_size(
	unsigned int _symmetry_group_index) const
{
	const unsigned int dimension = 3;
	assert(_symmetry_group_index < num_reflection_symmetry_groups_);

	unsigned int index = reflection_symmetry_group_variable_start_index_
		+ num_reflection_symmetry_group_variables_ * _symmetry_group_index
		+ 0;
	return std::make_pair(index, dimension);
}

std::pair<Index, Index> MeshCuboidNonLinearSolver::get_reflection_symmetry_group_variable_t_index_size(
	unsigned int _symmetry_group_index) const
{
	const unsigned int dimension = 3;
	assert(_symmetry_group_index < num_reflection_symmetry_groups_);

	unsigned int index = reflection_symmetry_group_variable_start_index_
		+ num_reflection_symmetry_group_variables_ * _symmetry_group_index
		+ dimension;
	return std::make_pair(index, 1);
}

std::pair<Index, Index> MeshCuboidNonLinearSolver::get_rotation_symmetry_group_variable_n_index_size(
	unsigned int _symmetry_group_index) const
{
	const unsigned int dimension = 3;
	assert(_symmetry_group_index < num_rotation_symmetry_groups_);

	unsigned int index = rotation_symmetry_group_variable_start_index_
		+ num_rotation_symmetry_group_variables_ * _symmetry_group_index
		+ 0;
	return std::make_pair(index, dimension);
}

std::pair<Index, Index> MeshCuboidNonLinearSolver::get_rotation_symmetry_group_variable_t_index_size(
	unsigned int _symmetry_group_index) const
{
	const unsigned int dimension = 3;
	assert(_symmetry_group_index < num_rotation_symmetry_groups_);

	unsigned int index = rotation_symmetry_group_variable_start_index_
		+ num_rotation_symmetry_group_variables_ * _symmetry_group_index
		+ dimension;
	return std::make_pair(index, dimension);
}

NLPVectorExpression MeshCuboidNonLinearSolver::create_cuboid_corner_variable(
	unsigned int _cuboid_index, unsigned int _corner_index) const
{
	return create_vector_variable(get_cuboid_corner_variable_index_size(
		_cuboid_index, _corner_index));
}

NLPVectorExpression MeshCuboidNonLinearSolver::create_cuboid_axis_variable(
	unsigned int _cuboid_index, unsigned int _axis_index) const
{
	return create_vector_variable(get_cuboid_axis_variable_index_size(
		_cuboid_index, _axis_index));
}

NLPVectorExpression MeshCuboidNonLinearSolver::create_reflection_symmetry_group_variable_n(
	unsigned int _symmetry_group_index) const
{
	return create_vector_variable(get_reflection_symmetry_group_variable_n_index_size(
		_symmetry_group_index));
}

NLPVectorExpression MeshCuboidNonLinearSolver::create_reflection_symmetry_group_variable_t(
	unsigned int _symmetry_group_index) const
{
	return create_vector_variable(get_reflection_symmetry_group_variable_t_index_size(
		_symmetry_group_index));
}

NLPVectorExpression MeshCuboidNonLinearSolver::create_rotation_symmetry_group_variable_n(
	unsigned int _symmetry_group_index) const
{
	return create_vector_variable(get_rotation_symmetry_group_variable_n_index_size(
		_symmetry_group_index));
}

NLPVectorExpression MeshCuboidNonLinearSolver::create_rotation_symmetry_group_variable_t(
	unsigned int _symmetry_group_index) const
{
	return create_vector_variable(get_rotation_symmetry_group_variable_t_index_size(
		_symmetry_group_index));
}

NLPEigenQuadFunction* MeshCuboidNonLinearSolver::create_quadratic_energy_function(
	const Eigen::MatrixXd& _quadratic_term,
	const Eigen::VectorXd& _linear_term,
	const double _constant_term)
{
	if (_quadratic_term.rows() != num_total_cuboid_corner_variables()
		|| _quadratic_term.cols() != num_total_cuboid_corner_variables()
		|| _linear_term.rows() != num_total_cuboid_corner_variables())
	{
		std::cerr << "Error: " << std::endl;
		return NULL;
	}

	quadratic_term_ = Eigen::MatrixXd::Zero(num_total_variables(), num_total_variables());
	linear_term_ = Eigen::VectorXd::Zero(num_total_variables());
	constant_term_ = _constant_term;

	quadratic_term_.block(0, 0,
		num_total_cuboid_corner_variables(), num_total_cuboid_corner_variables()) = _quadratic_term;
	linear_term_.segment(0, num_total_cuboid_corner_variables()) = _linear_term;

	//
	add_reflection_symmetry_group_energy_functions(quadratic_term_, linear_term_, constant_term_);
	add_rotation_symmetry_group_energy_functions(quadratic_term_, linear_term_, constant_term_);
	//

	NLPEigenQuadFunction *function = new NLPEigenQuadFunction(
		quadratic_term_, linear_term_, constant_term_);

	return function;
}

void MeshCuboidNonLinearSolver::add_reflection_symmetry_group_energy_functions(
	Eigen::MatrixXd& _quadratic_term,
	Eigen::VectorXd& _linear_term,
	double &_constant_term)
{
	assert(_quadratic_term.rows() == num_total_variables());
	assert(_quadratic_term.cols() == num_total_variables());
	assert(_linear_term.rows() == num_total_variables());

	std::vector<ANNpointArray> cuboid_ann_points;
	std::vector<ANNkd_tree *> cuboid_ann_kd_tree;

	create_cuboid_sample_point_ann_trees(cuboid_ann_points, cuboid_ann_kd_tree);

	for (unsigned int symmetry_group_index = 0; symmetry_group_index < num_reflection_symmetry_groups_;
		++symmetry_group_index)
	{
		add_reflection_symmetry_group_energy_functions(symmetry_group_index,
			cuboid_ann_points, cuboid_ann_kd_tree,
			_quadratic_term, _linear_term, _constant_term);
	}

	delete_cuboid_sample_point_ann_trees(cuboid_ann_points, cuboid_ann_kd_tree);
}

void MeshCuboidNonLinearSolver::add_reflection_symmetry_group_energy_functions(
	const unsigned int _symmetry_group_index,
	const std::vector<ANNpointArray>& _cuboid_ann_points,
	const std::vector<ANNkd_tree *>& _cuboid_ann_kd_tree,
	Eigen::MatrixXd& _quadratic_term,
	Eigen::VectorXd& _linear_term,
	double &_constant_term)
{
	const MeshCuboidSymmetryGroup* symmetry_group = reflection_symmetry_groups_[_symmetry_group_index];
	assert(symmetry_group);

	std::list<MeshCuboidSymmetryGroup::WeightedPointPair> sample_point_pairs;
	symmetry_group->get_symmetric_sample_point_pairs(cuboids_,
		_cuboid_ann_points, _cuboid_ann_kd_tree, neighbor_distance_,
		sample_point_pairs);

	// Eq (1):
	// min {(I - nn^T)(x - y)}^2. Let d = (x - y).
	// Since (I - nn^T)^2 = (I - nn^T),
	// => min d^T(I - nn^T)d = d^Td - d^Tnn^Td = d^Td - n^Tdd^Tn.

	// Eq (2):
	// min {n^T(x + y) - 2t}^2.

	Eigen::Matrix3d A1 = Eigen::Matrix3d::Zero();
	Eigen::MatrixXd A2 = Eigen::MatrixXd::Zero(3 + 1, 3 + 1);

	Real sum_weight = 0.0;
	for (std::list<MeshCuboidSymmetryGroup::WeightedPointPair>::iterator it = sample_point_pairs.begin();
		it != sample_point_pairs.end(); ++it)
	{
		assert((*it).weight_ >= 0);
		sum_weight += (*it).weight_;
	}

	if (sum_weight == 0)
		return;

	for (std::list<MeshCuboidSymmetryGroup::WeightedPointPair>::iterator it = sample_point_pairs.begin();
		it != sample_point_pairs.end(); ++it)
	{
		Real weight = (*it).weight_ / sum_weight;

		Eigen::Vector3d sum_p, diff_p;
		for (int i = 0; i < 3; ++i)
		{
			sum_p[i] = (*it).p1_[i] + (*it).p2_[i];
			diff_p[i] = (*it).p1_[i] - (*it).p2_[i];
		}

		A1 += (weight * (-1) * diff_p * diff_p.transpose());

		_constant_term += (diff_p.transpose() * diff_p);

		Eigen::VectorXd b(3 + 1);
		b << sum_p, -2;
		A2 += (weight * b * b.transpose());
	}

	std::pair<Index, Index> index_size_pair;
	index_size_pair = get_reflection_symmetry_group_variable_n_index_size(_symmetry_group_index);

	_quadratic_term.block<3, 3>(index_size_pair.first, index_size_pair.first) +=
		symmetry_energy_term_weight_ * A1;
	
	// NOTE:
	// Variables 'n' and 't' are adjacent in the variable list.
	_quadratic_term.block<3 + 1, 3 + 1>(index_size_pair.first, index_size_pair.first) +=
		symmetry_energy_term_weight_ * A2;
}

void MeshCuboidNonLinearSolver::add_rotation_symmetry_group_energy_functions(
	Eigen::MatrixXd& _quadratic_term,
	Eigen::VectorXd& _linear_term,
	double &_constant_term)
{
	assert(_quadratic_term.rows() == num_total_variables());
	assert(_quadratic_term.cols() == num_total_variables());
	assert(_linear_term.rows() == num_total_variables());

	for (unsigned int symmetry_group_index = 0; symmetry_group_index < num_rotation_symmetry_groups_;
		++symmetry_group_index)
	{
		const MeshCuboidSymmetryGroup* symmetry_group = rotation_symmetry_groups_[symmetry_group_index];
		assert(symmetry_group);

		std::pair<Index, Index> index_size_pair;
		index_size_pair = get_rotation_symmetry_group_variable_t_index_size(symmetry_group_index);

		// Minimize \| t \|_2^2.
		_quadratic_term.block<3, 3>(index_size_pair.first, index_size_pair.first) +=
			Eigen::MatrixXd::Identity(index_size_pair.second, index_size_pair.second);
	}
}

void MeshCuboidNonLinearSolver::add_constraints(NLPFormulation &_formulation)
{
	add_cuboid_constraints(_formulation);
	add_reflection_symmetry_group_constraints(_formulation);
	add_rotation_symmetry_group_constraints(_formulation);
}

void MeshCuboidNonLinearSolver::add_cuboid_constraints(NLPFormulation &_formulation)
{
	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids_; ++cuboid_index)
	{
		add_cuboid_constraints(cuboid_index, _formulation);
	}
}

void MeshCuboidNonLinearSolver::add_reflection_symmetry_group_constraints(NLPFormulation &_formulation)
{
	for (unsigned int symmetry_group_index = 0; symmetry_group_index < num_reflection_symmetry_groups_;
		++symmetry_group_index)
	{
		add_reflection_symmetry_group_constraints(symmetry_group_index, _formulation);
	}
}

void MeshCuboidNonLinearSolver::add_rotation_symmetry_group_constraints(NLPFormulation &_formulation)
{
	for (unsigned int symmetry_group_index = 0; symmetry_group_index < num_rotation_symmetry_groups_;
		++symmetry_group_index)
	{
		add_rotation_symmetry_group_constraints(symmetry_group_index, _formulation);
	}
}

void MeshCuboidNonLinearSolver::add_cuboid_constraints(
	const unsigned int _cuboid_index,
	NLPFormulation &_formulation)
{
	const unsigned int dimension = 3;
	
	for (unsigned int axis_index = 0; axis_index < dimension; ++axis_index)
	{
		NLPVectorExpression axis_variable_1 = create_cuboid_axis_variable(_cuboid_index, axis_index);

		// Unit vector constraint.
		_formulation.add_constraint(NLPVectorExpression::dot_product(
			axis_variable_1, axis_variable_1), 1, 1);

		// Orthogonality constraint.
		if (dimension >= 2)
		{
			NLPVectorExpression axis_variable_2 = create_cuboid_axis_variable(_cuboid_index,
				(axis_index + 1) % dimension);
			assert(axis_variable_1.dimension() == axis_variable_2.dimension());

			// NOTE:
			// Relaxing equality constraint.
			// The equality constraints cause "too few degrees" errors.
			_formulation.add_constraint(NLPVectorExpression::dot_product(
				axis_variable_1, axis_variable_2), -1.0E-12, 1.0E-12);
		}
	}

	for (int axis_index = 0; axis_index < dimension; ++axis_index)
	{
		// NOTE:
		// Implemented only for 3 dimension.
		assert(dimension == 3);

		for (int axis_1 = 0; axis_1 < 2; ++axis_1)
		{
			for (int axis_2 = 0; axis_2 < 2; ++axis_2)
			{
				std::bitset<3> bits;
				bits[(axis_index + 1) % 3] = axis_1;
				bits[(axis_index + 2) % 3] = axis_2;

				bits[axis_index] = false;
				int corner_index_1 = bits.to_ulong();

				bits[axis_index] = true;
				int corner_index_2 = bits.to_ulong();

				NLPVectorExpression corner_variable_1 = create_cuboid_corner_variable(_cuboid_index, corner_index_1);
				NLPVectorExpression corner_variable_2 = create_cuboid_corner_variable(_cuboid_index, corner_index_2);

				// A cuboid edge is orthogonal with two axes.
				for (int other_axis_index = 0; other_axis_index < dimension; ++other_axis_index)
				{
					if (other_axis_index == axis_index) continue;
					NLPVectorExpression other_axis_variable = create_cuboid_axis_variable(_cuboid_index, other_axis_index);

					// n^T(x - y) = 0.
					NLPExpression expression = NLPVectorExpression::dot_product(other_axis_variable,
						corner_variable_1 - corner_variable_2);

					// NOTE:
					// Relaxing equality constraint.
					// The equality constraints cause "too few degrees" errors.
					_formulation.add_constraint(expression, -1.0E-12, 1.0E-12);
				}
			}
		}
	}
}

void MeshCuboidNonLinearSolver::add_reflection_symmetry_group_constraints(
	const unsigned int _symmetry_group_index,
	NLPFormulation &_formulation)
{
	const MeshCuboidReflectionSymmetryGroup *symmetry_group = reflection_symmetry_groups_[_symmetry_group_index];
	assert(symmetry_group);

	const unsigned int dimension = 3;
	NLPVectorExpression n_variable = create_reflection_symmetry_group_variable_n(_symmetry_group_index);


	std::vector< unsigned int > single_cuboid_indices;
	symmetry_group->get_single_cuboid_indices(cuboids_, single_cuboid_indices);

	for (std::vector< unsigned int >::const_iterator it = single_cuboid_indices.begin();
		it != single_cuboid_indices.end(); ++it)
	{
		unsigned int cuboid_index = (*it);
		
		NLPVectorExpression axis_variable = create_cuboid_axis_variable(cuboid_index,
			symmetry_group->get_aligned_global_axis_index());

		// The symmetry axis should be identical with the cuboid axis.
		NLPVectorExpression expression = n_variable - axis_variable;
		_formulation.add_constraint(expression, 0, 0);
	}


	std::vector< std::pair<unsigned int, unsigned int> > pair_cuboid_indices;
	symmetry_group->get_pair_cuboid_indices(cuboids_, pair_cuboid_indices);

	for (std::vector< std::pair<unsigned int, unsigned int> >::const_iterator it = pair_cuboid_indices.begin();
		it != pair_cuboid_indices.end(); ++it)
	{
		add_cuboid_reflection_constraints(
			(*it).first, (*it).second,
			_symmetry_group_index,
			symmetry_group->get_aligned_global_axis_index(),
			_formulation);
	}


	if (single_cuboid_indices.empty())
	{
		// Unit vector constraint.
		// NOTE:
		// If there is an identical cuboid axis,
		// the unit vector constraint is added from the cuboid constraint.
		NLPExpression expression = NLPVectorExpression::dot_product(
			n_variable, n_variable);
		_formulation.add_constraint(expression, 1, 1);
	}
}

void MeshCuboidNonLinearSolver::add_reflection_constraints(
	const NLPVectorExpression& _x_variable,
	const NLPVectorExpression& _y_variable,
	const NLPVectorExpression& _n_variable,
	const NLPVectorExpression& _t_variable,
	NLPFormulation &_formulation)
{
	const unsigned int dimension = _x_variable.dimension();
	assert(_y_variable.dimension() == dimension);
	assert(_n_variable.dimension() == dimension);
	assert(_t_variable.dimension() == 1);

	// n^T(x + y) - 2t = 0.
	NLPExpression expression_1 = NLPVectorExpression::dot_product(_n_variable,
		_x_variable + _y_variable)
		+ (_t_variable[0] * -2);

	// NOTE:
	// Relaxing equality constraint.
	// The equality constraints cause "too few degrees" errors.
	_formulation.add_constraint(expression_1, -1.0E-12, 1.0E-12);


	// n^T(x - y).
	NLPExpression temp = NLPVectorExpression::dot_product(_n_variable,
		_x_variable - _y_variable);

	// (I - nn^T)(x - y) = (x - y) - n(n^T(x - y)) = 0.
	NLPVectorExpression expression_2 = (_x_variable - _y_variable) - (_n_variable * temp);

	// NOTE:
	// Relaxing equality constraint.
	// The equality constraints cause "too few degrees" errors.
	_formulation.add_constraint(expression_2, -1.0E-12, 1.0E-12);
}

void MeshCuboidNonLinearSolver::add_cuboid_reflection_constraints(
	const unsigned int _cuboid_index_1,
	const unsigned int _cuboid_index_2,
	const unsigned int _symmetry_group_index,
	const unsigned int _reflection_axis_index,
	NLPFormulation &_formulation)
{
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

		add_reflection_constraints(
			create_cuboid_corner_variable(_cuboid_index_1, corner_index_1),
			create_cuboid_corner_variable(_cuboid_index_2, corner_index_2),
			create_reflection_symmetry_group_variable_n(_symmetry_group_index),
			create_reflection_symmetry_group_variable_t(_symmetry_group_index),
			_formulation);
	}
}

void MeshCuboidNonLinearSolver::add_rotation_symmetry_group_constraints(
	const unsigned int _symmetry_group_index,
	NLPFormulation &_formulation)
{
	const MeshCuboidRotationSymmetryGroup *symmetry_group = rotation_symmetry_groups_[_symmetry_group_index];
	assert(symmetry_group);

	const unsigned int dimension = 3;
	NLPVectorExpression n_variable = create_rotation_symmetry_group_variable_n(_symmetry_group_index);
	NLPVectorExpression t_variable = create_rotation_symmetry_group_variable_t(_symmetry_group_index);


	std::vector< unsigned int > single_cuboid_indices;
	symmetry_group->get_single_cuboid_indices(cuboids_, single_cuboid_indices);

	for (std::vector< unsigned int >::const_iterator it = single_cuboid_indices.begin();
		it != single_cuboid_indices.end(); ++it)
	{
		unsigned int cuboid_index = (*it);

		NLPVectorExpression axis_variable = create_cuboid_axis_variable(cuboid_index,
			symmetry_group->get_aligned_global_axis_index());

		// The symmetry axis should be identical with the cuboid axis.
		NLPVectorExpression expression_1 = n_variable - axis_variable;
		_formulation.add_constraint(expression_1, 0, 0);

		NLPVectorExpression cuboid_center_variable(dimension);
		for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
			cuboid_center_variable += create_cuboid_corner_variable(cuboid_index, corner_index);
		cuboid_center_variable *= (1.0 / MeshCuboid::k_num_corners);

		//  cross(n, (c - t)) = 0.
		NLPVectorExpression expression_2 = NLPVectorExpression::cross_product(
			n_variable, cuboid_center_variable - t_variable);
		_formulation.add_constraint(NLPVectorExpression::dot_product(expression_2, expression_2), -1.0E-12, 1.0E-12);
	}

	// Ignore pairs of cuboids.


	if (single_cuboid_indices.empty())
	{
		// Unit vector constraint.
		// NOTE:
		// If there is an identical cuboid axis,
		// the unit vector constraint is added from the cuboid constraint.
		NLPExpression expression = NLPVectorExpression::dot_product(
			n_variable, n_variable);
		_formulation.add_constraint(expression, 1, 1);
	}
}

bool MeshCuboidNonLinearSolver::compute_initial_values(const Eigen::VectorXd &_input,
	Eigen::VectorXd &_output)
{
	if (_input.rows() != num_total_cuboid_corner_variables())
	{
		std::cerr << "Error: " << std::endl;
		return false;
	}

	_output = Eigen::VectorXd::Zero(num_total_variables());
	_output.segment(0, num_total_cuboid_corner_variables()) = _input;

	compute_cuboid_axis_values(_output);
	compute_reflection_symmetry_group_values(_output);
	compute_rotation_symmetry_group_values(_output);

	return true;
}

void MeshCuboidNonLinearSolver::compute_cuboid_axis_values(Eigen::VectorXd &_values)
{
	assert(_values.rows() >= num_total_variables());

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids_; ++cuboid_index)
	{
		const MeshCuboid* cuboid = cuboids_[cuboid_index];
		assert(cuboid);
		for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
		{
			const MyMesh::Normal &axis = cuboid->get_bbox_axis(axis_index);

			std::pair<Index, Index> index_size_pair = get_cuboid_axis_variable_index_size(cuboid_index, axis_index);
			assert(index_size_pair.second == 3);
			for (unsigned int i = 0; i < index_size_pair.second; ++i)
				_values[index_size_pair.first + i] = axis[i];
		}
	}
}

void MeshCuboidNonLinearSolver::compute_reflection_symmetry_group_values(Eigen::VectorXd &_values)
{
	assert(_values.rows() >= num_total_variables());

	for (unsigned int symmetry_group_index = 0; symmetry_group_index < num_reflection_symmetry_groups_;
		++symmetry_group_index)
	{
		const MeshCuboidReflectionSymmetryGroup* symmetry_group = reflection_symmetry_groups_[symmetry_group_index];
		assert(symmetry_group);

		MyMesh::Normal n; double t;
		symmetry_group->get_reflection_plane(n, t);

		std::pair<Index, Index> index_size_pair;
		index_size_pair = get_reflection_symmetry_group_variable_n_index_size(symmetry_group_index);
		assert(index_size_pair.second == 3);
		for (unsigned int i = 0; i < index_size_pair.second; ++i)
			_values[index_size_pair.first + i] = n[i];

		index_size_pair = get_reflection_symmetry_group_variable_t_index_size(symmetry_group_index);
		assert(index_size_pair.second == 1);
		_values[index_size_pair.first] = t;
	}
}

void MeshCuboidNonLinearSolver::compute_rotation_symmetry_group_values(Eigen::VectorXd &_values)
{
	assert(_values.rows() >= num_total_variables());

	for (unsigned int symmetry_group_index = 0; symmetry_group_index < num_rotation_symmetry_groups_;
		++symmetry_group_index)
	{
		const MeshCuboidRotationSymmetryGroup* symmetry_group = rotation_symmetry_groups_[symmetry_group_index];
		assert(symmetry_group);

		MyMesh::Normal n; MyMesh::Point t;
		symmetry_group->get_rotation_axis(n, t);

		std::pair<Index, Index> index_size_pair;
		index_size_pair = get_rotation_symmetry_group_variable_n_index_size(symmetry_group_index);
		assert(index_size_pair.second == 3);
		for (unsigned int i = 0; i < index_size_pair.second; ++i)
			_values[index_size_pair.first + i] = n[i];

		index_size_pair = get_rotation_symmetry_group_variable_t_index_size(symmetry_group_index);
		assert(index_size_pair.second == 3);
		for (unsigned int i = 0; i < index_size_pair.second; ++i)
			_values[index_size_pair.first + i] = t[i];
	}
}

void MeshCuboidNonLinearSolver::optimize(
	const Eigen::MatrixXd& _quadratic_term,
	const Eigen::VectorXd& _linear_term,
	const double _constant_term,
	Eigen::VectorXd* _init_values_vec)
{
	NLPEigenQuadFunction* function = create_quadratic_energy_function(
		_quadratic_term, _linear_term, _constant_term);
	assert(function);
	NLPFormulation formulation(function);
	add_constraints(formulation);


	std::cout << "Energy: (";

	if (_init_values_vec)
	{
		Eigen::VectorXd all_init_values_vec;
		bool ret = compute_initial_values(*_init_values_vec, all_init_values_vec);
		std::cout << "initial = " << function->eval(all_init_values_vec) << ", ";
		formulation.set_values(all_init_values_vec.data());
	}


	// ---- //
	// Create a new instance of your nlp
	//  (use a SmartPtr, not raw)
	SmartPtr<TNLP> mynlp = new IPOPTSolver(&formulation);
	//SmartPtr<TNLP> mynlp = new HS071_NLP();

	// Create a new instance of IpoptApplication
	//  (use a SmartPtr, not raw)
	// We are using the factory, since this allows us to compile this
	// example with an Ipopt Windows DLL
	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	//app->RethrowNonIpoptException(true);

	// Change some options
	// Note: The following choices are only examples, they might not be
	//       suitable for your optimization problem.
	app->Options()->SetNumericValue("tol", 1e-8);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
	app->Options()->SetIntegerValue("print_level", 0);
	app->Options()->SetStringValue("fixed_variable_treatment", "relax_bounds");
	// app->Options()->SetStringValue("output_file", "ipopt.out");
	// The following overwrites the default name (ipopt.opt) of the
	// options file
	// app->Options()->SetStringValue("option_file_name", "ipopt.opt");

	// Initialize the IpoptApplication and process the options
	ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Solve_Succeeded) {
		std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
		assert(false);
	}

	// Ask Ipopt to solve the problem
	status = app->OptimizeTNLP(mynlp);

	if (status == Solve_Succeeded) {
		//std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
	}
	else {
		std::cout << std::endl << std::endl << "*** The problem FAILED!"
			<< " (" << static_cast<int>(status) << ")" << std::endl;
		//assert(false);
	}

	// As the SmartPtrs go out of scope, the reference count
	// will be decremented and the objects will automatically
	// be deleted.
	// ---- //


	// DEBUG.
	//formulation.print_constraint_evaluations();

	std::vector< Number > output;
	formulation.get_values(output);
	assert(output.size() == num_total_variables());
	std::cout << "final = " << function->eval(&(output[0])) << ")" << std::endl;


	// Update symmetry groups.
	update(output);
}

void MeshCuboidNonLinearSolver::update(const std::vector< Number >& _values)
{
	update_cuboids(_values);
	update_reflection_symmetry_groups(_values);
	update_rotation_symmetry_groups(_values);
}

void MeshCuboidNonLinearSolver::update_cuboids(const std::vector< Number >& _values)
{
	const unsigned int dimension = 3;

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids_; ++cuboid_index)
	{
		MeshCuboid *cuboid = cuboids_[cuboid_index];
		assert(cuboid);

		MyMesh::Point new_bbox_center(0.0);
		std::array<MyMesh::Point, MeshCuboid::k_num_corners> new_bbox_corners;
		std::array<MyMesh::Normal, dimension> new_bbox_axes;

		for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
		{
			std::pair<Index, Index> index_size_pair;
			index_size_pair = get_cuboid_corner_variable_index_size(cuboid_index, corner_index);
			assert(index_size_pair.second == dimension);
			for (unsigned int i = 0; i < index_size_pair.second; ++i)
				new_bbox_corners[corner_index][i] = _values[index_size_pair.first + i];

			new_bbox_center += new_bbox_corners[corner_index];
		}
		new_bbox_center = new_bbox_center / MeshCuboid::k_num_corners;

		for (unsigned int axis_index = 0; axis_index < dimension; ++axis_index)
		{
			std::pair<Index, Index> index_size_pair;
			index_size_pair = get_cuboid_axis_variable_index_size(cuboid_index, axis_index);
			assert(index_size_pair.second == dimension);
			for (unsigned int i = 0; i < index_size_pair.second; ++i)
				new_bbox_axes[axis_index][i] = _values[index_size_pair.first + i];

			new_bbox_axes[axis_index].normalize();


			// Flip axes.
			MyMesh::Point minus_axis_direction(0.0), plus_axis_direction(0.0);

			for (int axis_1 = 0; axis_1 < 2; ++axis_1)
			{
				for (int axis_2 = 0; axis_2 < 2; ++axis_2)
				{
					std::bitset<3> bits;
					bits[(axis_index + 1) % 3] = axis_1;
					bits[(axis_index + 2) % 3] = axis_2;

					bits[axis_index] = false;
					int minus_corner_index = bits.to_ulong();
					minus_axis_direction += new_bbox_corners[minus_corner_index];

					bits[axis_index] = true;
					int plus_corner_index = bits.to_ulong();
					plus_axis_direction += new_bbox_corners[plus_corner_index];
				}
			}

			// Compare the axis with the direction from the center to the average of corner points on +/- side.
			minus_axis_direction = (minus_axis_direction / 4.0 - new_bbox_center).normalized();
			plus_axis_direction = (plus_axis_direction / 4.0 - new_bbox_center).normalized();

			if (dot(new_bbox_axes[axis_index], minus_axis_direction)
				> dot(new_bbox_axes[axis_index], plus_axis_direction))
			{
				// Flip.
				new_bbox_axes[axis_index] = -new_bbox_axes[axis_index];
			}
		}

		cuboid->set_bbox_center(new_bbox_center);
		cuboid->set_bbox_corners(new_bbox_corners);
		cuboid->set_bbox_axes(new_bbox_axes, false);

		// Update cuboid surface points.
		for (unsigned int point_index = 0; point_index < cuboid->num_cuboid_surface_points();
			++point_index)
		{
			MeshCuboidSurfacePoint *cuboid_surface_point = cuboid->get_cuboid_surface_point(point_index);

			MyMesh::Point new_point(0.0);
			for (unsigned int corner_index = 0; corner_index < MeshCuboid::k_num_corners; ++corner_index)
				new_point += cuboid_surface_point->corner_weights_[corner_index] * new_bbox_corners[corner_index];

			cuboid_surface_point->point_ = new_point;
		}
	}
}

void MeshCuboidNonLinearSolver::update_reflection_symmetry_groups(const std::vector< Number >& _values)
{
	const unsigned int dimension = 3;

	for (unsigned int symmetry_group_index = 0; symmetry_group_index < num_reflection_symmetry_groups_;
		++symmetry_group_index)
	{
		MyMesh::Normal n; double t;

		std::pair<Index, Index> index_size_pair;
		index_size_pair = get_reflection_symmetry_group_variable_n_index_size(symmetry_group_index);
		assert(index_size_pair.second == dimension);
		for (unsigned int i = 0; i < index_size_pair.second; ++i)
			n[i] = _values[index_size_pair.first + i];

		index_size_pair = get_reflection_symmetry_group_variable_t_index_size(symmetry_group_index);
		assert(index_size_pair.second == 1);
		t = _values[index_size_pair.first];

		reflection_symmetry_groups_[symmetry_group_index]->set_reflection_plane(n, t);
	}
}

void MeshCuboidNonLinearSolver::update_rotation_symmetry_groups(const std::vector< Number >& _values)
{
	const unsigned int dimension = 3;

	for (unsigned int symmetry_group_index = 0; symmetry_group_index < num_rotation_symmetry_groups_;
		++symmetry_group_index)
	{
		MyMesh::Normal n; MyMesh::Point t;

		std::pair<Index, Index> index_size_pair;
		index_size_pair = get_rotation_symmetry_group_variable_n_index_size(symmetry_group_index);
		assert(index_size_pair.second == dimension);
		for (unsigned int i = 0; i < index_size_pair.second; ++i)
			n[i] = _values[index_size_pair.first + i];

		index_size_pair = get_rotation_symmetry_group_variable_t_index_size(symmetry_group_index);
		assert(index_size_pair.second == dimension);
		for (unsigned int i = 0; i < index_size_pair.second; ++i)
			t[i] = _values[index_size_pair.first + i];

		rotation_symmetry_groups_[symmetry_group_index]->set_rotation_axis(n, t);
	}
}

void MeshCuboidNonLinearSolver::create_cuboid_sample_point_ann_trees(
	std::vector<ANNpointArray>& _cuboid_ann_points,
	std::vector<ANNkd_tree *>& _cuboid_ann_kd_tree) const
{
	delete_cuboid_sample_point_ann_trees(_cuboid_ann_points, _cuboid_ann_kd_tree);
	_cuboid_ann_points.resize(num_cuboids_);
	_cuboid_ann_kd_tree.resize(num_cuboids_);

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids_; ++cuboid_index)
	{
		_cuboid_ann_kd_tree[cuboid_index] = NULL;
		_cuboid_ann_points[cuboid_index] = NULL;

		MeshCuboid *cuboid = cuboids_[cuboid_index];
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

		_cuboid_ann_kd_tree[cuboid_index] = ICP::create_kd_tree(cuboid_sample_points,
			_cuboid_ann_points[cuboid_index]);
		assert(_cuboid_ann_points[cuboid_index]);
		assert(_cuboid_ann_kd_tree[cuboid_index]);
	}
}

void MeshCuboidNonLinearSolver::delete_cuboid_sample_point_ann_trees(
	std::vector<ANNpointArray>& _cuboid_ann_points,
	std::vector<ANNkd_tree *>& _cuboid_ann_kd_tree) const
{
	if (_cuboid_ann_points.empty() && _cuboid_ann_kd_tree.empty())
		return;

	assert(_cuboid_ann_points.size() == num_cuboids_);
	assert(_cuboid_ann_kd_tree.size() == num_cuboids_);

	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids_; ++cuboid_index)
	{
		if (_cuboid_ann_points[cuboid_index]) annDeallocPts(_cuboid_ann_points[cuboid_index]);
		if (_cuboid_ann_kd_tree[cuboid_index]) delete _cuboid_ann_kd_tree[cuboid_index];
	}

	_cuboid_ann_points.clear();
	_cuboid_ann_kd_tree.clear();
}
