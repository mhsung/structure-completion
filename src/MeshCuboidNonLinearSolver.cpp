#include "MeshCuboidNonLinearSolver.h"

#include <bitset>

#include <Eigen/Eigenvalues> 


MeshCuboidNonLinearSolver::MeshCuboidNonLinearSolver(
	const std::vector<MeshCuboid *>& _cuboids,
	const std::vector<MeshCuboidSymmetryGroup *>& _symmetry_groups)
	: cuboids_(_cuboids)
	, symmetry_groups_(_symmetry_groups)
	, num_cuboid_corner_variables_(MeshCuboidAttributes::k_num_attributes)
	, num_cuboid_axis_variables_(3 * 3)
	, num_symmetry_group_variables_(4)
	, constant_term_(0)
{
	num_cuboids_ = _cuboids.size();
	//assert(num_cuboids_ > 0);

	num_symmetry_groups_ = _symmetry_groups.size();
	//assert(num_symmetry_groups_ > 0);

	cuboid_corner_variable_start_index_ = 0;
	cuboid_axis_variable_start_index_ = cuboid_corner_variable_start_index_
		+ num_total_cuboid_corner_variables();
	symmetry_group_variable_start_index_ = cuboid_axis_variable_start_index_
		+ num_total_cuboid_axis_variables();
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

unsigned int MeshCuboidNonLinearSolver::num_total_symmetry_gtoup_variables() const
{
	return num_symmetry_group_variables_ * num_symmetry_groups_;
}

unsigned int MeshCuboidNonLinearSolver::num_total_variables() const
{
	return num_total_cuboid_corner_variables()
		+ num_total_cuboid_axis_variables()
		+ num_total_symmetry_gtoup_variables();
}

unsigned int MeshCuboidNonLinearSolver::get_cuboid_corner_variable_start_index(
	unsigned int _cuboid_index) const
{
	unsigned int index = cuboid_corner_variable_start_index_
		+ num_cuboid_corner_variables_ * _cuboid_index;
	return index;
}

unsigned int MeshCuboidNonLinearSolver::get_cuboid_corner_variable_start_index(
	unsigned int _cuboid_index, unsigned int _corner_index) const
{
	unsigned int index = cuboid_corner_variable_start_index_
		+ num_cuboid_corner_variables_ * _cuboid_index
		+ 3 * _corner_index;
	return index;
}

unsigned int MeshCuboidNonLinearSolver::get_cuboid_axis_variable_start_index(
	unsigned int _cuboid_index) const
{
	unsigned int index = cuboid_axis_variable_start_index_
		+ num_cuboid_axis_variables_ * _cuboid_index;
	return index;
}

unsigned int MeshCuboidNonLinearSolver::get_cuboid_axis_variable_start_index(
	unsigned int _cuboid_index, unsigned int axis_index) const
{
	unsigned int index = cuboid_axis_variable_start_index_
		+ num_cuboid_axis_variables_ * _cuboid_index
		+ 3 * axis_index;
	return index;
}

unsigned int MeshCuboidNonLinearSolver::get_symmetry_group_variable_start_index(
	unsigned int _symmetry_group_index) const
{
	unsigned int index = symmetry_group_variable_start_index_
		+ num_symmetry_group_variables_ * _symmetry_group_index;
	return index;
}

NLPEigenQuadFunction* MeshCuboidNonLinearSolver::create_quadratic_function(
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

	NLPEigenQuadFunction *function = new NLPEigenQuadFunction(
		quadratic_term_, linear_term_, constant_term_);

	return function;
}

void MeshCuboidNonLinearSolver::add_cuboid_constraints(NLPFormulation &_formulation)
{
	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids_; ++cuboid_index)
	{
		add_cuboid_constraints(
			get_cuboid_corner_variable_start_index(cuboid_index),
			get_cuboid_axis_variable_start_index(cuboid_index),
			_formulation);
	}
}

void MeshCuboidNonLinearSolver::add_symmetry_group_constraints(NLPFormulation &_formulation)
{
	for (unsigned int symmetry_grou_index = 0; symmetry_grou_index < num_symmetry_groups_;
		++symmetry_grou_index)
	{
		add_symmetry_group_constraints(
			symmetry_grou_index, symmetry_groups_[symmetry_grou_index], _formulation);
	}
}

void MeshCuboidNonLinearSolver::add_symmetry_group_constraints(
	const unsigned int _symmetry_group_variable_index,
	const MeshCuboidSymmetryGroup *_symmetry_group,
	NLPFormulation &_formulation)
{
	assert(_symmetry_group);

	const unsigned int dimension = 3;
	const unsigned int n_variable_index = _symmetry_group_variable_index;


	std::vector< std::pair<unsigned int, unsigned int> > pair_cuboid_indices;
	_symmetry_group->get_pair_cuboid_indices(cuboids_, pair_cuboid_indices);

	for (std::vector< std::pair<unsigned int, unsigned int> >::const_iterator it = pair_cuboid_indices.begin();
		it != pair_cuboid_indices.end(); ++it)
	{
		add_cuboid_reflection_constraints((*it).first, (*it).second,
			_symmetry_group_variable_index,
			_symmetry_group->get_reflection_axis_index(),
			_formulation);
	}

	// Unit vector constraint.
	NLPExpression expression = NLPVectorExpression::dot_product(dimension,
		n_variable_index, n_variable_index);
	_formulation.add_constraint(expression, 1, 1);
}

void MeshCuboidNonLinearSolver::add_cuboid_constraints(
	const unsigned int _cuboid_corner_variable_index,
	const unsigned int _cuboid_axis_variable_index,
	NLPFormulation &_formulation)
{
	// Assume that corner point variables have successive indices
	// from '_corner_start_index' as follows:
	// [v1x v1y v1z v2x v2y v2z ... vnx vny xnz]

	// Also, assume that axis variables have successive indices
	// from '_axis_start_index' as follows:
	// [n1x n1y n1z n2x n2y n2z n3x n3y n3z]

	const unsigned int dimension = 3;

	for (unsigned int axis_index = 0; axis_index < dimension; ++axis_index)
	{
		int normal_index = _cuboid_axis_variable_index + dimension * axis_index;
		int next_normal_index = _cuboid_axis_variable_index + dimension * ((axis_index + 1) % dimension);

		// Unit vector constraint.
		_formulation.add_constraint(NLPVectorExpression::dot_product(dimension,
			normal_index, normal_index), 1, 1);

		// Orthogonality constraint.
		if (normal_index != next_normal_index)
		{
			_formulation.add_constraint(NLPVectorExpression::dot_product(dimension,
				normal_index, next_normal_index), 0, 0);
		}
	}


	for (int axis_index = 0; axis_index < dimension; ++axis_index)
	{
		// NOTICE:
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
				int corner_index_1 = _cuboid_corner_variable_index + dimension * bits.to_ulong();

				bits[axis_index] = true;
				int corner_index_2 = _cuboid_corner_variable_index + dimension * bits.to_ulong();

				//// A cuboid edge is parallel with a axis.
				//int normal_index = _normal_start_index + dimension * axis_index;
				//std::vector<NLPExpression> expr_1, expr_2;
				//NLPVectorExpression::cross_product(dimension, corner_index_1, normal_index, expr_1);
				//NLPVectorExpression::cross_product(dimension, corner_index_2, normal_index, expr_2);
				//assert(expr_1.size() == dimension);
				//assert(expr_2.size() == dimension);

				//for (int i = 0; i < dimension; ++i)
				//{
				//	NLPExpression expression = expr_1[i] - expr_2[i];
				//	// NOTICE:
				//	// A little bit relax the equality constraint.
				//	_formulation.add_constraint(expression, -1.0E-12, 1.0E-12);
				//}

				// A cuboid edge is orthogonal with two axes.
				for (int i = 0; i < 3; ++i)
				{
					if (i == axis_index) continue;
					int normal_index = _cuboid_axis_variable_index + dimension * i;

					NLPExpression expression
						= NLPVectorExpression::dot_product(dimension, corner_index_1, normal_index)
						- NLPVectorExpression::dot_product(dimension, corner_index_2, normal_index);
					_formulation.add_constraint(expression, 0, 0);
				}
			}
		}
	}
}

void MeshCuboidNonLinearSolver::add_reflection_constraints(
	const unsigned int _x_variable_index,
	const unsigned int _y_variable_index,
	const unsigned int _symmetry_group_variable_index,
	NLPFormulation &_formulation)
{
	// '(x, y)' is a symmetric point pair (3-dimensional),
	// 'n' is a unit normal vector of reflection plane (3-dimensional),
	// and 't' is a scalar value s.t. np = t for any 'p' on the reflection plane.

	const unsigned int dimension = 3;

	const unsigned int n_variable_index = _symmetry_group_variable_index;
	const unsigned int t_variable_index = _symmetry_group_variable_index + dimension;


	// n^T(x + y) - 2t = 0.
	NLPExpression expression_1
		= NLPVectorExpression::dot_product(dimension, n_variable_index, _x_variable_index)
		+ NLPVectorExpression::dot_product(dimension, n_variable_index, _y_variable_index)
		+ NLPExpression(-2, t_variable_index);

	_formulation.add_constraint(expression_1, 0, 0);

	// (I - nn^T)(x - y) = 0.
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		// n^T(x - y).
		NLPExpression expression_2
			= NLPVectorExpression::dot_product(dimension, n_variable_index, _x_variable_index)
			- NLPVectorExpression::dot_product(dimension, n_variable_index, _y_variable_index);

		// -nn^T(x - y) for each axis.
		expression_2 *= NLPExpression(-1, n_variable_index + axis_index);

		// (x - y) - nn^T(x - y) for each axis.
		expression_2 += NLPExpression(+1, _x_variable_index + axis_index);
		expression_2 += NLPExpression(-1, _y_variable_index + axis_index);

		_formulation.add_constraint(expression_2, 0, 0);
	}
}

void MeshCuboidNonLinearSolver::add_cuboid_reflection_constraints(
	const unsigned int _cuboid_index_1,
	const unsigned int _cuboid_index_2,
	const unsigned int _symmetry_group_variable_index,
	const unsigned int _reflection_axis_index,
	NLPFormulation &_formulation)
{
	// Assume that corner point variables have successive indices
	// from '_corner_start_index' as follows:
	// [v1x v1y v1z v2x v2y v2z ... vnx vny xnz]

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

		add_reflection_constraints(
			get_cuboid_corner_variable_start_index(_cuboid_index_1, corner_index_1),
			get_cuboid_corner_variable_start_index(_cuboid_index_2, corner_index_2),
			_symmetry_group_variable_index, _formulation);
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
			Eigen::Vector3d axis_vec;
			for (unsigned int i = 0; i < 3; ++i)
				axis_vec[i] = axis[i];

			_values.segment(get_cuboid_axis_variable_start_index(cuboid_index, axis_index), 3) = axis_vec;
		}
	}
}

void MeshCuboidNonLinearSolver::compute_symmetry_group_values(Eigen::VectorXd &_values)
{
	assert(_values.rows() >= num_total_variables());

	for (unsigned int symmetry_group_index = 0; symmetry_group_index < num_symmetry_groups_;
		++symmetry_group_index)
	{
		const MeshCuboidSymmetryGroup* symmetry_group = symmetry_groups_[symmetry_group_index];
		assert(symmetry_group);

		MyMesh::Normal n; double t;
		symmetry_group->get_reflection_plane(n, t);
		Eigen::VectorXd plane_vec = Eigen::VectorXd(num_symmetry_group_variables_);
		plane_vec << n[0], n[1], n[2], t;

		_values.segment(get_symmetry_group_variable_start_index(symmetry_group_index),
			num_symmetry_group_variables_) = plane_vec;
	}
}

Eigen::VectorXd MeshCuboidNonLinearSolver::solve_quadratic_programming_with_constraints(
	const Eigen::MatrixXd& _quadratic_term,
	const Eigen::VectorXd& _linear_term,
	const double _constant_term,
	Eigen::VectorXd* _init_values_vec)
{
	NLPEigenQuadFunction* function = create_quadratic_function(
		_quadratic_term, _linear_term, _constant_term);
	assert(function);
	NLPFormulation formulation(function);

	add_cuboid_constraints(formulation);
	add_symmetry_group_constraints(formulation);


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
	app->Options()->SetNumericValue("tol", 1e-7);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
	app->Options()->SetIntegerValue("print_level", 0);
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
		std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
		assert(false);
	}

	// As the SmartPtrs go out of scope, the reference count
	// will be decremented and the objects will automatically
	// be deleted.
	// ---- //


	std::vector< Number > output;
	formulation.get_values(output);
	assert(output.size() == num_total_variables());
	std::cout << "final = " << function->eval(&(output[0])) << ")" << std::endl;

	Eigen::VectorXd output_vec(num_total_cuboid_corner_variables());
	for (int i = 0; i < num_total_cuboid_corner_variables(); ++i)
		output_vec[i] = output[i];

	return output_vec;
}

