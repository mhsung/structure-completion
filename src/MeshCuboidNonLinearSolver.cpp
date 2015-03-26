#include "MeshCuboidNonLinearSolver.h"

#include <bitset>

#include <Eigen/Eigenvalues> 

// IPOPT.
#include "NLPFormulation.h"
#include "NLPEigenQuadFunction.h"
#include "NLPVectorExpression.h"
#include "IPOPTSolver.h"
#include "IpIpoptApplication.hpp"


void add_cuboid_constraints(
	const unsigned int _corner_variable_start_index,
	const unsigned int _axis_variable_start_index,
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
		int normal_index = _axis_variable_start_index + dimension * axis_index;
		int next_normal_index = _axis_variable_start_index + dimension * ((axis_index + 1) % dimension);

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
				int corner_index_1 = _corner_variable_start_index + dimension * bits.to_ulong();

				bits[axis_index] = true;
				int corner_index_2 = _corner_variable_start_index + dimension * bits.to_ulong();

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
					int normal_index = _axis_variable_start_index + dimension * i;

					NLPExpression expression
						= NLPVectorExpression::dot_product(dimension, corner_index_1, normal_index)
						- NLPVectorExpression::dot_product(dimension, corner_index_2, normal_index);
					_formulation.add_constraint(expression, 0, 0);
				}
			}
		}
	}
}

void add_reflection_constraints(
	const unsigned int _x_index, const unsigned int _y_index,
	const unsigned int _reflection_variable_start_index,
	NLPFormulation &_formulation)
{
	// '(x, y)' is a symmetric point pair (3-dimensional),
	// 'n' is a unit normal vector of reflection plane (3-dimensional),
	// and 't' is a scalar value s.t. np = t for any 'p' on the reflection plane.

	const unsigned int dimension = 3;

	const unsigned int n_index = _reflection_variable_start_index;
	const unsigned int t_index = _reflection_variable_start_index + dimension;

	// n^Tn = 1.
	//NLPExpression expression_1 = NLPVectorExpression::dot_product(dimension,
	//	_n_index, _n_index);

	//_formulation.add_constraint(expression_1, 1, 1);

	// n^T(x + y) - 2t = 0.
	NLPExpression expression_1
		= NLPVectorExpression::dot_product(dimension, n_index, _x_index)
		+ NLPVectorExpression::dot_product(dimension, n_index, _y_index)
		+ NLPExpression(-2, t_index);

	_formulation.add_constraint(expression_1, 0, 0);

	// (I - nn^T)(x - y) = 0.
	for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
	{
		// n^T(x - y).
		NLPExpression expression_2
			= NLPVectorExpression::dot_product(dimension, n_index, _x_index)
			- NLPVectorExpression::dot_product(dimension, n_index, _y_index);

		// -nn^T(x - y) for each axis.
		expression_2 *= NLPExpression(-1, n_index + axis_index);

		// (x - y) - nn^T(x - y) for each axis.
		expression_2 += NLPExpression(+1, _x_index + axis_index);
		expression_2 += NLPExpression(-1, _y_index + axis_index);

		_formulation.add_constraint(expression_2, 0, 0);
	}
}

void add_cuboid_reflection_constraints(
	const unsigned int _corner_variable_start_index_1,
	const unsigned int _corner_variable_start_index_2,
	const unsigned int _reflection_variable_start_index,
	const unsigned _reflection_axis_index,
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
			_corner_variable_start_index_1 + dimension * corner_index_1,
			_corner_variable_start_index_2 + dimension * corner_index_2,
			_reflection_variable_start_index, _formulation);
	}
}

void add_cuboid_reflection_corner_points(
	const MeshCuboid *cuboid_1, const MeshCuboid *cuboid_2,
	const std::vector<MeshCuboid *>& _cuboids,
	const unsigned int _cuboid_index_1,
	const unsigned int _cuboid_index_2,
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

void compute_reflection_plane(
	const std::list< std::pair < MyMesh::Point, MyMesh::Point > > &_point_pairs,
	Eigen::Vector3d &_n, double &_t)
{
	assert(!_point_pairs.empty());

	Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
	Eigen::Vector3d b = Eigen::Vector3d::Zero();

	for (std::list< std::pair < MyMesh::Point, MyMesh::Point > >::const_iterator it = _point_pairs.begin();
		it != _point_pairs.end(); ++it)
	{
		Eigen::Vector3d x, y;
		for (int i = 0; i < 3; ++i)
		{
			x[i] = (*it).first[i];
			y[2] = (*it).second[i];
		}

		Eigen::Vector3d d = (x - y);
		A += (d * d.transpose());
		b += (0.5 * (x + y));
	}

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(A);
	_n = es.eigenvectors().col(3 - 1);
	_t = (_n.transpose() * b);
	_t /= _point_pairs.size();
}

Eigen::VectorXd solve_quadratic_programming_with_constraints(
	const std::vector<MeshCuboid *>& _cuboids,
	const Eigen::MatrixXd& _quadratic_term,
	const Eigen::VectorXd& _linear_term,
	const double _constant_term,
	Eigen::VectorXd* _init_values_vec)
{
	unsigned int num_variables = 0;

	const unsigned int cuboid_corner_variable_start_index = 0;
	const unsigned int num_cuboid_corner_variables = _quadratic_term.cols();
	num_variables += num_cuboid_corner_variables;

	assert(_quadratic_term.rows() == num_cuboid_corner_variables);
	assert(_linear_term.rows() == num_cuboid_corner_variables);
	assert((num_cuboid_corner_variables % MeshCuboidAttributes::k_num_attributes) == 0);

	const int num_cuboids = num_cuboid_corner_variables / MeshCuboidAttributes::k_num_attributes;
	assert(num_cuboids == _cuboids.size());
	assert(num_cuboids > 0);


	// TEST.
	// Symmetry axes.
	const int num_symmetry_axes = 0;
	//


	// Additional variables.
	const unsigned int cuboid_axis_variable_start_index = num_variables;
	const unsigned int num_cuboid_axis_variables = (3 * 3 * num_cuboids);
	num_variables += num_cuboid_axis_variables;

	const unsigned int symmetry_axis_variables_start_index = num_variables;
	const unsigned int num_symmetry_axis_variables = (4 * num_symmetry_axes);
	num_variables += num_symmetry_axis_variables;

	Eigen::MatrixXd new_quadratic_term = Eigen::MatrixXd::Zero(num_variables, num_variables);
	Eigen::VectorXd new_linear_term = Eigen::VectorXd::Zero(num_variables);

	new_quadratic_term.block(0, 0, num_cuboid_corner_variables, num_cuboid_corner_variables) = _quadratic_term;
	new_linear_term.segment(0, num_cuboid_corner_variables) = _linear_term;


	// Create the energy function.
	NLPEigenQuadFunction *function = new NLPEigenQuadFunction(
		//_quadratic_term, _linear_term, _constant_term);
		new_quadratic_term, new_linear_term, _constant_term);
	NLPFormulation formulation(function);


	// Add cuboid constraints.
	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		unsigned int corner_variable_start_index = cuboid_corner_variable_start_index
			+ cuboid_index * MeshCuboidAttributes::k_num_attributes;
		unsigned int axis_variable_start_index = cuboid_axis_variable_start_index
			+ 3 * 3 * cuboid_index;
		add_cuboid_constraints(corner_variable_start_index, axis_variable_start_index, formulation);
	}
	

	if (_init_values_vec)
	{
		assert((*_init_values_vec).size() == num_cuboid_corner_variables);
		Eigen::VectorXd new_init_values_vec = Eigen::VectorXd::Zero(num_variables);
		new_init_values_vec.segment(0, num_cuboid_corner_variables) = (*_init_values_vec);
		
		// Initialize axis values.
		for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
		{
			unsigned int axis_start_index = num_cuboid_corner_variables + 3 * 3 * cuboid_index;
			for (unsigned int axis_index = 0; axis_index < 3; ++axis_index)
			{
				const MyMesh::Normal &axis = _cuboids[cuboid_index]->get_bbox_axis(axis_index);
				for (unsigned int i = 0; i < 3; ++i)
				{
					new_init_values_vec(axis_start_index + 3 * axis_index + i) = axis[i];
				}
				
			}
		}

		std::cout << "initial = " << function->eval(new_init_values_vec) << ", ";
		formulation.set_values(new_init_values_vec.data());
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
	assert(output.size() == num_variables);
	std::cout << "final = " << function->eval(&(output[0])) << ")" << std::endl;

	//Eigen::Map<Eigen::VectorXd> output_vec(&(output[0]), num_corner_vars);
	Eigen::VectorXd output_vec(num_cuboid_corner_variables);

	for (int i = 0; i < num_cuboid_corner_variables; ++i)
		output_vec[i] = output[i];

	return output_vec;
}

