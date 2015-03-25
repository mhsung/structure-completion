#include "MeshCuboidNonLinearSolver.h"

#include <bitset>

// IPOPT.
#include "NLPFormulation.h"
#include "NLPEigenQuadFunction.h"
#include "NLPVectorExpression.h"
#include "IPOPTSolver.h"
#include "IpIpoptApplication.hpp"


void add_cuboid_constraints(const unsigned int _corner_start_index,
	const unsigned int _axis_start_index,
	NLPFormulation &_formulation)
{
	const unsigned int dimension = 3;

	for (unsigned int axis_index = 0; axis_index < dimension; ++axis_index)
	{
		int normal_index = _axis_start_index + dimension * axis_index;
		int next_normal_index = _axis_start_index + dimension * ((axis_index + 1) % dimension);

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
				int corner_index_1 = _corner_start_index + dimension * bits.to_ulong();

				bits[axis_index] = true;
				int corner_index_2 = _corner_start_index + dimension * bits.to_ulong();

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
					int normal_index = _axis_start_index + dimension * i;

					NLPExpression expression
						= NLPVectorExpression::dot_product(dimension, corner_index_1, normal_index)
						- NLPVectorExpression::dot_product(dimension, corner_index_2, normal_index);
					_formulation.add_constraint(expression, 0, 0);
				}
			}
		}
	}
}

Eigen::VectorXd solve_quadratic_programming_with_constraints(
	const std::vector<MeshCuboid *>& _cuboids,
	const Eigen::MatrixXd& _quadratic_term,
	const Eigen::VectorXd& _linear_term,
	const double _constant_term,
	Eigen::VectorXd* _init_values_vec)
{
	const int num_corner_vars = _quadratic_term.cols();
	assert(_quadratic_term.rows() == num_corner_vars);
	assert(_linear_term.rows() == num_corner_vars);
	assert((num_corner_vars % MeshCuboidAttributes::k_num_attributes) == 0);

	const int num_cuboids = num_corner_vars / MeshCuboidAttributes::k_num_attributes;
	assert(num_cuboids == _cuboids.size());
	assert(num_cuboids > 0);

	// Add three 3-dimensional axes values for each cuboid.
	const int num_vars = num_corner_vars + (3 * 3 * num_cuboids);

	Eigen::MatrixXd new_quadratic_term = Eigen::MatrixXd::Zero(num_vars, num_vars);
	Eigen::VectorXd new_linear_term = Eigen::VectorXd::Zero(num_vars);

	new_quadratic_term.block(0, 0, num_corner_vars, num_corner_vars) = _quadratic_term;
	new_linear_term.segment(0, num_corner_vars) = _linear_term;


	NLPEigenQuadFunction *function = new NLPEigenQuadFunction(
		//_quadratic_term, _linear_term, _constant_term);
		new_quadratic_term, new_linear_term, _constant_term);
	NLPFormulation formulation(function);


	// Add cuboid constraints.
	for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
	{
		add_cuboid_constraints(0 + cuboid_index * MeshCuboidAttributes::k_num_attributes,
			num_corner_vars + 3 * 3 * cuboid_index, formulation);
	}
	

	if (_init_values_vec)
	{
		assert((*_init_values_vec).size() == num_corner_vars);
		Eigen::VectorXd new_init_values_vec = Eigen::VectorXd::Zero(num_vars);
		new_init_values_vec.segment(0, num_corner_vars) = (*_init_values_vec);
		
		// Initialize axis values.
		for (unsigned int cuboid_index = 0; cuboid_index < num_cuboids; ++cuboid_index)
		{
			unsigned int axis_start_index = num_corner_vars + 3 * 3 * cuboid_index;
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
	assert(output.size() == num_vars);
	std::cout << "final = " << function->eval(&(output[0])) << ")" << std::endl;

	//Eigen::Map<Eigen::VectorXd> output_vec(&(output[0]), num_corner_vars);
	Eigen::VectorXd output_vec(num_corner_vars);

	for (int i = 0; i < num_corner_vars; ++i)
		output_vec[i] = output[i];

	return output_vec;
}

