// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: hs071_nlp.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "IPOPTSolver.h"

#include <cassert>
#include <iostream>

using namespace Ipopt;

// constructor
IPOPTSolver::IPOPTSolver(NLPFormulation *_formulation)
	: formulation_(_formulation)
{
	assert(formulation_);
}

//destructor
IPOPTSolver::~IPOPTSolver()
{}

// returns the size of the problem
bool IPOPTSolver::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
	Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	n = formulation_->num_variables();
	m = formulation_->num_contraints();

	nnz_jac_g = formulation_->nnz_constraint_gradients();
	nnz_h_lag = formulation_->nnz_hessian();

	// use the C style indexing (0-based)
	index_style = TNLP::C_STYLE;

	return true;
}

// returns the variable bounds
bool IPOPTSolver::get_bounds_info(Index n, Number* x_l, Number* x_u,
	Index m, Number* g_l, Number* g_u)
{
	assert(n == formulation_->num_variables());
	assert(m == formulation_->num_contraints());

	formulation_->get_variable_bounds(x_l, x_u);
	formulation_->get_constraint_bounds(g_l, g_u);

	return true;
}

// returns the initial point for the problem
bool IPOPTSolver::get_starting_point(Index n, bool init_x, Number* x,
	bool init_z, Number* z_L, Number* z_U,
	Index m, bool init_lambda,
	Number* lambda)
{
	// Here, we assume we only have starting values for x, if you code
	// your own NLP, you can provide starting values for the dual variables
	// if you wish
	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);

	formulation_->get_values(x);

	return true;
}

// returns the value of the objective function
bool IPOPTSolver::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
	assert(n == formulation_->num_variables());

	obj_value = formulation_->eval_function(x);

	return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IPOPTSolver::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	assert(n == formulation_->num_variables());

	formulation_->eval_function_gredient(x, grad_f);

	return true;
}

// return the value of the constraints: g(x)
bool IPOPTSolver::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
	assert(n == formulation_->num_variables());
	assert(m == formulation_->num_contraints());

	formulation_->eval_constraints(x, g);

	return true;
}

// return the structure or values of the jacobian
bool IPOPTSolver::eval_jac_g(Index n, const Number* x, bool new_x,
	Index m, Index nele_jac, Index* iRow, Index *jCol,
	Number* values)
{
	assert(n == formulation_->num_variables());
	assert(m == formulation_->num_contraints());

	formulation_->eval_constraint_gradients(x, iRow, jCol, values);

	return true;
}

//return the structure or values of the hessian
bool IPOPTSolver::eval_h(Index n, const Number* x, bool new_x,
	Number obj_factor, Index m, const Number* lambda,
	bool new_lambda, Index nele_hess, Index* iRow,
	Index* jCol, Number* values)
{
	assert(n == formulation_->num_variables());
	assert(m == formulation_->num_contraints());

	formulation_->eval_hessian(x, obj_factor, lambda, iRow, jCol, values);

	return true;
}

void IPOPTSolver::finalize_solution(SolverReturn status,
	Index n, const Number* x, const Number* z_L, const Number* z_U,
	Index m, const Number* g, const Number* lambda,
	Number obj_value,
	const IpoptData* ip_data,
	IpoptCalculatedQuantities* ip_cq)
{
	// here is where we would store the solution to variables, or write to a file, etc
	// so we could use the solution.

	formulation_->set_values(x);

	/*
	// For this example, we write the solution to the console
	std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
	for (Index i = 0; i < n; i++) {
		std::cout << "x[" << i << "] = " << x[i] << std::endl;
	}

	std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
	for (Index i = 0; i < n; i++) {
		std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
	}
	for (Index i = 0; i < n; i++) {
		std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
	}

	std::cout << std::endl << std::endl << "Objective value" << std::endl;
	std::cout << "f(x*) = " << obj_value << std::endl;

	std::cout << std::endl << "Final value of the constraints:" << std::endl;
	for (Index i = 0; i < m; i++) {
		std::cout << "g(" << i << ") = " << g[i] << std::endl;
	}
	*/
}
