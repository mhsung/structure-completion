#include "NLPFormulation.h"

#include <stdio.h>
#include <string.h>
#include <cassert>
#include <iostream>
#include <sstream>

NLPFunction::NLPFunction(const Index _num_vars)
	: num_vars_(_num_vars)
{

}

NLPConstraint::NLPConstraint(const Index _num_vars,
	const Number _lower_bound, const Number _upper_bound)
	: num_vars_(_num_vars)
	, lower_bound_(_lower_bound)
	, upper_bound_(_upper_bound)
{

}

NLPSparseFunction::NLPSparseFunction(const Index _num_vars,
	const NLPExpression &_expression)
	: NLPFunction(_num_vars)
	, expression_(_expression)
{
	bool ret;
	ret = _expression.get_sparse_gradient(_num_vars, gradient_.indices_, gradient_.expressions_);
	assert(ret);
	ret = _expression.get_sparse_hessian(_num_vars, hessian_.indices_, hessian_.expressions_);
	assert(ret);
}

NLPSparseFunction::~NLPSparseFunction()
{

}

Ipopt::Number NLPSparseFunction::eval(const Number* _x) const
{
	return expression_.eval(_x);
}

void NLPSparseFunction::eval_gradient(const Number* _x, Number* _output) const
{
	// Initialize all output values to zero.
	memset(_output, 0, num_vars_ * sizeof(Number));

	unsigned int nnz_gradient = gradient_.expressions_.size();
	assert(gradient_.indices_.size() == nnz_gradient);

	for (unsigned int i = 0; i < nnz_gradient; ++i)
	{
		Index index = gradient_.indices_[i];
		const NLPExpression &expression = gradient_.expressions_[i];

		assert(index <= num_variables());
		_output[index] = expression.eval(_x);
	}
}

void NLPSparseFunction::eval_hessian(const Number* _x,
	const Number _weight, Number* _output) const
{
	// NOTE:
	// Make dense Hessian matrix.

	assert(_output != NULL);
	// Assume that output values are initialized.

	unsigned int nnz_funcion_hessian = hessian_.expressions_.size();
	assert(hessian_.indices_.size() == nnz_funcion_hessian);

	for (unsigned int i = 0; i < nnz_funcion_hessian; ++i)
	{
		Index index_i = hessian_.indices_[i].first;
		Index index_j = hessian_.indices_[i].second;
		const NLPExpression &expression = hessian_.expressions_[i];

		assert(index_i <= num_variables());
		assert(index_j <= num_variables());
		assert(index_i >= index_j);

		Index index = index_i * (index_i + 1) / 2 + index_j;
		_output[index] += _weight * expression.eval(_x);
	}
}

NLPSparseConstraint::NLPSparseConstraint(const Index _num_vars,
	const Number _lower_bound, const Number _upper_bound,
	const NLPExpression &_expression)
	: NLPConstraint(_num_vars, _lower_bound, _upper_bound)
	, expression_(_expression)
{
	bool ret;
	ret = _expression.get_sparse_gradient(_num_vars, gradient_.indices_, gradient_.expressions_);
	assert(ret);
	ret = _expression.get_sparse_hessian(_num_vars, hessian_.indices_, hessian_.expressions_);
	assert(ret);
}

NLPSparseConstraint::~NLPSparseConstraint()
{

}

Ipopt::Number NLPSparseConstraint::nnz_gradients() const
{
	unsigned int nnz_gradient = gradient_.expressions_.size();
	assert(gradient_.indices_.size() == nnz_gradient);
	return nnz_gradient;
}

Ipopt::Number NLPSparseConstraint::eval(const Number* _x) const
{
	return expression_.eval(_x);
}

void NLPSparseConstraint::eval_gradient(const Number* _x,
	std::vector < std::pair < Index, Number > > &_output) const
{
	unsigned int nnz_gradient = gradient_.expressions_.size();
	assert(gradient_.indices_.size() == nnz_gradient);

	_output.clear();
	_output.reserve(nnz_gradient);

	for (unsigned int i = 0; i < nnz_gradient; ++i)
	{
		Index index = gradient_.indices_[i];
		const NLPExpression &expression = gradient_.expressions_[i];
		assert(index <= num_vars_);

		// NOTE:
		// If '_x' is NULL, ignore the output value.
		Number value = 0;
		if (_x) value = expression.eval(_x);

		_output.push_back(std::make_pair(index, value));
	}
}

void NLPSparseConstraint::eval_hessian(const Number* _x,
	const Number _weight, Number* _output) const
{
	// NOTE:
	// Make dense Hessian matrix.

	assert(_output != NULL);
	// Assume that output values are initialized.

	unsigned int nnz_funcion_hessian = hessian_.expressions_.size();
	assert(hessian_.indices_.size() == nnz_funcion_hessian);

	for (unsigned int i = 0; i < nnz_funcion_hessian; ++i)
	{
		Index index_i = hessian_.indices_[i].first;
		Index index_j = hessian_.indices_[i].second;
		const NLPExpression &expression = hessian_.expressions_[i];

		assert(index_i <= num_variables());
		assert(index_j <= num_variables());
		assert(index_i >= index_j);

		Index index = index_i * (index_i + 1) / 2 + index_j;
		_output[index] += _weight * expression.eval(_x);
	}
}

NLPFormulation::NLPFormulation(NLPFunction *_function)
{
	assert(_function);
	num_vars_ = _function->num_variables();

	functions_.push_back(_function);

	lower_bound_.resize(num_vars_, -NLP_BOUND_INFINITY);
	upper_bound_.resize(num_vars_, NLP_BOUND_INFINITY);
	values_.resize(num_vars_, 0);
}

NLPFormulation::NLPFormulation(const std::vector<NLPFunction *> &_functions)
{
	assert(!_functions.empty());
	functions_ = _functions;
	num_vars_ = functions_.front()->num_variables();

	for (std::vector<NLPFunction*>::iterator it = functions_.begin();
		it != functions_.end(); ++it)
	{
		assert((*it)->num_variables() == num_vars_);
	}

	lower_bound_.resize(num_vars_, -NLP_BOUND_INFINITY);
	upper_bound_.resize(num_vars_, NLP_BOUND_INFINITY);
	values_.resize(num_vars_, 0);
}

NLPFormulation::~NLPFormulation()
{
	for (std::vector<NLPFunction*>::iterator it = functions_.begin();
		it != functions_.end(); ++it)
	{
		assert(*it);
		delete (*it);
	}

	for (std::vector< NLPConstraint* >::iterator it = constraints_.begin();
		it != constraints_.end(); ++it)
	{
		assert(*it);
		delete (*it);
	}
}

Number NLPFormulation::nnz_constraint_gradients() const
{
	unsigned int count = 0;

	const unsigned int num_constraints = constraints_.size();
	for (unsigned int i = 0; i < num_constraints; i++)
	{
		assert(constraints_[i]);
		count += constraints_[i]->nnz_gradients();
	}

	return count;
}

Number NLPFormulation::nnz_hessian() const
{
	// NOTE:
	// Make dense Hessian matrix.
	return num_vars_ * (num_vars_ + 1) / 2;
}

bool NLPFormulation::set_variable_bounds(const std::vector< Number > &_lower_bound,
	const std::vector< Number > &_upper_bound)
{
	if (_lower_bound.size() != num_vars_)
	{
		std::cerr << "Error: The number of variables in the lower bound is not the same: ("
			<< _lower_bound .size() << " != " << num_vars_ << ")" << std::endl;
		return false;
	}
	else if (_upper_bound.size() != num_vars_)
	{
		std::cerr << "Error: The number of variables in the upper bound is not the same: ("
			<< _upper_bound.size() << " != " << num_vars_ << ")" << std::endl;
		return false;
	}

	lower_bound_.clear();
	upper_bound_.clear();
	lower_bound_ = _lower_bound;
	upper_bound_ = _upper_bound;

	return true;
}

void NLPFormulation::set_values(const Number* _x)
{
	for (Index i = 0; i < num_vars_; i++)
	{
		values_[i] = _x[i];
	}
}

bool NLPFormulation::set_values(const std::vector< Number > &_values)
{
	if (_values.size() != num_vars_)
	{
		std::cerr << "Error: The number of variables in the starting point is not the same: ("
			<< _values.size() << " != " << num_vars_ << ")" << std::endl;
		return false;
	}

	values_.clear();
	values_ = _values;

	return true;
}

void NLPFormulation::get_variable_bounds(Number* _lower_bound, Number* _upper_bound) const
{
	for (Index i = 0; i < num_vars_; i++)
	{
		_lower_bound[i] = lower_bound_[i];
		_upper_bound[i] = upper_bound_[i];
	}
}

void NLPFormulation::get_values(Number* _values) const
{
	for (Index i = 0; i < num_vars_; i++)
	{
		_values[i] = values_[i];
	}
}

void NLPFormulation::get_values(std::vector< Number > &_values) const
{
	_values.clear();
	_values = values_;
}

bool NLPFormulation::add_constraint(const NLPExpression &_expression,
	Number _lower_bound, Number _upper_bound)
{
	NLPSparseConstraint *new_constraint = new NLPSparseConstraint(
		num_vars_, _lower_bound, _upper_bound, _expression);
	constraints_.push_back(new_constraint);
	return true;
}

bool NLPFormulation::add_constraint(const NLPVectorExpression &_vector_expression,
	Number _lower_bound, Number _upper_bound)
{
	for (unsigned int i = 0; i < _vector_expression.dimension(); ++i)
	{
		NLPSparseConstraint *new_constraint = new NLPSparseConstraint(
			num_vars_, _lower_bound, _upper_bound, _vector_expression[i]);
		constraints_.push_back(new_constraint);
	}
	return true;
}

void NLPFormulation::get_constraint_bounds(Number* _lower_bound, Number* _upper_bound) const
{
	const unsigned int num_constraints = constraints_.size();
	for (unsigned int i = 0; i < num_constraints; i++)
	{
		assert(constraints_[i]);
		_lower_bound[i] = constraints_[i]->lower_bound_;
		_upper_bound[i] = constraints_[i]->upper_bound_;
	}
}

Number NLPFormulation::eval_function(const Number* _x) const
{
	Number output = 0;
	for (std::vector<NLPFunction*>::const_iterator it = functions_.begin();
		it != functions_.end(); ++it)
	{
		assert(*it);
		output += (*it)->eval(_x);
	}
	return output;
}

void NLPFormulation::eval_function_gredient(const Number* _x, Number* _output) const
{
	memset(_output, 0, num_vars_ * sizeof(Number));

	for (std::vector<NLPFunction*>::const_iterator it = functions_.begin();
		it != functions_.end(); ++it)
	{
		assert(*it);
		Number* temp = new Number[num_vars_];
		memset(temp, 0, num_vars_ * sizeof(Number));

		(*it)->eval_gradient(_x, temp);

		for (int i = 0; i < num_vars_; ++i)
			_output[i] += temp[i];
	}
}

void NLPFormulation::eval_function_hessian(const Number* _x,
	const Number _obj_factor, Number* _output) const
{
	// NOTE:
	// When evaluating Hessian, assume that '_output' is initialized.

	for (std::vector<NLPFunction*>::const_iterator it = functions_.begin();
		it != functions_.end(); ++it)
	{
		assert(*it);
		(*it)->eval_hessian(_x, _obj_factor, _output);
	}
}

void NLPFormulation::eval_constraints(const Number* _x, Number* _output) const
{
	const unsigned int num_constraints = constraints_.size();
	for (unsigned int i = 0; i < num_constraints; i++)
	{
		assert(constraints_[i]);
		_output[i] = constraints_[i]->eval(_x);
	}
}

void NLPFormulation::eval_constraint_gradients(const Number* _x,
	Index* _constraint_indices, Index *_variable_indices, Number* _output) const
{
	unsigned int count = 0;

	const unsigned int num_constraints = constraints_.size();
	for (unsigned int i = 0; i < num_constraints; i++)
	{
		assert(constraints_[i]);
		std::vector < std::pair < Index, Number > > constraint_output;
		constraints_[i]->eval_gradient(_x, constraint_output);

		for (unsigned int j = 0; j < constraint_output.size(); j++)
		{
			if (_output == NULL)
			{
				// Return the structure of the constraint gradients.
				_constraint_indices[count] = i;
				_variable_indices[count] = constraint_output[j].first;
			}
			else
			{
				// Return the values of the constraint gradients.
				_output[count] = constraint_output[j].second;
			}
			++count;
		}
	}
}

void NLPFormulation::eval_constraint_hessian(const Number* _x,
	const Number* _lambda, Number* _output) const
{
	const unsigned int num_constraints = constraints_.size();
	for (unsigned int i = 0; i < num_constraints; i++)
	{
		assert(constraints_[i]);
		constraints_[i]->eval_hessian(_x, _lambda[i], _output);
	}
}

void NLPFormulation::eval_hessian(const Number* _x,
	const Number _obj_factor, const Number* _lambda,
	Index* _variable_indices_i, Index *_variable_indices_j, Number* _output) const
{
	// NOTE:
	// Make dense Hessian matrix.

	if (_output == NULL)
	{
		// Return the structure of the constraint gradients.

		Index count = 0;
		for (Index index_i = 0; index_i < num_vars_; index_i++) {
			for (Index index_j = 0; index_j <= index_i; index_j++) {
				_variable_indices_i[count] = index_i;
				_variable_indices_j[count] = index_j;
				count++;
			}
		}
	}
	else
	{
		// Return the values of the constraint gradients.
                
        // Initialize all output values to zero.
        memset(_output, 0, nnz_hessian() * sizeof(Number));

		eval_function_hessian(_x, _obj_factor, _output);
		eval_constraint_hessian(_x, _lambda, _output);
	}
}

void NLPFormulation::print_constraint_evaluations() const
{
	const unsigned int num_constraints = constraints_.size();
	for (unsigned int constraint_index = 0; constraint_index < num_constraints; ++constraint_index)
	{
		std::cout << "[" << constraint_index << "]; ";
		print_constraint_evaluations(constraints_[constraint_index]);
	}
}

void NLPFormulation::print_constraint_evaluations(const NLPConstraint* _constraint) const
{
	assert(_constraint);
	Number output = _constraint->eval(&(values_[0]));

	if (output < _constraint->lower_bound_)
	{
		std::cout << "(Out of Range) " << output << " < "
			<< _constraint->lower_bound_ << "(L) <= "
			<< _constraint->upper_bound_ << "(H)" << std::endl;
	}
	else if (output > _constraint->upper_bound_)
	{
		std::cout << "(Out of Range) " << _constraint->lower_bound_ << "(L) <= "
			<< _constraint->upper_bound_ << "(H) < "
			<< output << std::endl;
	}
	else
	{
		std::cout << _constraint->lower_bound_ << "(L) <= "
			<< output << " <= "
			<< _constraint->upper_bound_ << "(H)" << std::endl;
	}
}
