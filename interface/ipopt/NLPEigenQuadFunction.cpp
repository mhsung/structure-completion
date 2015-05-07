#include "NLPEigenQuadFunction.h"

#include <stdio.h>
#include <string.h>
#include <cassert>


NLPEigenQuadFunction::NLPEigenQuadFunction(const Eigen::MatrixXd &_quadratic_term_,
	const Eigen::VectorXd &_linear_term, const double &_constant_term)
	: NLPFunction(0)
	, quadratic_term_(_quadratic_term_)
	, linear_term_(_linear_term)
	, constant_term_(_constant_term)
{
	num_vars_ = quadratic_term_.cols();
	assert(quadratic_term_.rows() == num_vars_);
	assert(linear_term_.rows() == num_vars_);
}

NLPEigenQuadFunction::~NLPEigenQuadFunction()
{

}

Ipopt::Number NLPEigenQuadFunction::eval(const Number* _x) const
{
	// FIXME:
	// Change the following code in more efficient way.
	Eigen::VectorXd x(num_vars_);
	for (int i = 0; i < num_vars_; ++i)
		x[i] = static_cast<double>(_x[i]);

	double output = 0;
	output += x.transpose() * quadratic_term_ * x;
	output += 2 * linear_term_.transpose() * x;
	output += constant_term_;

	return static_cast<Number>(output);
}

Ipopt::Number NLPEigenQuadFunction::eval(const Eigen::VectorXd &_x) const
{
	assert(_x.rows() == num_vars_);

	double output = 0;
	output += _x.transpose() * quadratic_term_ * _x;
	output += 2 * linear_term_.transpose() * _x;
	output += constant_term_;

	return static_cast<Number>(output);
}

void NLPEigenQuadFunction::eval_gradient(const Number* _x, Number* _output) const
{
	// FIXME:
	// Change the following code in more efficient way.
	Eigen::VectorXd x(num_vars_);
	for (int i = 0; i < num_vars_; ++i)
		x[i] = static_cast<double>(_x[i]);

	Eigen::VectorXd output(num_vars_);
	output = 2 * (quadratic_term_ * x + linear_term_);

	for (int i = 0; i < num_vars_; ++i)
		_output[i] = static_cast<Number>(output[i]);
}

void NLPEigenQuadFunction::eval_hessian(const Number* _x,
	const Number _weight, Number* _output) const
{
	for (Index index_i = 0; index_i < num_vars_; ++index_i)
	{
		for (Index index_j = 0; index_j <= index_i; ++index_j)
		{
			assert(index_i <= quadratic_term_.rows());
			assert(index_j <= quadratic_term_.cols());
			assert(index_i >= index_j);

			Index index = index_i * (index_i + 1) / 2 + index_j;
			Number value = 0.5 * (quadratic_term_(index_i, index_j) +
				quadratic_term_(index_j, index_i));
			_output[index] += _weight * (2 * value);
		}
	}
}
