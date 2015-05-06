// Author: Minhyuk Sung
// Mar. 2015

#ifndef __NLP_EIGEN_QUAD_FUNCTION_H__
#define __NLP_EIGEN_QUAD_FUNCTION_H__

#include <list>
#include <string>
#include <vector>

#include "NLPFormulation.h"

#include <Eigen/Core>

class NLPSparseFunction;
class NLPSparseConstraint;
class NLPFormulation;


class NLPEigenQuadFunction : public NLPFunction
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	NLPEigenQuadFunction(const Eigen::MatrixXd &_quadratic_term,
		const Eigen::VectorXd &_linear_term, const double &_constant_term);
	~NLPEigenQuadFunction();

	virtual Number eval(const Number* _x) const;
	virtual void eval_gradient(const Number* _x, Number* _output) const;
	virtual void eval_hessian(const Number* _x, const Number _weight,
		Number* _output) const;

	virtual Number eval(const Eigen::VectorXd &_x) const;

private:
	const Eigen::MatrixXd quadratic_term_;
	const Eigen::VectorXd linear_term_;
	const double constant_term_;
};

#endif	// __NLP_EIGEN_QUAD_FUNCTION_H__
