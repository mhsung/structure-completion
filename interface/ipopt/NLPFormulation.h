// Author: Minhyuk Sung
// Mar. 2015

#ifndef __NLP_FORMULATION_H__
#define __NLP_FORMULATION_H__

#include <list>
#include <string>
#include <vector>

#include "NLPExpression.h"
#include "NLPVectorExpression.h"

#define NLP_BOUND_INFINITY	1e19

class NLPSparseFunction;
class NLPSparseConstraint;
class NLPFormulation;


class NLPFunction
{
public:
	NLPFunction(const Index _num_vars);
	~NLPFunction() {};

	virtual Number num_variables() const { return num_vars_; }
	virtual Number eval(const Number* _x) const = 0;
	virtual void eval_gradient(const Number* _x, Number* _output) const = 0;
	virtual void eval_hessian(const Number* _x, const Number _weight,
		Number* _output) const = 0;

	Index num_vars_;
};

class NLPConstraint
{
public:
	NLPConstraint(const Index _num_vars,
		const Number _lower_bound, const Number _upper_bound);
	~NLPConstraint() {};
	
	virtual Number num_variables() const { return num_vars_; }
	virtual Number nnz_gradients() const = 0;
	virtual Number eval(const Number* _x) const = 0;
	virtual void eval_gradient(const Number* _x,
		std::vector < std::pair < Index, Number > > &_output) const = 0;
	virtual void eval_hessian(const Number* _x, const Number _weight,
		Number* _output) const = 0;

	Index num_vars_;
	Number lower_bound_;
	Number upper_bound_;
};

class NLPSparseFunction : public NLPFunction
{
public:
	NLPSparseFunction(const Index _num_vars, const NLPExpression &_expression);
	~NLPSparseFunction();

	virtual Number eval(const Number* _x) const;
	virtual void eval_gradient(const Number* _x, Number* _output) const;
	virtual void eval_hessian(const Number* _x, const Number _weight,
		Number* _output) const;

private:
	NLPExpression expression_;
	NLPSparseGradientExpression gradient_;
	NLPSparseHessianExpression hessian_;
};

class NLPSparseConstraint : public NLPConstraint
{
public:
	NLPSparseConstraint(const Index _num_vars,
		const Number _lower_bound, const Number _upper_bound,
		const NLPExpression &_expression);
	~NLPSparseConstraint();
	
	virtual Number nnz_gradients() const;
	virtual Number eval(const Number* _x) const;
	virtual void eval_gradient(const Number* _x,
		std::vector < std::pair < Index, Number > > &_output) const;
	virtual void eval_hessian(const Number* _x, const Number _weight,
		Number* _output) const;

private:
	NLPExpression expression_;
	NLPSparseGradientExpression gradient_;
	NLPSparseHessianExpression hessian_;
};


class NLPFormulation
{
public:
	NLPFormulation(NLPFunction *_function);
	NLPFormulation(const std::vector<NLPFunction *> &_functions);
	~NLPFormulation();

	Index num_variables() const { return num_vars_; }
	Index num_contraints() const { return constraints_.size(); };

	Number nnz_constraint_gradients() const;
	Number nnz_hessian() const;

	void get_variable_bounds(Number* _lower_bound, Number* _upper_bound) const;
	void get_values(Number* _values) const;
	void get_values(std::vector< Number > &_values) const;
	void get_constraint_bounds(Number* _lower_bound, Number* _upper_bound) const;

	bool set_variable_bounds(const std::vector< Number > &_lower_bound,
		const std::vector< Number > &_upper_bound);
	void set_values(const Number* _x);
	bool set_values(const std::vector< Number > &_values);

	bool add_constraint(const NLPExpression &_expression,
		Number _lower_bound, Number _upper_bound);

	bool add_constraint(const NLPVectorExpression &_vector_expression,
		Number _lower_bound, Number _upper_bound);

	Number eval_function(const Number* _x) const;
	void eval_function_gredient(const Number* _x, Number* _output) const;
	void eval_function_hessian(const Number* _x, const Number _obj_factor,
		Number* _output) const;

	void eval_constraints(const Number* _x, Number* _output) const;
	void eval_constraint_gradients(const Number* _x,
		Index* _constraint_indices, Index *_variable_indices, Number* _output) const;	
	void eval_constraint_hessian(const Number* _x, const Number* _lambda,
		Number* _output) const;

	void eval_hessian(const Number* _x,
		const Number _obj_factor, const Number* _lambda,
		Index* _variable_indices_i, Index *_variable_indices_j, Number* _output) const;

	void print_constraint_evaluations() const;
	void print_constraint_evaluations(const NLPConstraint* _constraint) const;


private:
	Index num_vars_;
	std::vector< NLPFunction *> functions_;
	std::vector< Number > lower_bound_;
	std::vector< Number > upper_bound_;
	std::vector< Number > values_;
	std::vector< NLPConstraint* > constraints_;
};


#endif	// __NLP_FORMULATION_H__
