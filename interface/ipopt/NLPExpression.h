// Author: Minhyuk Sung
// Mar. 2015

#ifndef __NLP_EXPRESSION_H__
#define __NLP_EXPRESSION_H__

#include <list>
#include <string>
#include <vector>

#include "IpTNLP.hpp"

using namespace Ipopt;
//typedef double Number;
//typedef int Index;

class NLPTerm;
class NLPExpression;


class NLPTerm
{
public:
	NLPTerm();
	NLPTerm(const Number _coeff);
	NLPTerm(const Number _coeff, const Index _index);
	NLPTerm(const Number _coeff, const std::vector<Index> &_indices);
	NLPTerm(const NLPTerm &_term);
	~NLPTerm();


	void set_coeff(Number _coeff);

	// FIXME:
	// Should be replaced with *= operator.
	void multiply_coeff(Number _coeff);

	void multiply_var(Index _index, int _power);

	bool empty() const;

	Number eval(const Number* _x) const;

	bool get_gradient(const Index _num_vars, std::vector< NLPTerm > &_gradient_terms) const;
	bool get_hessian(const Index _num_vars,
		std::vector< std::vector< NLPTerm > > &_hessian_terms) const;

	bool get_sparse_gradient(const Index _num_vars, std::vector< Index > &_gradient_indices,
		std::vector< NLPTerm > &_gradient_terms) const;
	bool get_sparse_hessian(const Index _num_vars,
		std::vector< std::pair<Index, Index> > &_hessian_indices,
		std::vector< NLPTerm > &_hessian_terms) const;

	std::string to_string();


	NLPTerm operator-();

	NLPTerm& operator*=(const Number& rhs);
	friend NLPTerm operator*(NLPTerm lhs, const Number& rhs)
	{
		return lhs *= rhs;
	}

	NLPTerm& operator*=(const NLPTerm& rhs);
	friend NLPTerm operator*(NLPTerm lhs, const NLPTerm& rhs)
	{
		return lhs *= rhs;
	}

private:
	Number coeff_;
	std::list< std::pair<Index, int> > vars_;	// (index, power).
};


class NLPExpression
{
public:
	NLPExpression();
	NLPExpression(const Number _coeff, const Index _index);
	NLPExpression(const Number _coeff, const std::vector<Index> &_indices);
	NLPExpression(const NLPTerm &_term);
	NLPExpression(const NLPExpression &_expression);
	~NLPExpression();


	// FIXME:
	// Should be replaced with += operator.
	void add_term(const NLPTerm &_term);

	Number eval(const Number* _x) const;

	bool get_gradient(const Index _num_vars,
		std::vector< NLPExpression > &_gradient_expressions) const;
	bool get_hessian(const Index _num_vars,
		std::vector< std::vector<NLPExpression > > &_hessian_expressions) const;

	bool get_sparse_gradient(const Index _num_vars,
		std::vector< Index > &_gradient_indices,
		std::vector< NLPExpression > &_gradient_expressions) const;
	bool get_sparse_hessian(const Index _num_vars,
		std::vector< std::pair<Index, Index> > &_hessian_indices,
		std::vector< NLPExpression > &_hessian_expressions) const;

	std::string to_string();


	NLPExpression operator-();

	NLPExpression& operator+=(const Number& rhs);
	friend NLPExpression operator+(NLPExpression lhs, const Number& rhs)
	{
		return lhs += rhs;
	}

	NLPExpression& operator+=(const NLPTerm& rhs);
	friend NLPExpression operator+(NLPExpression lhs, const NLPTerm& rhs)
	{
		return lhs += rhs;
	}

	NLPExpression& operator+=(const NLPExpression& rhs);
	friend NLPExpression operator+(NLPExpression lhs, const NLPExpression& rhs)
	{
		return lhs += rhs;
	}


	NLPExpression& operator-=(const Number& rhs);
	friend NLPExpression operator-(NLPExpression lhs, const Number& rhs)
	{
		return lhs -= rhs;
	}

	NLPExpression& operator-=(const NLPTerm& rhs);
	friend NLPExpression operator-(NLPExpression lhs, const NLPTerm& rhs)
	{
		return lhs -= rhs;
	}

	NLPExpression& operator-=(const NLPExpression& rhs);
	friend NLPExpression operator-(NLPExpression lhs, const NLPExpression& rhs)
	{
		return lhs -= rhs;
	}


	NLPExpression& operator*=(const Number& rhs);
	friend NLPExpression operator*(NLPExpression lhs, const Number& rhs)
	{
		return lhs *= rhs;
	}

	NLPExpression& operator*=(const NLPTerm& rhs);
	friend NLPExpression operator*(NLPExpression lhs, const NLPTerm& rhs)
	{
		return lhs *= rhs;
	}

	NLPExpression& operator*=(const NLPExpression& rhs);
	friend NLPExpression operator*(NLPExpression lhs, const NLPExpression& rhs)
	{
		return lhs *= rhs;
	}

private:
	std::list< NLPTerm > terms_;
};

typedef std::vector< NLPExpression > NLPGradientExpression;
typedef std::vector< std::vector< NLPExpression > > NLPHessianExpression;

typedef struct NLPSparseGradientExpression
{
	std::vector< Index > indices_;
	std::vector< NLPExpression > expressions_;
};

typedef struct NLPSparseHessianExpression
{
	std::vector< std::pair<Index, Index> > indices_;
	std::vector< NLPExpression > expressions_;
};

#endif	// __NLP_EXPRESSION_H__
