// Author: Minhyuk Sung
// Mar. 2015

#ifndef __NLP_VECTOR_EXPRESSION_H__
#define __NLP_VECTOR_EXPRESSION_H__

#include "NLPExpression.h"

#include <Eigen/Core>


class NLPVectorExpression
{
public:
	NLPVectorExpression(const Index _dimension);
	~NLPVectorExpression();

	inline Index dimension() const { return expressions_.size(); }
	void get_segment(const Index _index, const Index _sub_dimension) const;
	void set_segment(const Index _index, const NLPVectorExpression &_segment);

	static NLPExpression dot_product(const NLPVectorExpression& _vector_1,
		const NLPVectorExpression& _vector_2);

	static NLPExpression dot_product(const NLPVectorExpression& _vector,
		const Eigen::VectorXd& _constant_vector);

	static NLPExpression dot_product(const Eigen::VectorXd& _constant_vector,
		const NLPVectorExpression& _vector)
	{
		return dot_product(_vector, _constant_vector);
	}

	static NLPVectorExpression cross_product(const NLPVectorExpression& _vector_1,
		const NLPVectorExpression& _vector_2);

	static NLPVectorExpression cross_product(const NLPVectorExpression& _vector,
		const Eigen::VectorXd& _constant_vector);

	static NLPVectorExpression cross_product(const Eigen::VectorXd& _constant_vector,
		const NLPVectorExpression& _vector)
	{
		return -cross_product(_vector, _constant_vector);
	}


	NLPVectorExpression operator-();

	NLPVectorExpression& operator+=(const NLPVectorExpression& rhs);
	friend NLPVectorExpression operator+(NLPVectorExpression lhs, const NLPVectorExpression& rhs)
	{
		return lhs += rhs;
	}

	NLPVectorExpression& operator+=(const Eigen::VectorXd& rhs);
	friend NLPVectorExpression operator+(NLPVectorExpression lhs, const Eigen::VectorXd& rhs)
	{
		return lhs += rhs;
	}

	friend NLPVectorExpression operator+(const Eigen::VectorXd& lhs, NLPVectorExpression rhs)
	{
		return (rhs) += lhs;
	}

	NLPVectorExpression& operator-=(const NLPVectorExpression& rhs);
	friend NLPVectorExpression operator-(NLPVectorExpression lhs, const NLPVectorExpression& rhs)
	{
		return lhs -= rhs;
	}

	NLPVectorExpression& operator-=(const Eigen::VectorXd& rhs);
	friend NLPVectorExpression operator-(NLPVectorExpression lhs, const Eigen::VectorXd& rhs)
	{
		return lhs -= rhs;
	}

	friend NLPVectorExpression operator-(const Eigen::VectorXd& lhs, NLPVectorExpression rhs)
	{
		return (-rhs) += lhs;
	}

	NLPVectorExpression& operator*=(const Number& rhs);
	friend NLPVectorExpression operator*(NLPVectorExpression lhs, const Number& rhs)
	{
		return lhs *= rhs;
	}

	friend NLPVectorExpression operator*(const Number& lhs, NLPVectorExpression rhs)
	{
		return (rhs) *= lhs;
	}

	NLPVectorExpression& operator*=(const NLPExpression& rhs);
	friend NLPVectorExpression operator*(NLPVectorExpression lhs, const NLPExpression& rhs)
	{
		return lhs *= rhs;
	}

	friend NLPVectorExpression operator*(const NLPExpression& lhs, NLPVectorExpression rhs)
	{
		return (rhs) *= lhs;
	}


	NLPExpression& operator[](const Index i);

	const NLPExpression& operator[](const Index i) const;


private:
	std::vector< NLPExpression > expressions_;
};

#endif	// __NLP_VECTOR_EXPRESSION_H__
