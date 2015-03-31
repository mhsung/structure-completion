// Author: Minhyuk Sung
// Mar. 2015

#ifndef __NLP_VECTOR_EXPRESSION_H__
#define __NLP_VECTOR_EXPRESSION_H__

#include "NLPExpression.h"

#include "NLPExpression.h"


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

	static NLPVectorExpression cross_product(const NLPVectorExpression& _vector_1,
		const NLPVectorExpression& _vector_2);


	NLPVectorExpression& operator+=(const NLPVectorExpression& rhs);
	friend NLPVectorExpression operator+(NLPVectorExpression lhs, const NLPVectorExpression& rhs)
	{
		return lhs += rhs;
	}

	NLPVectorExpression& operator-=(const NLPVectorExpression& rhs);
	friend NLPVectorExpression operator-(NLPVectorExpression lhs, const NLPVectorExpression& rhs)
	{
		return lhs -= rhs;
	}

	NLPVectorExpression& operator*=(const Number& rhs);
	friend NLPVectorExpression operator*(NLPVectorExpression lhs, const Number& rhs)
	{
		return lhs *= rhs;
	}

	NLPVectorExpression& operator*=(const NLPExpression& rhs);
	friend NLPVectorExpression operator*(NLPVectorExpression lhs, const NLPExpression& rhs)
	{
		return lhs *= rhs;
	}

	NLPExpression& operator[](const Index i);

	const NLPExpression& operator[](const Index i) const;


private:
	std::vector< NLPExpression > expressions_;
};

#endif	// __NLP_VECTOR_EXPRESSION_H__
