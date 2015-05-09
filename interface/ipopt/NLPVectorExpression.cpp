#include "NLPVectorExpression.h"

#include <cassert>
#include <iostream>


NLPVectorExpression::NLPVectorExpression(const Index _dimension)
{
	expressions_.resize(_dimension, NLPExpression());
}

NLPVectorExpression::~NLPVectorExpression()
{

}

void NLPVectorExpression::get_segment(const Index _index, const Index _sub_dimension) const
{
	assert(_index + _sub_dimension <= dimension());

	NLPVectorExpression segment(_sub_dimension);
	for (Index i = 0; i < _sub_dimension; ++i)
	{
		segment.expressions_[i] = this->expressions_[_index + i];
	}
}

void NLPVectorExpression::set_segment(const Index _index, const NLPVectorExpression &_segment)
{
	const Index sub_dimension = _segment.dimension();
	assert(_index + sub_dimension <= dimension());

	for (Index i = 0; i < sub_dimension; ++i)
	{
		this->expressions_[i] = _segment.expressions_[_index + i];
	}
}

NLPVectorExpression NLPVectorExpression::operator-()
{
	NLPVectorExpression other(*this);
	for (std::vector< NLPExpression >::iterator it = other.expressions_.begin();
		it != other.expressions_.end(); ++it)
	{
		(*it) *= (-1);
	}
	return other;
}

NLPVectorExpression& NLPVectorExpression::operator+=(const NLPVectorExpression& rhs)
{
	assert(this->dimension() == rhs.dimension());
	std::vector< NLPExpression >::iterator it = this->expressions_.begin();
	std::vector< NLPExpression >::const_iterator jt = rhs.expressions_.begin();
	for (; it != this->expressions_.end() && jt != rhs.expressions_.end(); ++it, ++jt)
	{
		(*it) += (*jt);
	}
	return *this;
}

NLPVectorExpression& NLPVectorExpression::operator+=(const Eigen::VectorXd& rhs)
{
	assert(this->dimension() == rhs.rows());
	for (Index i = 0; i < dimension(); ++i)
	{
		expressions_[i] += static_cast<Number>(rhs[i]);
	}
	return *this;
}

NLPVectorExpression& NLPVectorExpression::operator-=(const NLPVectorExpression& rhs)
{
	assert(this->dimension() == rhs.dimension());
	std::vector< NLPExpression >::iterator it = this->expressions_.begin();
	std::vector< NLPExpression >::const_iterator jt = rhs.expressions_.begin();
	for (; it != this->expressions_.end() && jt != rhs.expressions_.end(); ++it, ++jt)
	{
		(*it) -= (*jt);
	}
	return *this;
}

NLPVectorExpression& NLPVectorExpression::operator-=(const Eigen::VectorXd& rhs)
{
	assert(this->dimension() == rhs.rows());
	for (Index i = 0; i < dimension(); ++i)
	{
		expressions_[i] -= static_cast<Number>(rhs[i]);
	}
	return *this;
}

NLPVectorExpression& NLPVectorExpression::operator*=(const Number& rhs)
{
	for (std::vector< NLPExpression >::iterator it = this->expressions_.begin();
		it != this->expressions_.end(); ++it)
	{
		(*it) *= rhs;
	}
	return *this;
}

NLPVectorExpression& NLPVectorExpression::operator*=(const NLPExpression& rhs)
{
	for (std::vector< NLPExpression >::iterator it = this->expressions_.begin();
		it != this->expressions_.end(); ++it)
	{
		(*it) *= rhs;
	}
	return *this;
}

NLPExpression& NLPVectorExpression::operator[](const Index i)
{
	assert(i < dimension());
	return expressions_[i];
}

const NLPExpression& NLPVectorExpression::operator[](const Index i) const
{
	assert(i < dimension());
	return expressions_[i];
}

NLPExpression NLPVectorExpression::dot_product(
	const NLPVectorExpression& _vector_1, const NLPVectorExpression& _vector_2)
{
	assert(_vector_1.dimension() == _vector_2.dimension());
	const Index dimension = _vector_1.dimension();

	NLPExpression output;
	for (unsigned int i = 0; i < dimension; ++i)
		output += (_vector_1.expressions_[i] * _vector_2.expressions_[i]);

	return output;
}

NLPExpression NLPVectorExpression::dot_product(
	const NLPVectorExpression& _vector, const Eigen::VectorXd& _constant_vector)
{
	assert(_vector.dimension() == _constant_vector.rows());
	const Index dimension = _vector.dimension();

	NLPExpression output;
	for (unsigned int i = 0; i < dimension; ++i)
		output += (_vector.expressions_[i] * _constant_vector[i]);

	return output;
}

NLPVectorExpression NLPVectorExpression::cross_product(
	const NLPVectorExpression& _vector_1, const NLPVectorExpression& _vector_2)
{
	assert(_vector_1.dimension() == _vector_2.dimension());
	const Index dimension = _vector_1.dimension();

	// NOTE:
	// Implemented only for 3 dimension.
	if (dimension != 3)
	{
		do {
			std::cout << "Error: Cross product is implemented only for 3-dimension." << std::endl;
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	NLPVectorExpression output(dimension);

	output[0] = _vector_1.expressions_[1] * _vector_2.expressions_[2]
		- _vector_1.expressions_[2] * _vector_2.expressions_[1];

	output[1] = _vector_1.expressions_[2] * _vector_2.expressions_[0]
		- _vector_1.expressions_[0] * _vector_2.expressions_[2];

	output[2] = _vector_1.expressions_[0] * _vector_2.expressions_[1]
		- _vector_1.expressions_[1] * _vector_2.expressions_[0];

	return output;
}

NLPVectorExpression NLPVectorExpression::cross_product(
	const NLPVectorExpression& _vector, const Eigen::VectorXd& _scalar_vector)
{
	assert(_vector.dimension() == _scalar_vector.rows());
	const Index dimension = _vector.dimension();

	// NOTE:
	// Implemented only for 3 dimension.
	if (dimension != 3)
	{
		do {
			std::cout << "Error: Cross product is implemented only for 3-dimension." << std::endl;
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}

	NLPVectorExpression output(dimension);

	output[0] = _vector.expressions_[1] * _scalar_vector[2]
		- _vector.expressions_[2] * _scalar_vector[1];

	output[1] = _vector.expressions_[2] * _scalar_vector[0]
		- _vector.expressions_[0] * _scalar_vector[2];

	output[2] = _vector.expressions_[0] * _scalar_vector[1]
		- _vector.expressions_[1] * _scalar_vector[0];

	return output;
}
