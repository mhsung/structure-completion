// Author: Minhyuk Sung
// Mar. 2015

#ifndef __NLP_VECTOR_EXPRESSION_H__
#define __NLP_VECTOR_EXPRESSION_H__

#include "NLPExpression.h"


class NLPVectorExpression
{
public:
	static NLPExpression dot_product(const Index _dimension,
		const Index _vec_var_start_index_1, const Index _vec_var_start_index_2);

	static void cross_product(const Index _dimension,
		const Index _vec_var_start_index_1, const Index _vec_var_start_index_2,
		std::vector<NLPExpression>& _output);
};

#endif	// __NLP_VECTOR_EXPRESSION_H__
