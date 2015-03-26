#include "NLPVectorExpression.h"

#include <cassert>

NLPExpression NLPVectorExpression::dot_product(const Index _dimension,
	const Index _vec_var_start_index_1, const Index _vec_var_start_index_2)
{
	NLPExpression output;

	for (int i = 0; i < _dimension; ++i)
	{
		std::vector < Index > indices;
		indices.push_back(_vec_var_start_index_1 + i);
		indices.push_back(_vec_var_start_index_2 + i);
		output.add_term(NLPTerm(1.0, indices));
	}

	return output;
}

void NLPVectorExpression::cross_product(const Index _dimension,
	const Index _vec_var_start_index_1, const Index _vec_var_start_index_2,
	std::vector<NLPExpression>& _output)
{
	// NOTICE:
	// Implemented only for 3 dimension.
	assert(_dimension == 3);
	_output.clear();
	_output.resize(_dimension);

	// (coeff, v1_index, v2_index);
	const int cross_product_terms[3][2][3] = {
		{ { 1, 0, 1 }, { -1, 1, 0 } },
		{ { 1, 1, 2 }, { -1, 2, 1 } },
		{ { 1, 2, 0 }, { -1, 0, 2 } }
	};

	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			Number coeff = static_cast<Number>(cross_product_terms[i][j][0]);
			std::vector < Index > indices;
			indices.push_back(_vec_var_start_index_1 + cross_product_terms[i][j][1]);
			indices.push_back(_vec_var_start_index_2 + cross_product_terms[i][j][2]);
			_output[i].add_term(NLPTerm(coeff, indices));
		}
	}
}
