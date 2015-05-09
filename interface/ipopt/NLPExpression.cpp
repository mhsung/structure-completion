#include "NLPExpression.h"

#include <stdio.h>
#include <string.h>
#include <cassert>
#include <iostream>
#include <sstream>

NLPTerm::NLPTerm()
	: coeff_(0)
{

}

NLPTerm::NLPTerm(const Number _coeff)
{
	set_coeff(_coeff);
}

NLPTerm::NLPTerm(const Number _coeff, const Index _index)
{
	set_coeff(_coeff);
	multiply_var(_index, 1);
}

NLPTerm::NLPTerm(const Number _coeff, const std::vector<Index> &_indices)
{
	set_coeff(_coeff);
	for (std::vector<Index>::const_iterator index_it = _indices.begin();
		index_it != _indices.end(); ++index_it)
		multiply_var(*index_it, 1);
}

NLPTerm::NLPTerm(const NLPTerm &_term)
	: coeff_(_term.coeff_)
	, vars_(_term.vars_)
{

}

NLPTerm::~NLPTerm()
{

}

bool NLPTerm::empty() const
{
	return (coeff_ == 0);
}

void NLPTerm::set_coeff(Number _coeff)
{
	coeff_ = _coeff;
}

void NLPTerm::multiply_coeff(Number _coeff)
{
	coeff_ *= _coeff;
}

void NLPTerm::multiply_var(Index _index, int _power)
{
	if (_power <= 0)
		return;

	bool exist = false;
	for (std::list< std::pair<Index, int> >::iterator var_it = vars_.begin();
		var_it != vars_.end(); ++var_it)
	{
		if ((*var_it).first == _index)
		{
			(*var_it).second += _power;
			exist = true;
			break;
		}
	}

	if (!exist)
	{
		vars_.push_back(std::make_pair(_index, _power));
	}
}

Number NLPTerm::eval(const Number* _x) const
{
	Number ret = coeff_;
	for (std::list< std::pair< Index, int > >::const_iterator var_it = vars_.begin();
		var_it != vars_.end(); ++var_it)
	{
		Index index = (*var_it).first;
		int power = (*var_it).second;
		for (int count = 0; count < power; ++count)
			ret *= _x[index];
	}
	return ret;
}

bool NLPTerm::get_gradient(const Index _num_vars, std::vector<NLPTerm> &_gradient_terms) const
{
	_gradient_terms.clear();
	_gradient_terms.resize(_num_vars);

	for (std::list< std::pair< Index, int > >::const_iterator i_var_it = vars_.begin();
		i_var_it != vars_.end(); ++i_var_it)
	{
		Index index_i = (*i_var_it).first;
		if (index_i >= _num_vars)
		{
			std::cerr << "Error: A variable index is equal or greater than the number of variables: ("
				<< index_i << " >= " << _num_vars << ")" << std::endl;
			return false;
		}

		int power_i = (*i_var_it).second;
		if (power_i <= 0) continue;

		// Copy coefficient.
		_gradient_terms[index_i].set_coeff(coeff_);

		// d(c*x^n)/dx = c*n*x^(n-1).
		_gradient_terms[index_i].multiply_coeff(power_i);
		_gradient_terms[index_i].multiply_var(index_i, power_i - 1);

		for (std::list< std::pair< Index, int > >::const_iterator j_var_it = vars_.begin();
			j_var_it != vars_.end(); ++j_var_it)
		{
			if (i_var_it == j_var_it) continue;

			Index index_j = (*j_var_it).first;
			if (index_j >= _num_vars)
			{
				std::cerr << "Error: A variable index is equal or greater than the number of variables: ("
					<< index_j << " >= " << _num_vars << ")" << std::endl;
				return false;
			}

			int power_j = (*j_var_it).second;
			if (power_j <= 0) continue;

			// d(c*x^n)/dy = c*x^n.
			_gradient_terms[index_i].multiply_var(index_j, power_j);
		}
	}

	return true;
}

bool NLPTerm::get_hessian(const Index _num_vars,
	std::vector< std::vector< NLPTerm > > &_hessian_terms) const
{
	bool ret;

	std::vector<NLPTerm> gradient_terms;
	ret = get_gradient(_num_vars, gradient_terms);
	if (!ret) return false;
	assert(gradient_terms.size() == _num_vars);

	_hessian_terms.clear();
	_hessian_terms.resize(_num_vars);
	for (Index index = 0; index < _num_vars; ++index)
	{
		ret = gradient_terms[index].get_gradient(_num_vars, _hessian_terms[index]);
		if (!ret) return false;
	}

	return true;
}

bool NLPTerm::get_sparse_gradient(const Index _num_vars, std::vector< Index > &_gradient_indices,
	std::vector< NLPTerm > &_gradient_terms) const
{
	_gradient_indices.clear();
	_gradient_terms.clear();

	_gradient_indices.reserve(_num_vars);
	_gradient_terms.reserve(_num_vars);

	for (std::list< std::pair< Index, int > >::const_iterator var_it_i = vars_.begin();
		var_it_i != vars_.end(); ++var_it_i)
	{
		Index index_i = (*var_it_i).first;
		if (index_i >= _num_vars)
		{
			std::cerr << "Error: A variable index is equal or greater than the number of variables: ("
				<< index_i << " >= " << _num_vars << ")" << std::endl;
			return false;
		}

		int power_i = (*var_it_i).second;
		if (power_i <= 0) continue;

		NLPTerm term_i;

		// Copy coefficient.
		term_i.set_coeff(coeff_);

		// d(c*x^n)/dx = c*n*x^(n-1).
		term_i.multiply_coeff(power_i);
		term_i.multiply_var(index_i, power_i - 1);

		for (std::list< std::pair< Index, int > >::const_iterator var_it_j = vars_.begin();
			var_it_j != vars_.end(); ++var_it_j)
		{
			if (var_it_i == var_it_j) continue;

			Index index_j = (*var_it_j).first;
			if (index_j >= _num_vars)
			{
				std::cerr << "Error: A variable index is equal or greater than the number of variables: ("
					<< index_j << " >= " << _num_vars << ")" << std::endl;
				return false;
			}

			int power_j = (*var_it_j).second;
			if (power_j <= 0) continue;

			// d(c*x^n)/dy = c*x^n.
			term_i.multiply_var(index_j, power_j);
		}

		if (!term_i.empty())
		{
			_gradient_indices.push_back(index_i);
			_gradient_terms.push_back(term_i);
		}
	}

	return true;
}

bool NLPTerm::get_sparse_hessian(const Index _num_vars,
	std::vector< std::pair<Index, Index> > &_hessian_indices,
	std::vector< NLPTerm > &_hessian_terms) const
{
	bool ret;

	std::vector<Index> gradient_indices_i;
	std::vector<NLPTerm> gradient_terms_i;
	ret = get_sparse_gradient(_num_vars, gradient_indices_i, gradient_terms_i);
	if (!ret) return false;

	unsigned int num_gradient_termss_i = gradient_terms_i.size();
	assert(gradient_indices_i.size() == num_gradient_termss_i);


	_hessian_indices.clear();
	_hessian_terms.clear();

	int max_num_pairs = _num_vars*(_num_vars + 1) / 2;
	_hessian_indices.reserve(max_num_pairs);
	_hessian_terms.reserve(max_num_pairs);

	for (unsigned int i = 0; i < num_gradient_termss_i; ++i)
	{
		std::vector<Index> gradient_indices_j;
		std::vector<NLPTerm> gradient_terms_j;

		ret = gradient_terms_i[i].get_sparse_gradient(_num_vars, gradient_indices_j, gradient_terms_j);
		if (!ret) return false;

		unsigned int num_gradient_termss_j = gradient_terms_j.size();
		assert(gradient_indices_j.size() == num_gradient_termss_j);

		Index index_i = gradient_indices_i[i];
		for (unsigned int j = 0; j < num_gradient_termss_j; ++j)
		{
			Index index_j = gradient_indices_j[j];
			NLPTerm term = gradient_terms_j[j];
			if (index_i >= index_j)
			{
				_hessian_indices.push_back(std::make_pair(index_i, index_j));
				_hessian_terms.push_back(term);
			}
		}
	}

	return true;
}

std::string NLPTerm::to_string()
{
	std::stringstream sstr;

	if (coeff_ != 1 || vars_.empty())
		sstr << coeff_ << " ";

	for (std::list< std::pair<Index, int> >::const_iterator var_it = vars_.begin();
		var_it != vars_.end(); ++var_it)
	{
		Index index = (*var_it).first;
		int power = (*var_it).second;
		if (power <= 0) continue;

		sstr << "x" << index;
		if (power >= 2) sstr << "^" << power;
		sstr << " ";
	}

	return sstr.str();
}

NLPTerm NLPTerm::operator-()
{
	NLPTerm other(*this);
	other.coeff_ *= -1;
	return other;
}

NLPTerm& NLPTerm::operator*=(const Number& rhs)
{
	coeff_ *= rhs;
	return *this;
}

NLPTerm& NLPTerm::operator*=(const NLPTerm& rhs)
{
	this->multiply_coeff(rhs.coeff_);
	for (std::list< std::pair<Index, int> >::const_iterator it = rhs.vars_.begin();
		it != rhs.vars_.end(); ++it)
	{
		this->multiply_var((*it).first, (*it).second);
	}
	return *this;
}

NLPExpression::NLPExpression()
{

}

NLPExpression::NLPExpression(const Number _coeff, const Index _index)
{
	add_term(NLPTerm(_coeff, _index));
}

NLPExpression::NLPExpression(const Number _coeff, const std::vector<Index> &_indices)
{
	add_term(NLPTerm(_coeff, _indices));
}

NLPExpression::NLPExpression(const NLPTerm &_term)
{
	add_term(_term);
}

NLPExpression::NLPExpression(const NLPExpression &_expression)
	: terms_(_expression.terms_)
{

}

NLPExpression::~NLPExpression()
{

}

void NLPExpression::add_term(const NLPTerm &_term)
{
	if (!_term.empty())
		terms_.push_back(_term);
}

Number NLPExpression::eval(const Number* _x) const
{
	Number ret = 0;
	for (std::list< NLPTerm >::const_iterator it = terms_.begin(); it != terms_.end(); ++it)
	{
		ret += (*it).eval(_x);
	}
	return ret;
}

bool NLPExpression::get_gradient(const Index _num_vars,
	std::vector< NLPExpression > &_gradient_expressions) const
{
	bool ret;

	_gradient_expressions.clear();
	_gradient_expressions.resize(_num_vars, NLPExpression());

	for (std::list< NLPTerm >::const_iterator it = terms_.begin(); it != terms_.end(); ++it)
	{
		std::vector< NLPTerm > gradient_terms;
		ret = (*it).get_gradient(_num_vars, gradient_terms);
		if (!ret) return false;

		for (Index index = 0; index < _num_vars; ++index)
			_gradient_expressions[index].add_term(gradient_terms[index]);
	}

	return true;
}

bool NLPExpression::get_hessian(const Index _num_vars,
	std::vector< std::vector< NLPExpression > > &_hessian_expressions) const
{
	bool ret;

	_hessian_expressions.clear();
	_hessian_expressions.resize(_num_vars);
	for (Index index = 0; index < _num_vars; ++index)
		_hessian_expressions[index].resize(_num_vars, NLPExpression());

	for (std::list< NLPTerm >::const_iterator it = terms_.begin(); it != terms_.end(); ++it)
	{
		std::vector< std::vector< NLPTerm > > hessian_term;
		ret = (*it).get_hessian(_num_vars, hessian_term);
		if (!ret) return false;

		for (Index index_i = 0; index_i < _num_vars; ++index_i)
			for (Index index_j = 0; index_j < _num_vars; ++index_j)
				_hessian_expressions[index_i][index_j].add_term(hessian_term[index_i][index_j]);
	}

	return true;
}

bool NLPExpression::get_sparse_gradient(const Index _num_vars,
	std::vector< Index > &_gradient_indices,
	std::vector< NLPExpression > &_gradient_expressions) const
{
	bool ret;

	_gradient_indices.clear();
	_gradient_expressions.clear();

	_gradient_indices.reserve(_num_vars);
	_gradient_expressions.reserve(_num_vars);

	for (std::list< NLPTerm >::const_iterator it = terms_.begin(); it != terms_.end(); ++it)
	{
		std::vector< Index > gradient_terms_indices;
		std::vector< NLPTerm > gradient_terms;
		ret = (*it).get_sparse_gradient(_num_vars, gradient_terms_indices, gradient_terms);
		if (!ret) return false;

		unsigned int num_gradient_terms = gradient_terms.size();
		assert(gradient_terms_indices.size() == num_gradient_terms);

		for (unsigned int i = 0; i < num_gradient_terms; ++i)
		{
			Index index_i = gradient_terms_indices[i];
			NLPTerm term = gradient_terms[i];

			assert(_gradient_indices.size() == _gradient_expressions.size());
			bool exist = false;
			for (Index k = 0; k < _gradient_indices.size(); ++k)
			{
				if (index_i == _gradient_indices[k])
				{
					_gradient_expressions[k].add_term(term);
					exist = true;
					break;
				}
			}

			if (!exist)
			{
				_gradient_indices.push_back(index_i);
				_gradient_expressions.push_back(NLPExpression());
				_gradient_expressions.back().add_term(term);
			}
		}
	}

	return true;
}

bool NLPExpression::get_sparse_hessian(const Index _num_vars,
	std::vector< std::pair<Index, Index> > &_hessian_indices,
	std::vector< NLPExpression > &_hessian_expressions) const
{
	bool ret;

	_hessian_indices.clear();
	_hessian_expressions.clear();

	int max_num_pairs = _num_vars*(_num_vars + 1) / 2;
	_hessian_indices.reserve(max_num_pairs);
	_hessian_expressions.reserve(max_num_pairs);

	for (std::list< NLPTerm >::const_iterator it = terms_.begin(); it != terms_.end(); ++it)
	{
		std::vector< std::pair<Index, Index> > hessian_term_indices;
		std::vector< NLPTerm > hessian_term;
		ret = (*it).get_sparse_hessian(_num_vars, hessian_term_indices, hessian_term);
		if (!ret) return false;

		unsigned int num_hessian_termss = hessian_term.size();
		assert(hessian_term_indices.size() == num_hessian_termss);

		for (unsigned int i = 0; i < num_hessian_termss; ++i)
		{
			std::pair<Index, Index> index_i = hessian_term_indices[i];
			NLPTerm term = hessian_term[i];

			assert(_hessian_indices.size() == _hessian_expressions.size());
			bool exist = false;
			for (Index k = 0; k < _hessian_indices.size(); ++k)
			{
				if (index_i.first == _hessian_indices[k].first
					&& index_i.second == _hessian_indices[k].second)
				{
					_hessian_expressions[k].add_term(term);
					exist = true;
					break;
				}
			}

			if (!exist)
			{
				_hessian_indices.push_back(index_i);
				_hessian_expressions.push_back(NLPExpression());
				_hessian_expressions.back().add_term(term);
			}
		}
	}

	return true;
}

std::string NLPExpression::to_string()
{
	std::stringstream sstr;

	for (std::list< NLPTerm >::iterator it = terms_.begin(); it != terms_.end(); ++it)
	{
		if (it != terms_.begin())
			sstr << "+ ";
		sstr << (*it).to_string();
	}

	return sstr.str();
}

NLPExpression NLPExpression::operator-()
{
	NLPExpression other(*this);
	for (std::list< NLPTerm >::iterator it = other.terms_.begin();
		it != other.terms_.end(); ++it)
	{
		(*it) *= (-1);
	}
	return other;
}

NLPExpression& NLPExpression::operator+=(const Number& rhs)
{
	terms_.push_back(NLPTerm(rhs));
	return *this;
}

NLPExpression& NLPExpression::operator+=(const NLPTerm& rhs)
{
	if (!rhs.empty())
		terms_.push_back(rhs);
	return *this;
}

NLPExpression& NLPExpression::operator+=(const NLPExpression& rhs)
{
	terms_.insert(terms_.end(), rhs.terms_.begin(), rhs.terms_.end());
	return *this;
}

NLPExpression& NLPExpression::operator-=(const Number& rhs)
{
	terms_.push_back(NLPTerm(-rhs));
	return *this;
}

NLPExpression& NLPExpression::operator-=(const NLPTerm& rhs)
{
	if (!rhs.empty())
	{
		NLPTerm term(rhs);
		term *= (-1);
		terms_.push_back(term);
	}
	return *this;
}

NLPExpression& NLPExpression::operator-=(const NLPExpression& rhs)
{
	for (std::list< NLPTerm >::const_iterator it = rhs.terms_.begin();
		it != rhs.terms_.end(); ++it)
	{
		NLPTerm term(*it);
		term *= (-1);
		terms_.push_back(term);
	}
	return *this;
}

NLPExpression& NLPExpression::operator*=(const Number& rhs)
{
	for (std::list< NLPTerm >::iterator it = terms_.begin();
		it != terms_.end(); ++it)
	{
		(*it) *= rhs;
	}
	return *this;
}

NLPExpression& NLPExpression::operator*=(const NLPTerm& rhs)
{
	for (std::list< NLPTerm >::iterator it = terms_.begin();
		it != terms_.end(); ++it)
	{
		(*it) *= rhs;
	}
	return *this;
}

NLPExpression& NLPExpression::operator*=(const NLPExpression& rhs)
{
	std::list< NLPTerm > product_terms;

	for (std::list< NLPTerm >::iterator it = terms_.begin();
		it != terms_.end(); ++it)
	{
		for (std::list< NLPTerm >::const_iterator jt = rhs.terms_.begin();
			jt != rhs.terms_.end(); ++jt)
		{
			product_terms.push_back((*it) * (*jt));
		}
	}

	product_terms.swap(terms_);

	return *this;
}
