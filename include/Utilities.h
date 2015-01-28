#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <cmath>
#include <iostream>
#include <string>


#define DEBUG_TEST
#define NUMERIAL_ERROR_THRESHOLD	1.0E-6

void CHECK_NUMERICAL_ERROR(const std::string& _desc, const double& _error);
void CHECK_NUMERICAL_ERROR(const std::string& _desc, const double& _value_1, const double& _value_2);

#endif	// _UTILITIES_H_