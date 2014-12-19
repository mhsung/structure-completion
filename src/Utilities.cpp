#include "Utilities.h"

void CHECK_NUMERICAL_ERROR(const std::string& _desc, const double& _error)
{
	if (std::abs(_error) > NUMERIAL_ERROR_THRESHOLD)
	{
		std::cerr << "(" << _desc << "): Numerical Error (" << _error << ")" << std::endl;
		do {
			std::cout << '\n' << "Press the Enter key to continue.";
		} while (std::cin.get() != '\n');
	}
}

void CHECK_NUMERICAL_ERROR(const std::string& _desc, const double& _value_1, const double& _value_2)
{
	CHECK_NUMERICAL_ERROR(_desc, _value_1 - _value_2);
}