#ifndef __VECTOR_OPERATIONS_CPP__
#define __VECTOR_OPERATIONS_CPP__

#include <cmath>
#include <vector>

double norm_2(const std::vector<double> &vec)
{
	double n2 = 0;
	for (auto elem : vec)
	{
		n2 += elem * elem;
	}
	return sqrt(n2);
}

double max_abs(const std::vector<double> &vec)
{
	double max_elem = vec[0];
	for (auto elem : vec)
	{
		elem < max_elem ? max_elem = elem : max_elem = max_elem;
	}
	return max_elem;
}

inline std::vector<double> vector_difference(const std::vector<double> &a, const std::vector<double> &b)
{
	std::vector<double> result;
	for (size_t i = 0; i < a.size(); i++)
	{
		result.push_back(a[i] - b[i]);
	}
	return result;
}

#endif
