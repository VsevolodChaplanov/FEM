#ifndef __VECTOR_OPERATIONS_CPP__
#define __VECTOR_OPERATIONS_CPP__

#include "VecLib.h"

std::vector<double> NumProduct(const std::vector<double> &vector, const double a)
{
	std::size_t N = vector.size();
	std::vector<double> result(N);
	for(std::size_t i = 0; i < N; i++)
	{
		result[i] = vector[i] * a;
	}
	return result;
}

double DotProduct(const std::vector<double> &a, const std::vector<double> &b)
{
	std::size_t N = a.size();
	double result;
	for(std::size_t i = 0; i < N; i++)
	{
		result += a[i] * b[i];
	}
	return result;
}

#endif