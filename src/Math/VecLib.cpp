#ifndef __VECTOR_OPERATIONS_CPP__
#define __VECTOR_OPERATIONS_CPP__

#include "VecLib.h"
#include <vector>

template<typename T>
std::vector<T> NumProduct(const std::vector<T> &vector, const double a)
{
	std::size_t N = vector.size();
	std::vector<T> result(N);
	for(std::size_t i = 0; i < N; i++)
	{
		result[i] = vector[i] * a;
	}
	return result;
}

template<typename T>
T DotProduct(const std::vector<T> &a, const std::vector<T> &b)
{
	std::size_t N = a.size();
	T result;
	for(std::size_t i = 0; i < N; i++)
	{
		result += a[i] * b[i];
	}
	return result;
}

#endif