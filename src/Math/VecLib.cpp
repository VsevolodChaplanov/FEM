#include "VecLib.h"

template<typename T>
std::vector<T> NumProduct(const std::vector<T> &vector, const double a)
{
	std::size_t N = vector.size();
	std::vector<T> result(N);
	for(std::size_t i = 0; i < N, i++)
	{
		result[i] = vector[i] * a;
	}
	return result;
}

template<typename T>
T DotProduct(const std::vector<T> &a, const std::vector<T> &b)
{
	std::size_t N = vector.size();
	T result;
	for(std::size_t i = 0; i < N, i++)
	{
		result += a[i] * b[i];
	}
	return result;
}