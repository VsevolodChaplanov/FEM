#include "Builder1D.h"

Builder1D::Builder1D(const double &left_bound, const double &right_bound, const std::size_t &NN)
{
	double h = (right_bound - left_bound) / (double) NN;
	std::vector<double> mesh(NN+1);
	for (size_t i = 0; i < NN+1; i++)
	{
		mesh[i] = left_bound + i * h;
	}
	for (size_t i = 0; i < NN; i++)
	{
		LinElem* LinearElement = new LinElem({mesh[i], mesh[i+1]}, {i,i+1});
		// TODO ускорение
		// 			Задать сразу длинну вектора по NN
		Elements.push_back(LinearElement);
	}
}