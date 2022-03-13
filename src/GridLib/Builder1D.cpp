#ifndef __BUILDER_1D_CPP__
#define __BUILDER_1D_CPP__

#include "Builder1D.h"

Builder1D::Builder1D(const double &left_bound, const double &right_bound, const std::size_t &NN) : nn(NN)
{
	elems.reserve(nn);
	double h = (right_bound - left_bound) / (double) nn;
	std::vector<double> mesh(nn+1);
	for (size_t i = 0; i < nn+1; i++)
	{
		mesh[i] = left_bound + i * h;
	}
	for (size_t i = 0; i < nn; i++)
	{
		LinElem* LinearElement = new LinElem({mesh[i], mesh[i+1]}, {i,i+1});
		// TODO ускорение
		// 			Задать сразу длинну вектора по NN
		elems[i] = LinearElement;
		//Elements.push_back(LinearElement);
	}
}

Builder1D::Builder1D(const std::vector<double> &Mesh) : nn(Mesh.size() - 1)
{
	elems.reserve(nn);
	for (size_t i = 0; i < nn; i++)
	{
		LinElem* LinearElement = new LinElem({Mesh[i], Mesh[i+1]}, {i,i+1});
		elems[i] = LinearElement;
	}
}

std::size_t Builder1D::GetNElem()
{
	return nn;
}

#endif