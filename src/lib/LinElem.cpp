#ifndef __LINELEM_CPP__
#define __LINELEM_CPP__

#include "LinElem.h"
#include <cmath>

LinElem::LinElem(const std::vector<double> &verteces, 
	const std::vector<std::size_t> &GIdences) : GlobalIndices(GIdences), StartPoint(verteces[0]),
	Length(abs(verteces[1] - verteces[0]))
{ }

double LinElem::phi1(const double &_xi)
{
	return 1 - _xi;
}

double LinElem::phi2(const double &_xi)
{
	return _xi;
}

std::vector<std::vector<double>>* LinElem::GetMass(const std::size_t i, const std::size_t j)
{
	return &MassMatrix;
}

std::vector<std::vector<double>>* LinElem::GetStiffness(const std::size_t i, const std::size_t j)
{
	return &StiffnessMatrix;
}

std::vector<double>* LinElem::GetLumped(const std::size_t i, const std::size_t j)
{
	return &LumpedMassMatrix;
}

#endif