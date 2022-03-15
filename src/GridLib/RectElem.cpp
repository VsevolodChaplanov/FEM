#ifndef __RECTANGULAR_ELEMENT_CPP__
#define __RECTANGULAR_ELEMENT_CPP__

#include "RectElem.h"

double RectElem::phi1(const double &_xi, const double &_eta)
{
	return (1-_xi)*(1-_eta);
}

double RectElem::phi2(const double &_xi, const double &_eta)
{
	return _xi*(1-_eta);
}

double RectElem::phi3(const double &_xi, const double &_eta)
{
	return _xi*_eta;
}

double RectElem::phi4(const double &_xi, const double &_eta)
{
	return (1-_xi)*_eta;
}

std::vector<std::vector<double>>* RectElem::GetMass(const std::size_t i, const std::size_t j)
{
	return &MassMatrix;
}

std::vector<std::vector<double>>* RectElem::GetStiffness(const std::size_t i, const std::size_t j)
{
	return &StiffnessMatrix;
}

std::vector<double>* RectElem::GetLumped(const std::size_t i, const std::size_t j)
{
	return &LumpedMatrix;
}

#endif