#ifndef __LINELEMCPP__
#define __LINELEMCPP__

#include "LinElem.h"
#include <cmath>

LinElem::LinElem(const std::vector<double> &verteces, 
				 const std::vector<std::size_t> &GIdences) : GlobalIndices(GIdences), StartPoint(verteces[0]),
				 Length(abs(verteces[1] - verteces[0])) { }

double LinElem::phi1(const double &_xi)
{
	return 1 - _xi;
}

double LinElem::phi2(const double &_xi)
{
	return _xi;
}

#endif