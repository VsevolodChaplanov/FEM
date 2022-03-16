#ifndef __FINITEELEMENTSMESHBUILDER__
#define __FINITEELEMENTSMESHBUILDER__

#include <vector>
#include <cmath>
#include "FemGrid.h"
#include "IFiniteElem.h"
#include "IBoundaryElem.h"
#include "LinElem.h"
#include "IBoundaryElem.h"


class Builder
{
public:

	// N - number of finite elements along axis
	// left - left bound of the calculation area
	// rigth - rigth bound of the calculation area
	static FemGrid BuildLinear1DGrid(double left, double rigth, size_t N);
};

#endif