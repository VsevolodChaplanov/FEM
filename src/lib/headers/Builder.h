#ifndef __FINITEELEMENTSMESHBUILDER__
#define __FINITEELEMENTSMESHBUILDER__

#include <vector>
#include <cmath>
#include "FemGrid.h"

class Builder
{
public:

	// N - number of finite elements along axis
	// left - left bound of the calculation area
	// rigth - rigth bound of the calculation area
	// Sets to "left" bound boundary type 1
	// Sets to "right" bound boundary type 2
	static FemGrid BuildLinear1DGrid(double left, double rigth, size_t N);
};

#endif