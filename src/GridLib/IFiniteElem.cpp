#ifndef __INTERFACE_FINITE_ELEMS_CPP__
#define __INTERFACE_FINITE_ELEMS_CPP__

#include "IFiniteElem.h"
#include "LinElem.h"
#include "RectElem.h"
#include "CubeElem.h"

IFiniteElement* IFiniteElement::Factory(const std::vector<double> &vertices, const std::vector<std::size_t> &GIndices, const std::size_t ElementType)
{
	if (ElementType == 3)
	{
		return new LinElem(vertices, GIndices);
	} else if(ElementType == 8) 
	{
		return new RectElem(vertices, GIndices);
	} else if(ElementType == 11)
	{
		return new CubeElem(vertices, GIndices);
	}
}

#endif