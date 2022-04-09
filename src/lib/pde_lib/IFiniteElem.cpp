// #ifndef __INTERFACE_FINITE_ELEMS_CPP__
// #define __INTERFACE_FINITE_ELEMS_CPP__

#include "IFiniteElem.h"
#include "LinElem.h"
#include "TriangleElem.h"

IFiniteElement* IFiniteElement::Factory(const std::vector<double> &verices, const std::vector<size_t> &GIndices, size_t element_type)
{
	if (element_type == 1)
	{
		return new LinElem(verices, GIndices);
	}
	if (element_type == 2)
	{
		return new TriangleElem(verices, GIndices);
	}
	
	throw std::runtime_error("Can't assemble element of this element type");
}

IFiniteElement::~IFiniteElement() { }

// #endif