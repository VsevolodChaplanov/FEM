// #ifndef __INTERFACE_FINITE_ELEMS_CPP__
// #define __INTERFACE_FINITE_ELEMS_CPP__

#include "../headers/IFiniteElem.h"
#include "../headers/LinElem.h"

IFiniteElement* IFiniteElement::Factory(const std::vector<double> &verices, const std::vector<size_t> &GIndices, ElementVTK_Type ElementType)
{
	if (ElementType == ElementVTK_Type::VTK_LINE)
	{
		return new LinElem(verices, GIndices);
	}
	return new LinElem(verices, GIndices);
}

IFiniteElement::~IFiniteElement() { }

// #endif