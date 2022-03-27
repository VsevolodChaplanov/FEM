// #ifndef __INTERFACE_FINITE_ELEMS_CPP__
// #define __INTERFACE_FINITE_ELEMS_CPP__

#include "../headers/IFiniteElem.h"
#include "../headers/LinElem.h"

IFiniteElement* IFiniteElement::Factory(const std::vector<double> &verices, const std::vector<size_t> &GIndices, size_t element_type)
{
	if (element_type == 1)
	{
		return new LinElem(verices, GIndices);
	}
	throw std::runtime_error("Can't assemble element of this element type");
}

IFiniteElement::~IFiniteElement() { }

// #endif