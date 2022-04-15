#include "IBoundaryElement.h"
#include "LinearPointBoundaryElement.h"
#include "LinearLineBoundaryElement.h"

IBoundaryElement* IBoundaryElement::Factory(const std::vector<double> &vertices, const std::vector<size_t> &GIndices, size_t element_type, size_t bound_type)
{
	if (element_type == 0)
	{
		return new LinearPointBoundaryElement(vertices, GIndices, bound_type);
	} else if (element_type == 1)
	{
		return new LinearLineBoundaryElement(vertices, GIndices, bound_type);
	}
	
	throw std::runtime_error("Can't assemble boundary element of this element type");
}

IBoundaryElement::~IBoundaryElement() { }