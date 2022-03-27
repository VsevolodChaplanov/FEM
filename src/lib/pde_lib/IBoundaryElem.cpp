#include "../headers/IBoundaryElem.h"
#include "../headers/PointBoundElem.h"

IBoundaryElement* IBoundaryElement::Factory(const std::vector<double> &vertices, const std::vector<size_t> &GIndices, size_t element_type, size_t bound_type)
{
	if (element_type == 0)
	{
		return new PointBoundaryElement(vertices, GIndices, bound_type);
	}
	throw std::runtime_error("Can't assemble boundary element of this element type");
}

IBoundaryElement::~IBoundaryElement() { }