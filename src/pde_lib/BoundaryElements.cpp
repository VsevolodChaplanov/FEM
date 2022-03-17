#ifndef __BOUNDARY_ELEMS_CPP__
#define __BOUNDARY_ELEMS_CPP__

#include "../headers/IBoundaryElem.h"

PointBoundaryElement::PointBoundaryElement(const std::vector<double> &vertex, const std::vector<size_t> &g_index, size_t boundary_type) : vertex(vertex)
{
	this->dim = 0;
	this->n_basis = 1;
	this->bound_type = boundary_type;
	this->global_indices = g_index;
}

double PointBoundaryElement::get_mass(size_t i, size_t j)
{
	return mass_matrix[0];
}

double PointBoundaryElement::get_stiffness(size_t i, size_t j)
{
	return stiffness_matrix[0];
}

double PointBoundaryElement::get_lumped(size_t i)
{
	return lumped_mass_matrix[0];
}

double* PointBoundaryElement::phys_to_param(double* point)
{
	// Or 0?
	return point;
}

double PointBoundaryElement::phi(double* point)
{
	for (size_t i = 0; i < sizeof(point) / sizeof(double); i++)
	{
		if (point[i] != vertex[i])
		{
			return 0;
		}
	}
	return 1;
}


#endif