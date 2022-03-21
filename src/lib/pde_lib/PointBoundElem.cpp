// #ifndef __BOUNDARY_ELEMS_CPP__
// #define __BOUNDARY_ELEMS_CPP__

#include "../headers/PointBoundElem.h"

PointBoundaryElement::PointBoundaryElement(const std::vector<double> &vertex, const std::vector<size_t> &g_index, size_t boundary_type) : dim(0),
	vertex(vertex),
	n_basis(1),
	bound_type(boundary_type)
{
	this->global_indices = g_index;
}

std::vector<size_t> PointBoundaryElement::get_global_indices_of_boundtype(size_t boundary_type) const
{
	if (this->bound_type == boundary_type)
	{
		return this->global_indices;
	}
	return {};	
}

double PointBoundaryElement::get_mass(size_t i, size_t j) const
{
	return mass_matrix[0];
}

double PointBoundaryElement::get_stiffness(size_t i, size_t j) const
{
	return stiffness_matrix[0];
}

double PointBoundaryElement::get_lumped(size_t i) const
{
	return lumped_mass_matrix[0];
}

double* PointBoundaryElement::phys_to_param(double* point) const
{
	// Or 0?
	return point;
}

size_t PointBoundaryElement::get_dim() const
{
	return dim;
}

size_t PointBoundaryElement::get_number_basis_func() const
{
	return n_basis;
}

size_t PointBoundaryElement::get_bound_type() const
{
	return bound_type;
}

std::vector<double> PointBoundaryElement::get_vertex() const
{
	return vertex;
}

double PointBoundaryElement::phi(double* point) const
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


// #endif