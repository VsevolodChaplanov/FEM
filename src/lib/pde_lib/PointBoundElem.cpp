// #ifndef __BOUNDARY_ELEMS_CPP__
// #define __BOUNDARY_ELEMS_CPP__

#include "PointBoundElem.h"

PointBoundaryElement::PointBoundaryElement(const std::vector<double> &vertex, const std::vector<size_t> &g_index, size_t boundary_type) :
	vertex(vertex),
	bound_type(boundary_type),
	global_indices(g_index)
{ }

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

const std::vector<size_t>& PointBoundaryElement::get_global_indices() const
{
	return global_indices;
}

void PointBoundaryElement::phys_to_param(const double* phys_in, double* param_out) const
{
	throw std::runtime_error("Can not calculate phys_to_param for boundary element");
}

void PointBoundaryElement::param_to_phys(const double* param_in, double* phys_out) const
{
	throw std::runtime_error("Can not calculate param_to_phys for boundary element");
}

const double* PointBoundaryElement::get_center_coordinates() const
{
	return &vertex[0];
}

double PointBoundaryElement::get_volume() const
{
	return 0;
}

size_t PointBoundaryElement::get_number_basis_func() const
{
	return 1;
}

size_t PointBoundaryElement::get_bound_type() const
{
	return bound_type;
}

size_t PointBoundaryElement::get_element_type() const
{
	return 0;
}

double PointBoundaryElement::phi(const double* point) const
{
	if (point[0] != vertex[0])
	{
		return 0;
	}
	return 1;
}

PointBoundaryElement::~PointBoundaryElement() { }


const std::vector<double> PointBoundaryElement::get_mass_matrix() const
{
	std::vector<double> local (mass_matrix.begin(), mass_matrix.end());
	return local;
}
const std::vector<double> PointBoundaryElement::get_stiffness_matrix() const
{
	std::vector<double> local (stiffness_matrix.begin(), stiffness_matrix.end());
	return local;
}
const std::vector<double> PointBoundaryElement::get_lumped_matrix() const
{
	std::vector<double> local (lumped_mass_matrix.begin(), lumped_mass_matrix.end());
	return local;
}

// #endif