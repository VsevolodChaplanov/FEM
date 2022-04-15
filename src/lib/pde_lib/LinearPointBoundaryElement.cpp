// #ifndef __BOUNDARY_ELEMS_CPP__
// #define __BOUNDARY_ELEMS_CPP__

#include "LinearPointBoundaryElement.h"

LinearPointBoundaryElement::LinearPointBoundaryElement(const std::vector<double> &vertices, const std::vector<size_t> &g_index) :
	global_indices(g_index),
	vertex(vertices)
{ }

LinearPointBoundaryElement::LinearPointBoundaryElement(const std::vector<double> &vertex, 
			const std::vector<size_t> &g_index,
			size_t bound_type) : LinearPointBoundaryElement(vertex, g_index)
{ 
	this->bound_type = bound_type;
}

double LinearPointBoundaryElement::get_mass(size_t i, size_t j) const
{
	return mass_matrix[0];
}

double LinearPointBoundaryElement::get_stiffness(size_t i, size_t j) const
{
	return stiffness_matrix[0];
}

double LinearPointBoundaryElement::get_lumped(size_t i) const
{
	return lumped_mass_matrix[0];
}

const std::vector<size_t>& LinearPointBoundaryElement::get_global_indices() const
{
	return global_indices;
}

void LinearPointBoundaryElement::phys_to_param(const double* phys_in, double* param_out) const
{
	throw std::runtime_error("Can not calculate phys_to_param for boundary element");
}

void LinearPointBoundaryElement::param_to_phys(const double* param_in, double* phys_out) const
{
	throw std::runtime_error("Can not calculate param_to_phys for boundary element");
}

const double* LinearPointBoundaryElement::get_center_coordinates() const
{
	return &vertex[0];
}

double LinearPointBoundaryElement::get_volume() const
{
	return 0;
}

size_t LinearPointBoundaryElement::get_number_basis_func() const
{
	return 1;
}

size_t LinearPointBoundaryElement::get_bound_type() const
{
	return bound_type;
}

size_t LinearPointBoundaryElement::get_element_type() const
{
	return 0;
}

double LinearPointBoundaryElement::phi(const double* point) const
{
	if (point[0] != vertex[0])
	{
		return 0;
	}
	return 1;
}

LinearPointBoundaryElement::~LinearPointBoundaryElement() { }


std::vector<double> LinearPointBoundaryElement::get_mass_matrix() const
{
	std::vector<double> local (mass_matrix.begin(), mass_matrix.end());
	return local;
}
std::vector<double> LinearPointBoundaryElement::get_stiffness_matrix() const
{
	std::vector<double> local (stiffness_matrix.begin(), stiffness_matrix.end());
	return local;
}
std::vector<double> LinearPointBoundaryElement::get_lumped_matrix() const
{
	std::vector<double> local (lumped_mass_matrix.begin(), lumped_mass_matrix.end());
	return local;
}

// #endif