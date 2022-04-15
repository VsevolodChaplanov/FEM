// #ifndef __LINEAR_FINITE_ELEMS_CPP__
// #define __LINEAR_FINITE_ELEMS_CPP__

#include "IFiniteElement.h"
#include "LinearLineElement.h"

LinearLineElement::LinearLineElement(const std::vector<double> &vertices, const std::vector<size_t> &GIdences) :
	volume(vector_lenght(vertices)),
	start_point(&vertices[0]),
	center(centre_vector(vertices)),
	global_indices(GIdences)
{
	this->mass_matrix = {volume / 3.0, volume / 6.0, volume / 6.0, volume / 3.0};
	this->stiffness_matrix = {1.0 / volume, - 1.0 / volume, - 1.0 / volume, 1.0 / volume};
	this->lumped_mass_matrix = {volume / 2.0, volume / 2.0};
	if (vertices.size() == 6)
	{
		this->dim = 3;
	} else if (vertices.size() == 4)
	{
		this->dim = 2;
	}
}


double LinearLineElement::phi1(const double* _param_point) const
{
	return 1 - _param_point[0];
}

double LinearLineElement::phi2(const double* _param_point) const
{
	return _param_point[0];
}

double LinearLineElement::get_mass(size_t i, size_t j) const
{
	return mass_matrix[i * 2 + j];
}

double LinearLineElement::get_stiffness(size_t i, size_t j) const
{
	return stiffness_matrix[i * 2 + j];
}

double LinearLineElement::get_lumped(size_t i) const
{
	return lumped_mass_matrix[i];
}

void LinearLineElement::phys_to_param(const double* phys_in, double* param_out) const
{
	for (size_t i = 0; i < dim; i++)
	{
		param_out[i] = (phys_in[i] - start_point[i]) / volume; 
	}
	
	param_out[0] = (phys_in[0] - start_point[0]) / volume;
}

void LinearLineElement::param_to_phys(const double* param_in, double* phys_out) const
{
	for (size_t i = 0; i < dim; i++)
	{
		phys_out[i] = param_in[i] * volume + start_point[i];
	}
}

const double* LinearLineElement::get_center_coordinates() const
{
	return &center[0];
}

double LinearLineElement::get_volume() const
{
	return volume;
}

size_t LinearLineElement::get_number_basis_func() const
{
	return 2;
}

const std::vector<size_t>& LinearLineElement::get_global_indices() const
{
	return global_indices;
}

size_t LinearLineElement::get_element_type() const
{
	return 1;
}

LinearLineElement::~LinearLineElement() { }

std::vector<double> LinearLineElement::get_mass_matrix() const
{
	std::vector<double> local (mass_matrix.begin(), mass_matrix.end());
	return local;
}
std::vector<double> LinearLineElement::get_stiffness_matrix() const
{
	std::vector<double> local (stiffness_matrix.begin(), stiffness_matrix.end());
	return local;
}
std::vector<double> LinearLineElement::get_lumped_matrix() const
{
	std::vector<double> local (lumped_mass_matrix.begin(), lumped_mass_matrix.end());
	return local;
}

// #endif