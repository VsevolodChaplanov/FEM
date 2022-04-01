// #ifndef __LINEAR_FINITE_ELEMS_CPP__
// #define __LINEAR_FINITE_ELEMS_CPP__

#include "IFiniteElem.h"
#include "LinElem.h"
#include "IFiniteElem.h"

LinElem::LinElem(const std::vector<double> &verteces, const std::vector<std::size_t> &GIdences) :
	length(verteces[1] - verteces[0]),
	det_j(verteces[1] - verteces[0]),
	start_point(verteces[0]),
	volume(length),
	center((verteces[1] + verteces[0]) / 2),
	global_indices(GIdences)
{
	this->mass_matrix = {length / 3, length / 6, length / 6, length / 3};
	this->stiffness_matrix = {1 / length, - 1 / length, - 1 / length, 1 / length};
	this->lumped_mass_matrix = {length / 2, length / 2};
}


double LinElem::phi1(const double* _param_point) const
{
	return 1 - _param_point[0];
}

double LinElem::phi2(const double* _param_point) const
{
	return _param_point[0];
}

double LinElem::get_mass(size_t i, size_t j) const
{
	return mass_matrix[i * 2 + j];
}

double LinElem::get_stiffness(size_t i, size_t j) const
{
	return stiffness_matrix[i * 2 + j];
}

double LinElem::get_lumped(size_t i) const
{
	return lumped_mass_matrix[i];
}

void LinElem::phys_to_param(const double* phys_in, double* param_out) const
{
	param_out[0] = (phys_in[0] - start_point) / length;
}

void LinElem::param_to_phys(const double* param_in, double* phys_out) const
{
	phys_out[0] = param_in[0] * length + start_point;
}

const double* LinElem::get_center_coordinates() const
{
	return &center;
}

double LinElem::get_volume() const
{
	return volume;
}

size_t LinElem::get_number_basis_func() const
{
	return 2;
}

const std::vector<size_t>& LinElem::get_global_indices() const
{
	return global_indices;
}

size_t LinElem::get_element_type() const
{
	return 1;
}

LinElem::~LinElem() { }

const std::vector<double> LinElem::get_mass_matrix() const
{
	std::vector<double> local (mass_matrix.begin(), mass_matrix.end());
	return local;
}
const std::vector<double> LinElem::get_stiffness_matrix() const
{
	std::vector<double> local (stiffness_matrix.begin(), stiffness_matrix.end());
	return local;
}
const std::vector<double> LinElem::get_lumped_matrix() const
{
	std::vector<double> local (lumped_mass_matrix.begin(), lumped_mass_matrix.end());
	return local;
}

// #endif