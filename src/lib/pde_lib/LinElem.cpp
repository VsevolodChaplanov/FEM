// #ifndef __LINEAR_FINITE_ELEMS_CPP__
// #define __LINEAR_FINITE_ELEMS_CPP__

#include "../headers/IFiniteElem.h"
#include "../headers/LinElem.h"
#include "../headers/IFiniteElem.h"

LinElem::LinElem(const std::vector<double> &verteces, const std::vector<std::size_t> &GIdences) : dim(1),
	n_basis(2),
	element_type(LinElem::VTK_LINE),
	length(verteces[1] - verteces[0]),
	start_point(verteces[0]),
	det_j(verteces[1] - verteces[0]),
	Volume(length),
	global_indices(GIdences)
{
	this->mass_matrix = {length / 3, length / 6, length / 6, length / 3};
	this->stiffness_matrix = {1 / length, - 1 / length, - 1 / length, 1 / length};
	this->lumped_mass_matrix = {length / 2, length / 2};
}


double LinElem::phi1(double* _param_point) const
{
	return 1 - _param_point[0];
}

double LinElem::phi2(double* _param_point) const
{
	return _param_point[0];
}

double LinElem::get_mass(size_t i, size_t j) const
{
	return mass_matrix[i * n_basis + j];
}

double LinElem::get_stiffness(size_t i, size_t j) const
{
	return stiffness_matrix[i * n_basis + j];
}

double LinElem::get_lumped(size_t i) const
{
	return lumped_mass_matrix[i];
}

double* LinElem::phys_to_param(double* point) const
{
	size_t N = sizeof(point) / sizeof(double);
	double* p_point = new double[N];
	p_point[0] = (point[0] - start_point) / length;
	return p_point;
}

double* LinElem::param_to_phys(double* p_point) const
{
	double* point = new double[1];
	point[0] = start_point + p_point[0] * length;
	return point;
}

double* LinElem::get_center_coordinates() const
{
	double* center = new double[1];
	double* param_center = new double{0.5};
	center = param_to_phys(param_center);
	return center;
}

double LinElem::get_volume() const
{
	return Volume;
}

size_t LinElem::get_number_basis_func() const
{
	return n_basis;
}

std::vector<size_t> LinElem::get_global_indices() const
{
	return global_indices;
}

size_t LinElem::get_element_type() const
{
	return element_type;
}

// #endif