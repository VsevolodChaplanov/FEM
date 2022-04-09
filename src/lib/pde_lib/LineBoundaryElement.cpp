#include "LineBoundaryElement.h"

LineBoundaryElement::LineBoundaryElement(const std::vector<double> &vertices, 
			const std::vector<size_t> &GIndices, 
			size_t bound_type) : bound_type(bound_type), global_indices(GIndices)
{
	if (vertices.size() == 6)
	{
		this->dim = 3;
	}
	J = {vertices[2] - vertices[0], vertices[3] - vertices[1]};
	volume = sqrt((vertices[2] - vertices[0]) * (vertices[2] - vertices[0]) + (vertices[3] - vertices[1]) * (vertices[3] - vertices[1]));
	this->mass_matrix = {volume / 3, volume / 6, volume / 6, volume / 3};
	this->stiffness_matrix = {1 / volume, - 1 / volume, - 1 / volume, 1 / volume};
	this->lumped_mass_matrix = {volume / 2, volume / 2};
	center = {(vertices[2] + vertices[0]) / 2, (vertices[3] - vertices[1]) / 2};
	start_point = {vertices[0], vertices[1]};
}

double LineBoundaryElement::get_mass(size_t i, size_t j) const
{
	return mass_matrix[i * 2 + j];
}
double LineBoundaryElement::get_stiffness(size_t i, size_t j) const
{
	return stiffness_matrix[i * j + j];
}
double LineBoundaryElement::get_lumped(size_t i) const
{
	return lumped_mass_matrix[i];
}
void LineBoundaryElement::phys_to_param(const double* phys_in, double* param_out) const
{
	param_out[0] = (phys_in[0] - start_point[0]) / J[0];
}
void LineBoundaryElement::param_to_phys(const double* param_in, double *phys_out) const
{
	phys_out[0] = start_point[0] + J[0] * param_in[0];
	phys_out[1] = start_point[1] + J[1] * param_in[0];
}
const double* LineBoundaryElement::get_center_coordinates() const
{
	return &center[0];
}
double LineBoundaryElement::get_volume() const
{
	return volume;
}
size_t LineBoundaryElement::get_number_basis_func() const
{
	return 2;
}
size_t LineBoundaryElement::get_element_type() const
{
	return 1;
}
const std::vector<size_t>& LineBoundaryElement::get_global_indices() const
{
	return global_indices;
}
size_t LineBoundaryElement::get_bound_type() const
{
	return bound_type;
}

const std::vector<double> LineBoundaryElement::get_mass_matrix() const
{
	std::vector<double> local(mass_matrix.begin(), mass_matrix.end());
	return local;
}
const std::vector<double> LineBoundaryElement::get_stiffness_matrix() const
{
	std::vector<double> local(stiffness_matrix.begin(), stiffness_matrix.end());
	return local;
}
const std::vector<double> LineBoundaryElement::get_lumped_matrix() const
{
	std::vector<double> local(lumped_mass_matrix.begin(), lumped_mass_matrix.end());
	return local;
}

LineBoundaryElement::~LineBoundaryElement()
{
}
