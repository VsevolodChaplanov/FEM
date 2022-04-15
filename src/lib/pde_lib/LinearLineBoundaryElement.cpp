#include "LinearLineBoundaryElement.h"

LinearLineBoundaryElement::LinearLineBoundaryElement(const std::vector<double> &vertices, 
			const std::vector<size_t> &GIndices) : volume(vector_lenght(vertices)),
			start_point(&vertices[0]),
			global_indices(GIndices)
{
	this->mass_matrix = {volume / 3, volume / 6, volume / 6, volume / 3};
	this->stiffness_matrix = {1 / volume, - 1 / volume, - 1 / volume, 1 / volume};
	this->lumped_mass_matrix = {volume / 2, volume / 2};
	if (vertices.size() == 6)
	{
		this->dim = 3;
		this->direction = {vertices[3] - vertices[0],
				vertices[4] - vertices[1],
				vertices[5] - vertices[2]};
	} else if (vertices.size() == 4)
	{
		this->dim = 2;
		this->direction = {vertices[3] - vertices[0],
				vertices[4] - vertices[1]};
	}
	
	this->center = centre_vector(vertices);
}


LinearLineBoundaryElement::LinearLineBoundaryElement(const std::vector<double> &vertices, 
			const std::vector<size_t> &GIndices, 
			size_t bound_type) : LinearLineBoundaryElement(vertices, GIndices)
{
	this->bound_type = bound_type;
}

double LinearLineBoundaryElement::get_mass(size_t i, size_t j) const
{
	return mass_matrix[i * 2 + j];
}
double LinearLineBoundaryElement::get_stiffness(size_t i, size_t j) const
{
	return stiffness_matrix[i * j + j];
}
double LinearLineBoundaryElement::get_lumped(size_t i) const
{
	return lumped_mass_matrix[i];
}
void LinearLineBoundaryElement::phys_to_param(const double* phys_in, double* param_out) const
{
	param_out[0] = (phys_in[0] - start_point[0]) / direction[0];
}
void LinearLineBoundaryElement::param_to_phys(const double* param_in, double *phys_out) const
{
	phys_out[0] = start_point[0] + direction[0] * param_in[0];
	phys_out[1] = start_point[1] + direction[1] * param_in[0];
}
const double* LinearLineBoundaryElement::get_center_coordinates() const
{
	return &center[0];
}
double LinearLineBoundaryElement::get_volume() const
{
	return volume;
}
size_t LinearLineBoundaryElement::get_number_basis_func() const
{
	return 2;
}
size_t LinearLineBoundaryElement::get_element_type() const
{
	return 1;
}
const std::vector<size_t>& LinearLineBoundaryElement::get_global_indices() const
{
	return global_indices;
}
size_t LinearLineBoundaryElement::get_bound_type() const
{
	return bound_type;
}

std::vector<double> LinearLineBoundaryElement::get_mass_matrix() const
{
	std::vector<double> local(mass_matrix.begin(), mass_matrix.end());
	return local;
}
std::vector<double> LinearLineBoundaryElement::get_stiffness_matrix() const
{
	std::vector<double> local(stiffness_matrix.begin(), stiffness_matrix.end());
	return local;
}
std::vector<double> LinearLineBoundaryElement::get_lumped_matrix() const
{
	std::vector<double> local(lumped_mass_matrix.begin(), lumped_mass_matrix.end());
	return local;
}

LinearLineBoundaryElement::~LinearLineBoundaryElement() { }
