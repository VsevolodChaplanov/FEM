#include "LinearTriangleElement.h"

LinearTriangleElement::LinearTriangleElement(const std::vector<double> &vertices, 
		const std::vector<size_t> &GIndices) : 
			start_point(&vertices[0]),
			global_indices(GIndices)
{
	if (vertices.size() == 9)
	{
		this->dim = 3;
	} else if (vertices.size() == 6)
	{
		this->dim = 2;
	}

	Jacobi = { vertices[dim * 1] - vertices[0], vertices[dim * 2] - vertices[0],
			vertices[dim * 1 + 1] - vertices[1], vertices[dim * 2 + 1] - vertices[1]};
	det_j = Jacobi[0] * Jacobi[3] - Jacobi[1] * Jacobi[2];
	volume = det_j / 2;

	mass_matrix = { det_j / 12, det_j / 24, det_j / 24,
			det_j / 24, det_j / 12, det_j / 24,
			det_j / 24, det_j / 24, det_j / 12 };

	stiffness_matrix = {1/2.0*((Jacobi[0]-Jacobi[1])*(Jacobi[0]-Jacobi[1])+(Jacobi[2]-Jacobi[3])*(Jacobi[2]-Jacobi[3]))/det_j, 1/2.0*((Jacobi[0]-Jacobi[1])*Jacobi[1]+(Jacobi[2]-Jacobi[3])*Jacobi[3])/det_j, 1/2.0*(Jacobi[0]*(-Jacobi[0]+Jacobi[1])+Jacobi[2]*(-Jacobi[2]+Jacobi[3]))/det_j,
			1/2.0*((Jacobi[0]-Jacobi[1])*Jacobi[1]+(Jacobi[2]-Jacobi[3])*Jacobi[3])/det_j, 1/2.0*(Jacobi[1]*Jacobi[1]+Jacobi[3]*Jacobi[3])/det_j, 1/2.0*(-Jacobi[0]*Jacobi[1]-Jacobi[2]*Jacobi[3])/det_j,
			1/2.0*(Jacobi[0]*(-Jacobi[0]+Jacobi[1])+Jacobi[2]*(-Jacobi[2]+Jacobi[3]))/det_j, 1/2.0*(-Jacobi[0]*Jacobi[1]-Jacobi[2]*Jacobi[3])/det_j, 1/2.0*(Jacobi[0]*Jacobi[0]+Jacobi[2]*Jacobi[2])/det_j};
	
	lumped_mass_matrix = {det_j/6, det_j/6, det_j/6};

	if (dim == 3)
	{
		center = { (vertices[0] + vertices[dim * 1] + vertices[dim * 2]) / 3, 
			(vertices[1] + vertices[dim * 1 + 1] + vertices[dim * 2 + 1]) / 3, 
			(vertices[2] + vertices[dim * 1 + 2] + vertices[dim * 2 + 2]) / 3 };
	} else if (dim == 2)
	{
		center = { (vertices[0] + vertices[dim * 1] + vertices[dim * 2]) / 3, 
			(vertices[1] + vertices[dim * 1 + 1] + vertices[dim * 2 + 1]) / 3 };
	}
}

size_t LinearTriangleElement::get_element_type() const
{
	return 2;
}

double LinearTriangleElement::get_mass(size_t i, size_t j) const
{
	return mass_matrix[i * dim + j];
}

double LinearTriangleElement::get_stiffness(size_t i, size_t j) const
{
	return stiffness_matrix[i * dim + j];
}

double LinearTriangleElement::get_lumped(size_t i) const
{
	return lumped_mass_matrix[i];
}

double LinearTriangleElement::get_volume() const
{
	return volume;
}

size_t LinearTriangleElement::get_number_basis_func() const
{
	return 3;
}

const std::vector<size_t>& LinearTriangleElement::get_global_indices() const
{
	return global_indices;
}

const double* LinearTriangleElement::get_center_coordinates() const
{
	return &center[0];
}

void LinearTriangleElement::phys_to_param(const double* phys_in, double* param_out) const
{
	param_out[0] = Jacobi[3] / det_j * (phys_in[0] - start_point[0]) - Jacobi[1] / det_j * (phys_in[1] - start_point[1]);
	param_out[1] = - Jacobi[2] / det_j * (phys_in[0] - start_point[0]) + Jacobi[0] / det_j * (phys_in[1] - start_point[1]);
}

void LinearTriangleElement::param_to_phys(const double* param_in, double *phys_out) const
{
	phys_out[0] = start_point[0] + Jacobi[0] * param_in[0] + Jacobi[1] * param_in[1];
	phys_out[1] = start_point[1] + Jacobi[2] * param_in[0] + Jacobi[3] * param_in[1];
}

std::vector<double> LinearTriangleElement::get_mass_matrix() const
{
	std::vector<double> local(mass_matrix.begin(), mass_matrix.end());
	return local;
}

std::vector<double> LinearTriangleElement::get_stiffness_matrix() const
{
	std::vector<double> local(stiffness_matrix.begin(), stiffness_matrix.end());
	return local;
}

std::vector<double> LinearTriangleElement::get_lumped_matrix() const
{
	std::vector<double> local(lumped_mass_matrix.begin(), lumped_mass_matrix.end());
	return local;
}

double LinearTriangleElement::phi1(const double* p_point) const
{
	return 1 - p_point[0] - p_point[1];
}

double LinearTriangleElement::phi2(const double* p_point) const
{
	return p_point[0];
}

double LinearTriangleElement::phi3(const double* p_point) const
{
	return p_point[1];
}

LinearTriangleElement::~LinearTriangleElement() { }