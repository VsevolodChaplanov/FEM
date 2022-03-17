#ifndef __LINELEM_CPP__
#define __LINELEM_CPP__

#include "../headers/IFiniteElem.h"

IFiniteElement* IFiniteElement::Factory(const std::vector<double> &verices, const std::vector<size_t> &GIndices, ElementVTK_Type ElementType)
{
	if (ElementType == ElementVTK_Type::VTK_LINE)
	{
		return new LinElem(verices, GIndices);
	}
	return new LinElem(verices, GIndices);
}

LinElem::LinElem(const std::vector<double> &verteces, 
	const std::vector<std::size_t> &GIdences) : global_indices(GIdences), 
	start_point(verteces[0]),
	length(abs(verteces[1] - verteces[0]))
{
	this->dim = 1;
	this->n_basis = 2;
	this->det_j = length;
	this->element_type = LinElem::VTK_LINE;
}

double LinElem::phi1(double* _param_point)
{
	return 1 - _param_point[0];
}

double LinElem::phi2(double* _param_point)
{
	return _param_point[0];
}

double LinElem::get_mass(size_t i, size_t j)
{
	return mass_matrix[i * n_basis + j];
}

double LinElem::get_stiffness(size_t i, size_t j)
{
	return stiffness_matrix[i * n_basis + j];
}

double LinElem::get_lumped(size_t i)
{
	return lumped_mass_matrix[i];
}

double* LinElem::phys_to_param(double* point)
{
	size_t N = sizeof(point) / sizeof(double);
	double* p_point = new double[N];
	p_point[0] = (point[0] - start_point) / length;
	return p_point;
}

double* LinElem::param_to_phys(double* p_point)
{
	double* point = new double[1];
	point[0] = start_point + p_point[0] * length;
	return point;
}

double* LinElem::get_center_coordinates()
{
	double* center = new double[1];
	double* param_center = new double{0.5};
	center = param_to_phys(param_center);
	return center;
}



#endif