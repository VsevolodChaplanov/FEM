#ifndef __TRIANGLE_FINITE_ELEMENT__
#define __TRIANGLE_FINITE_ELEMENT__

#include "IFiniteElem.h"

class TriangleElem : public IFiniteElement
{
private:

	double phi1(const double* p_point) const;
	double phi2(const double* p_point) const;
	double phi3(const double* p_point) const;

public:
	
	// gets: 3 vertives as {x0, y0, x1, y1, x2, y2}
	// 	3 global indices of this vertices
	TriangleElem(const std::vector<double> &vertices, const std::vector<size_t> &GIndices);
	double get_mass(size_t i, size_t j) const override;
	double get_stiffness(size_t i, size_t j) const override;
	double get_lumped(size_t i) const override;
	double* param_to_phys(double* p_point) const override;
	double* phys_to_param(double* point) const override;
	double* get_center_coordinates() const override;
	double get_volume() const override;
	size_t get_number_basis_func() const override;
	size_t get_element_type() const override;
	std::vector<size_t> get_global_indices() const override;

private:

	std::vector<double> jacobi_matrix;
	const double volume;
	const double det_j;
	const size_t dim;
	const size_t n_basis;
	const size_t element_type;
	std::vector<size_t> global_indices;
	std::vector<double> mass_matrix;
	std::vector<double> stiffness_matrix;
	std::vector<double> lumped_mass_matrix;
};



#endif