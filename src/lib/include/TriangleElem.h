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
	// Returns mass matrix element [i][j]
	double get_mass(size_t i, size_t j) const override;
	// Returns stiffness matrix element [i][j]						
	double get_stiffness(size_t i, size_t j) const override;
	// Returns lumpred mass matrix element [i][i]
	double get_lumped(size_t i) const override;
	// Converts point in physical space to parametric space
	void phys_to_param(const double* phys_in, double* param_out) const override;
	// Converts point from parametric space to physical space
	void param_to_phys(const double* param_in, double *phys_out) const override;
	// Returns coordinates of center of the finite element
	const double* get_center_coordinates() const override; 
	// Get volume of the element
	double get_volume() const override;
	// Get number of basis functions 
	size_t get_number_basis_func() const override;
	// Get element type
	size_t get_element_type() const override;
	// Get global indices
	const std::vector<size_t>& get_global_indices() const override;
	~TriangleElem();

	// Returns mass matrix element
	const std::vector<double> get_mass_matrix() const override;
	// Returns stiffness matrix element					
	const std::vector<double> get_stiffness_matrix() const override;
	// Returns lumpred mass matrix element
	const std::vector<double> get_lumped_matrix() const override;

private:

	// Determinant of Jacobi matrix
	double det_j;
	// Volume of the triangle elem
	double volume;
	// Reference point in physical space
	const double* start_point;
	// Dimension
	size_t dim = 2;
	// Coordinate of the center of the finite element
	std::vector<double> center;
	// Massive of basis functions
	std::array<std::function<double(const double*)>, 3> basis_functions;
	// Jacobi matrix
	std::array<double, 4> Jacobi;
	// Global indines of local vertexes
	std::vector<size_t> global_indices; // Тут оставил вектор, чтобы вечно не гонять array -> vector -> array
	// Local mass matrix
	std::array<double, 9> mass_matrix;
	// Local stiffness matrux
	std::array<double, 9> stiffness_matrix;
	// Lumped mass matrix of element
	std::array<double, 3> lumped_mass_matrix; // В русской литературе "Вектор нагрузки"
};



#endif