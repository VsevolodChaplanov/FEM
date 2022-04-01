#ifndef __POINT_BOUNDARY_ELEMENT__
#define __POINT_BOUNDARY_ELEMENT__

#include "IBoundaryElem.h"

// a0 - type "1" boundary condition
// an - type "2" boundary condition

class PointBoundaryElement : public IBoundaryElement
{
public:

	// Gets only 1 vertex
	// vertex - {x0,y0,z0}
	// g_index - {Global index of the vertex}
	// boundary_type - type of boundary element
	// for 1d case:
	//	1 - left bound
	//	2 - right bound
	PointBoundaryElement(const std::vector<double> &vertex, const std::vector<size_t> &g_index, size_t bound_type);
	// Returns mass matrix element with local indices [i,j]
	double get_mass(size_t i, size_t j) const; 
	// Returns stiffness matrix elements with local indices [i,j]
	double get_stiffness(size_t i, size_t j) const;
	// Returns [i,i] element of lumped mass matrix
	double get_lumped(size_t i) const;
	// Returns point in parametric space [x,y,z]
	void phys_to_param(const double* phys_in, double* param_out) const;
	//virtual double* phys_to_param(double* point) const;
	// Returns point in
	void param_to_phys(const double* param_in, double* phys_out) const;
	//virtual double* param_to_phys(const double* p_point) const;
	// Center coordinates
	const double* get_center_coordinates() const;
	// Get volume of the element
	double get_volume() const;
	// Get number of basis functions 
	size_t get_number_basis_func() const;
	// Get type of element
	size_t get_element_type() const;
	// Get global indices
	const std::vector<size_t>& get_global_indices() const;
	// Returns type of the boundary element
	size_t get_bound_type() const;
	~PointBoundaryElement();

	// Returns mass matrix element
	const std::vector<double> get_mass_matrix() const override;
	// Returns stiffness matrix element					
	const std::vector<double> get_stiffness_matrix() const override;
	// Returns lumpred mass matrix element
	const std::vector<double> get_lumped_matrix() const override;

private: 

	double phi(const double* point) const;

private:

	// Bound type
	const size_t bound_type;
	// Global indices of vertexes consisting of the boundary element
	std::vector<size_t> global_indices;
	// Mass matrix of boundary element
	std::array<double, 1> mass_matrix {1};
	// Stiffness matrix of boundary element;
	std::array<double, 1> stiffness_matrix {0};
	// lumped mass matrix of boundary element
	std::array<double, 1> lumped_mass_matrix {1};
	// Vertex storage to function phys_to_param;
	std::vector<double> vertex;
};

#endif