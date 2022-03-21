#ifndef __POINT_BOUNDARY_ELEMENT__
#define __POINT_BOUNDARY_ELEMENT__

#include "IBoundaryElem.h"

// a0 - type "1" boundary condition
// an - type "2" boundary condition
class PointBoundaryElement : public IBoundaryElement
{
public:

	enum bound_vert_types { left = 1,
	right = 2
	};
	// Gets only 1 vertex
	// vertex - {x0,y0,z0}
	// g_index - {Global index of the vertex}
	// boundary_type - type of boundary element
	// for 1d case:
	//	1 - left bound
	//	2 - right bound
	PointBoundaryElement(const std::vector<double> &vertex, const std::vector<size_t> &g_index, size_t boundary_type);
	// Returns mass matrix element with local indices [i,j]
	double get_mass(size_t i, size_t j) const override; 
	// Returns stiffness matrix element with local indices [i,j]
	double get_stiffness(size_t i, size_t j) const override;
	// Returns [i,i] element of lumped mass matrix
	double get_lumped(size_t i) const override;
	// Returns global indices according to boundary type
	std::vector<size_t> get_global_indices_of_boundtype(size_t boundary_type) const override;
	// Returns dimension of element
	size_t get_dim() const override;
	// Returns number of basis functions
	size_t get_number_basis_func() const override;
	// Returns bound type
	size_t get_bound_type() const override;
	// Returns coordinated of the vertex of the boundary element
	std::vector<double> get_vertex() const override;
	// Returns point in parametric space according to physical point
	double* phys_to_param(double* point) const override;

private: 

	double phi(double* point) const;

private:

	// Dimension of the boundary element
	const size_t dim;
	// Number of basis functions
	const size_t n_basis;
	// Bound type
	const size_t bound_type;
	// Global indices of vertexes consisting of the boundary element
	std::vector<size_t> global_indices;
	// Mass matrix of boundary element
	std::vector<double> mass_matrix {1};
	// Stiffness matrix of boundary element;
	std::vector<double> stiffness_matrix {0};
	// lumped mass matrix of boundary element
	std::vector<double> lumped_mass_matrix {1};
	// Vertex storage to function phys_to_param;
	std::vector<double> vertex;
};

#endif