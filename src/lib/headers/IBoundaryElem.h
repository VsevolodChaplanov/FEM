#ifndef __INTERFACE_BOUNDARY_ELEMS__
#define __INTERFACE_BOUNDARY_ELEMS__

#include <vector>

class IBoundaryElement
{
public:

	enum ElementVTK_Type {
		VTK_VERTEX,
		VTK_POLY_VERTEX,
		VTK_LINE,
		VTK_POLY_LINE,
		VTK_TRIANGLE,
		VTK_TRIANGLE_STRIP,
		VTK_POLYGON,
		VTK_PIXEL,
		VTK_QUAD
	};

	// Returns mass matrix element with local indices [i,j]
	virtual double get_mass(size_t i, size_t j) const = 0; 
	// Returns stiffness matrix elements with local indices [i,j]
	virtual double get_stiffness(size_t i, size_t j) const = 0;
	// Returns [i,i] element of lumped mass matrix
	virtual double get_lumped(size_t i) const = 0;
	// Returns global indices according to boundary type
	virtual std::vector<size_t> get_global_indices_of_boundtype(size_t boundary_type) const = 0;
	// Returns dimension of element
	virtual size_t get_dim() const = 0;
	// Returns number of basis functions
	virtual size_t get_number_basis_func() const = 0;
	// Returns bound type
	virtual size_t get_bound_type() const = 0;
	// Returns coordinated of the vertex of the boundary element
	virtual std::vector<double> get_vertex() const = 0;
	// Returns point in parametric space according to physical point
	virtual double* phys_to_param(double* point) const = 0;
};

#endif