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
	virtual double get_mass(size_t i, size_t j) = 0; 
	// Returns stiffness matrix elements with local indices [i,j]
	virtual double get_stiffness(size_t i, size_t j) = 0;
	// Returns [i,i] element of lumped mass matrix
	virtual double get_lumped(size_t i) = 0;
	// Returns global indices according to boundary type
	virtual std::vector<size_t> get_g_indices_for_belem_type(size_t boundary_type) = 0;
	// Returns point in parametric space according to physical point
	virtual double* phys_to_param(double* point) = 0;

public:

	size_t dim;
	size_t n_basis;
	size_t bound_type;
};

// a0 - type "1" boundary condition
// an - type "2" boundary condition
class PointBoundaryElement : public IBoundaryElement
{
public:

	enum bound_vert_types { left = 1,
	right = 2
	};

	// Gets only 1 vertice
	// vertex - {x0,y0,z0}
	// g_index - {Gloval index of the vertex}
	// boundary_type - type of boundary element
	// for 1d case:
	//	1 - left bound
	//	2 - right bound
	PointBoundaryElement(const std::vector<double> &vertex, const std::vector<size_t> &g_index, size_t boundary_type);
	// Returns mass matrix element with local indices [i,j]
	virtual double get_mass(size_t i, size_t j) override; 
	// Returns stiffness matrix element with local indices [i,j]
	virtual double get_stiffness(size_t i, size_t j) override;
	// Returns [i,i] element of lumped mass matrix
	virtual double get_lumped(size_t i) override;
	// Returns global indices according to boundary type
	virtual std::vector<size_t> get_g_indices_for_belem_type(size_t boundary_type) override;
	// Returns point in parametric space according to physical point
	virtual double* phys_to_param(double* point) override;

private:

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
	// Bound type 1 or 2
	size_t bound_type;

private: 

	double phi(double* point);
};

#endif