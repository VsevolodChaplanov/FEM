#ifndef __INTERFACE_FINITE_ELEMS__
#define __INTERFACE_FINITE_ELEMS__

#include <vector>
#include <functional>
#include <string>

// --------- Finite elements interface ---------
class IFiniteElement
{
public:
	// Enumerate vtk format types of finite elements
	enum ElementVTK_Type {
		VTK_VERTEX = 1,
		VTK_POLY_VERTEX,
		VTK_LINE,
		VTK_POLY_LINE,
		VTK_TRIANGLE,
		VTK_TRIANGLE_STRIP,
		VTK_POLYGON,
		VTK_PIXEL,
		VTK_QUAD,
		VTK_TETRA,
		VTK_VOXEL,
		VTK_HEXAHEDRON
	};
	
	static IFiniteElement* Factory(const std::vector<double> &vertices, const std::vector<std::size_t> &GIndices, ElementVTK_Type ElementType); 
	// Returns mass matrix element with local indices [i,j]
	virtual double get_mass(size_t i, size_t j) const = 0; 
	// Returns stiffness matrix elements with local indices [i,j]
	virtual double get_stiffness(size_t i, size_t j) const = 0;
	// Returns [i,i] element of lumped mass matrix
	virtual double get_lumped(size_t i) const = 0;
	// Returns point in parametric space [x,y,z]
	virtual double* phys_to_param(double* point) const = 0;
	// Returns point in
	virtual double* param_to_phys(double* p_point) const = 0;
	// Center coordinates
	virtual double* get_center_coordinates() const = 0;
	// Get volume of the element
	virtual double get_volume() const = 0;
	// Get number of basis functions 
	virtual size_t get_number_basis_func() const = 0;
	// Get type of element
	virtual size_t get_element_type() const = 0;
	// Get global indices
	virtual std::vector<size_t> get_global_indices() const = 0;
	virtual ~IFiniteElement() = 0;
};



#endif
