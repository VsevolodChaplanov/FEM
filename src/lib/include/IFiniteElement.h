#ifndef __INTERFACE_FINITE_ELEMS__
#define __INTERFACE_FINITE_ELEMS__

#include <vector>
#include <array>
#include <functional>
#include <string>
#include <cmath>
#include "VectorOperations.h"

// 0 - vertex
// 1 - line
//

// --------- Finite elements interface ---------
class IFiniteElement
{
public:

	// enum Elements {
	// 	Line,
	// 	Triangle
	// };
	
	// enum BoundaryElements {
	// 	Point,
	// 	Line
	// };

	static IFiniteElement* Factory(const std::vector<double> &vertices, const std::vector<std::size_t> &GIndices, size_t element_type); 
	// Returns mass matrix element with local indices [i,j]
	virtual double get_mass(size_t i, size_t j) const = 0; 
	// Returns stiffness matrix elements with local indices [i,j]
	virtual double get_stiffness(size_t i, size_t j) const = 0;
	// Returns [i,i] element of lumped mass matrix
	virtual double get_lumped(size_t i) const = 0;
	// Returns point in parametric space [x,y,z]
	virtual void phys_to_param(const double* phys_in, double* param_out) const = 0;
	//virtual double* phys_to_param(double* point) const = 0;
	// Returns point in
	virtual void param_to_phys(const double* param_in, double* phys_out) const = 0;
	//virtual double* param_to_phys(const double* p_point) const = 0;
	// Center coordinates
	virtual const double* get_center_coordinates() const = 0;
	// Get volume of the element
	virtual double get_volume() const = 0;
	// Get number of basis functions 
	virtual size_t get_number_basis_func() const = 0;
	// Get type of element
	virtual size_t get_element_type() const = 0;
	// Get global indices
	virtual const std::vector<size_t>& get_global_indices() const = 0;
	virtual ~IFiniteElement();

	// Returns mass matrix element
	virtual std::vector<double> get_mass_matrix() const = 0;
	// Returns stiffness matrix element					
	virtual std::vector<double> get_stiffness_matrix() const = 0;
	// Returns lumpred mass matrix element
	virtual std::vector<double> get_lumped_matrix() const = 0;
};



#endif
