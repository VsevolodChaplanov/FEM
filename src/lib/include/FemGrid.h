#ifndef __FINITELEMENTSGRID__
#define __FINITELEMENTSGRID__

#include <fstream>
#include <iostream>
#include <cmath>
#include "IFiniteElement.h"
#include "IBoundaryElement.h"
#include "LinearLineElement.h"
#include "LinearPointBoundaryElement.h"

class FemGrid
{
private:

	// container to store vertices
	// {x0,y0,z0,x1,y1,z1, ... , xn,yn,zn} -> length of container 3 times greater than number of vertices
	std::vector<double> vertices;
	// container which stores finite elements
	std::vector<IFiniteElement*> elements;
	// container which stores boundary elements 
	std::vector<IBoundaryElement*> boundary_elements;
	// Dimension of finite element grid
	const size_t dim;
	// Number of all vertices in fin. element grid
	// TODO Have to operate with numver of basis finction if calculation area
	const size_t vertices_number;
	// Number of all elements in fin. elment grid
	const size_t elements_number;
	// Volume of the claculation area
	double volume = 0;

public:

	// Construct finite elements mesh
	// dim - dimension of mesh
	// vertices {x0,y0,z0,x1,y1,z1, ... , xn,yn,zn}
	// elements - finite elements
	// boundary_elements - boundary elements
	FemGrid(size_t dim, const std::vector<double> &vertices, const std::vector<IFiniteElement*> &elements, const std::vector<IBoundaryElement*> &boundary_elements);
	// Return vector of indices of boundary elements of the specified type
	std::vector<size_t> boundary_element_indices(size_t boundary_element_type) const;
	// Return vector of indices of boundary elements of the specified type
	std::vector<size_t> boundary_element_indices(const std::function<bool(const double*)>) const;
	// Approximate analytical function along mesh
	std::vector<double> approximate(std::function<double(const double*)>) const;
	// Returns number of finite elements in f.el.mesh
	size_t get_elements_number() const;
	// Return number of boundary element plus number of elements
	size_t get_all_elements_number() const;
	// Returns number of vertices in f.el.mesh
	size_t get_vertices_number() const;
	// Returns massive of coordinates of the vertex with global index i
	const double* get_vertex(size_t i) const;
	// Return pointer to finite element with index i
	const IFiniteElement* get_element(size_t i) const;
	// Return pointer to boundary element with index i
	const IBoundaryElement* get_boundary_element(size_t i) const;
	// Calculates second norm of **difference** within approximation
	double norm2(const std::vector<double> &difference) const;
	// Save data as file of .vtk format
	void savevtk(const std::vector<double> &, const std::string &) const;
	~FemGrid();

private:

	// Dummy function to obtain number of required elements to function "savevtk"
	size_t _cell_size_() const;
};


// // Returns coordinates of center of finite element with index i
// const double* get_coord_of_center_of_finelem(size_t i) const;
#endif