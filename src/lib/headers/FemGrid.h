#ifndef __FINITELEMENTSGRID__
#define __FINITELEMENTSGRID__

#include <fstream>
#include <iostream>
#include "IFiniteElem.h"
#include "IBoundaryElem.h"

class FemGrid
{
public:

	// container to store vertices
	// {x0,y0,z0,x1,y1,z1, ... , xn,yn,zn} -> length of container 3 times greater than number of vertices
	std::vector<double> vertices;
	// container which stores finite elements
	std::vector<IFiniteElement*> elements;
	// container which stores boundary elements 
	std::vector<IBoundaryElement*> boundary_elements;

private:

	// Dimension of finite element grid
	size_t dim;
	// Number of all vertices in fin. element grid
	// It's equal vertices.size() / dim
	size_t vertices_number;
	// Number of all elements in fin. elment grid
	size_t elements_number;

public:

	// Construct finite elements mesh
	// dim - dimension of mesh
	// vertices {x0,y0,z0,x1,y1,z1, ... , xn,yn,zn}
	// elements - finite elements
	// boundary_elements - boundary elements
	FemGrid(size_t dim, std::vector<double> &vertices, std::vector<IFiniteElement*> &elements, std::vector<IBoundaryElement*> &boundary_elements);
	// Return vector of indices of boundary elements of the specified type
	std::vector<size_t> boundary_element_indices(size_t boundary_element_type);
	// Approximate analytical function along mesh
	std::vector<double> approximate(double (*analytical_func)(double *));
	// Returns coordinates of center of finite element with index i
	double* get_coord_of_center_of_finelem(size_t i);
	// Returns number of finite elements in f.el.mesh
	size_t get_elements_number();
	// Returns number of vertices in f.el.mesh
	size_t get_vertices_number();
	// Returns massive of coordinates of the vertex with global index i
	double* get_vertex(size_t i);
	// Calculates second norm of **difference** within approximation
	double norm2(const std::vector<double> &difference);
	// Save data as file of .vtk format
	void savevtk(const std::vector<double> &, const std::string &);
	~FemGrid();

private:

	// Dummy function to obtain number of required elements to function "savevtk"
	size_t _cell_size_();
};



#endif