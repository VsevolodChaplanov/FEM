#ifndef __FINITELEMENTSGRID__
#define __FINITELEMENTSGRID__

#include <fstream>
#include "IFiniteElem.h"
#include "IBoundaryElem.h"

class FemGrid
{
public:

	std::vector<double> vertices;
	std::vector<IFiniteElement> elements;
	std::vector<IBoundaryElement> boundary_elements;

private:

	size_t dim;
	size_t vertices_number;
	size_t elements_number;

public:

	FemGrid(size_t dim, std::vector<double> &vertices, std::vector<IFiniteElement> &elements, std::vector<IBoundaryElement> &boundary_element);
	size_t get_elements_number();
	size_t get_vertices_number();
	double* get_vertex(size_t i);
	~FemGrid();

	// Save data as file of .vtk format
	void savevtk(const std::vector<double> &, const std::string &);

private:

	// Dummy function to obtain number of required elements to function "savevtk"
	size_t _cell_size_();
};



#endif