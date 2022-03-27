#ifndef __INTERFACE_BOUNDARY_ELEMS__
#define __INTERFACE_BOUNDARY_ELEMS__

#include "IFiniteElem.h"

class IBoundaryElement : public IFiniteElement
{
public:

	/*
	* Returns pointer to boundary finite element
	* @param vertices vertices of the finite element
	* @param GIndices global indices of the vertices of the finite element
	* @param element_type type of boundary element 
	* @param bound_type type of bound element to determine difference between sides of the finite element
	*/
	static IBoundaryElement* Factory(const std::vector<double> &vertices, const std::vector<size_t> &GIndices, size_t element_type, size_t bound_type);
	// Returns type of the boundary element
	virtual size_t get_bound_type() const = 0;
	virtual ~IBoundaryElement();
};

// // Returns coordinated of the vertex of the boundary element
// virtual std::vector<double> get_vertex() const = 0;

#endif