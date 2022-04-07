// #ifndef __FINITE_ELEMENTS_MESH_BUILDER_CPP__
// #define __FINITE_ELEMENTS_MESH_BUILDER_CPP__

#include "Builder.h"
#include "IFiniteElem.h"
#include "IBoundaryElem.h"
#include "FemGrid.h"
#include "LinElem.h"


FemGrid Builder::BuildLinear1DGrid(double left, double right, size_t N)
{
	if (right < left)
	{
		throw std::invalid_argument( "Recieved wrong arguments for bounds" );
	}
	
	std::vector<IFiniteElement*> elements;
	std::vector<double> vertices(N + 1);
	std::vector<IBoundaryElement*> boundary_elements;
	size_t dim = 1;

	// Fill vertices vector by uniform grid vertices
	double h = (right - left) / (double) N;
	vertices[0] = left;
	for (size_t i = 1; i < N + 1; i++)
	{
		vertices[i] += i * h;
	}

	// Fill elements vector
	for (size_t i = 0; i < N; i++)
	{
		std::vector<double> temp_vert {vertices[i], vertices[i+1]};
		std::vector<size_t> temp_ind  {i, i+1};
		IFiniteElement* newelem = IFiniteElement::Factory(temp_vert, temp_ind, 1);
		elements.push_back(newelem);
	}

	boundary_elements.push_back(IBoundaryElement::Factory({vertices[0]}, {0}, 0, 1));
	boundary_elements.push_back(IBoundaryElement::Factory({vertices[N]}, {N}, 0, 2));

	// Making finite elements mesh
	FemGrid grid(dim, vertices, elements, boundary_elements);
	return grid;
}



// #endif