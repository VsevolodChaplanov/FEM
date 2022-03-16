#ifndef __FINITE_ELEMENTS_MESH_BUILDER_CPP__
#define __FINITE_ELEMENTS_MESH_BUILDER_CPP__

#include "Builder.h"

FemGrid Builder::BuildLinear1DGrid(double left, double right, size_t N)
{
	std::vector<IFiniteElement> elements(N);
	std::vector<double> vertices(N + 1);
	std::vector<IBoundaryElement> boundary_elements(2);
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
		// The IFiniteElements has Factory, which can be used here
		elements[i] = LinElem({vertices[i], vertices[i+1]}, {i, i+1});
	}
	
	// Fill boundary elements vector
	boundary_elements[0] = PointBoundaryElement({vertices[0]}, {0});
	boundary_elements[1] = PointBoundaryElement({vertices[N]}, {N});
	
	// Making finite elements mesh
	FemGrid grid(dim, vertices, elements, boundary_elements);

	return grid;
}

#endif