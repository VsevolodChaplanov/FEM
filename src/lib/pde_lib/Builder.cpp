// #ifndef __FINITE_ELEMENTS_MESH_BUILDER_CPP__
// #define __FINITE_ELEMENTS_MESH_BUILDER_CPP__

#include "../headers/Builder.h"
#include "../headers/IFiniteElem.h"
#include "../headers/IBoundaryElem.h"
#include "../headers/FemGrid.h"
#include "../headers/LinElem.h"


FemGrid Builder::BuildLinear1DGrid(double left, double right, size_t N)
{
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
		IFiniteElement* newelem = IFiniteElement::Factory(temp_vert, temp_ind, IFiniteElement::VTK_LINE);
		elements.push_back(newelem);
	}
	
	// TODO Factory for Boundary Elements

	// Fill boundary elements vector
	std::vector<double> temp_left_vert {vertices[0]};
	std::vector<size_t> temp_left_ind {0};
	std::vector<double> temp_right_vert {vertices[N]};
	std::vector<size_t> temp_right_ind {N};

	boundary_elements.push_back(new PointBoundaryElement(temp_left_vert, temp_left_ind, 1));
	boundary_elements.push_back(new PointBoundaryElement(temp_right_vert, temp_right_ind, 2));
	
	// Making finite elements mesh
	FemGrid grid(dim, vertices, elements, boundary_elements);
	return grid;
}

// #endif