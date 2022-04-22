// #ifndef __FINITE_ELEMENTS_MESH_BUILDER_CPP__
// #define __FINITE_ELEMENTS_MESH_BUILDER_CPP__

#include "Builder.h"
#include "IFiniteElement.h"
#include "IBoundaryElement.h"
#include "FemGrid.h"
#include "LinearLineElement.h"
#include "FiniteElemMeshParser.h"


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
	double h = (right - left) / static_cast<double> (N);
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

FemGrid Builder::BuildFromFile(const std::string &filename)
{
	VtkFEMParser* parser = new VtkFEMParser(filename);

	parser->load_mesh();

	size_t p_dim = 2; // Dummy var

	std::vector<double> vertices = parser->get_vertices();
	std::vector<std::vector<size_t>> cells = parser->get_cells();
	std::vector<size_t> cell_types = parser->get_cell_types();

	size_t Nelem = parser->get_elements_number();
	size_t Nvert = parser->get_vertices_number();

	delete parser;

	std::vector<IBoundaryElement*> bound_elems;
	std::vector<IFiniteElement*> elements;

	for (size_t i = 0; i < Nelem; i++)
	{
		if (cell_types[i] == 1)
		{
			std::vector<double> t_vert = {vertices[cells[i][0] * 3  + 0], vertices[cells[i][0] * 3  + 1], vertices[cells[i][0] * 3  + 2]};
			IBoundaryElement* bound_elem = IBoundaryElement::Factory(t_vert, {cells[i][0]}, 0, 1); // !!! Затычка
			bound_elems.push_back(bound_elem);
		}
		if (cell_types[i] == 3 && p_dim == 2)
		{
			std::vector<double> t_vert {
				vertices[cells[i][0] * 3  + 0], vertices[cells[i][0] * 3  + 1], vertices[cells[i][0] * 3  + 2],
				vertices[cells[i][1] * 3  + 0], vertices[cells[i][1] * 3  + 1], vertices[cells[i][1] * 3  + 2]
			};
			IBoundaryElement* bound_elem = IBoundaryElement::Factory(
				t_vert, {cells[i][0], cells[i][1]}, 1, 1
			); // !!! Затычка
			bound_elems.push_back(bound_elem);
		} 
		if (cell_types[i] == 5 && p_dim == 2)
		{
			std::vector<double> t_vert {
				vertices[cells[i][0] * 3  + 0], vertices[cells[i][0] * 3  + 1], vertices[cells[i][0] * 3  + 2],
				vertices[cells[i][1] * 3  + 0], vertices[cells[i][1] * 3  + 1], vertices[cells[i][1] * 3  + 2],
				vertices[cells[i][2] * 3  + 0], vertices[cells[i][2] * 3  + 1], vertices[cells[i][2] * 3  + 2]
			};
			IFiniteElement* felem = IFiniteElement::Factory(
				t_vert, {cells[i][0], cells[i][1], cells[i][2]}, 2
			);
			elements.push_back(felem);
		}	
	}

	FemGrid grid(3, vertices, elements, bound_elems);
	return grid;
}



// #endif