// #ifndef __FINITELEMENTSGRID_CPP__
// #define __FINITELEMENTSGRID_CPP__

#include "../headers/FemGrid.h"
#include "../headers/IBoundaryElem.h"
#include "../headers/LinElem.h"

FemGrid::FemGrid(size_t dim, const std::vector<double> &vertices, const std::vector<IFiniteElement*> &elements, const std::vector<IBoundaryElement*> &boundary_elements) : elements(elements),
	vertices(vertices),
	boundary_elements(boundary_elements),
	dim(dim),
	vertices_number(vertices.size() / dim),
	elements_number(elements.size()) 
	{ }

std::vector<size_t> FemGrid::boundary_element_indices(size_t boundary_element_type) const
{
	std::vector<size_t> global_indices;
	for (IBoundaryElement* elem : boundary_elements)
	{
		std::vector<size_t> temp = elem->get_global_indices_of_boundtype(boundary_element_type);
		for (auto index : temp)
		{
			global_indices.push_back(index);
		}
	}
	return global_indices;
}

std::vector<double> FemGrid::approximate(double (*analytical_func)(const double*)) const
{
	std::vector<double> appr_an_func;
	for (size_t i = 0; i < vertices_number; i++)
	{
		appr_an_func.push_back(analytical_func(get_vertex(i)));
	}
	return appr_an_func;
}

double FemGrid::norm2(const std::vector<double> &difference) const
{
	double sum = 0;
	double Volume = 0;
	for (IFiniteElement* element : elements)
	{
		Volume += element->get_volume();
	}
	
	for (IFiniteElement* element : elements)
	{
		for (size_t i = 0; i < element->get_number_basis_func(); i++)
		{
			sum += element->get_lumped(i) * difference[element->get_global_indices()[i]] * difference[element->get_global_indices()[i]];  // / elements[i].det_j
		}
	}
	// В формуле из pdf будто бы ошибка и там |D|_i -> если так то исправить добавлением в код коммента выше
	return sqrt(sum);
}

double* FemGrid::get_coord_of_center_of_finelem(size_t i) const
{
	return elements[i]->get_center_coordinates();
}

double* FemGrid::get_vertex(size_t i) const
{
	double* vertex = new double[dim];
	for (size_t j = 0; j < dim; j++)
	{
		vertex[j] = vertices[i * dim + j];
	}
	return vertex;
}

void FemGrid::savevtk(const std::vector<double> &solution, const std::string &filename) const
{
	std::ofstream File;
	File.open(filename);	

	// --------- Header of vtk file ---------
	// required version of vtk datafile
	File << "# vtk DataFile Version 3.0" << std::endl; 
	// name of the dataset
	File << "Finite Elements Method" << std::endl;
	// format of dataset
	File << "ASCII" << std::endl;
	// --------- Definition of the FEM type ---------
	File << "DATASET UNSTRUCTURED_GRID" << std::endl;

	// --------- Part to define finite elements mesh ---------
	File << "POINTS " << vertices_number << " " << "double" << std::endl;
	// --------- Vertices ---------
	for (size_t i = 0; i < vertices_number; i++)
	{
		if (dim <= 1. + 1.e-8)
		{
			File << vertices[dim * i] << " " << 0.0 << " " << 0.0 << std::endl;
			continue;
		}
		if (dim < 3 && dim > 1)
		{
			File << vertices[dim * i] << " " << vertices[dim * i + 1] << " " << 0.0 << std::endl;
			continue;
		}
		for (size_t j = 0; j < dim; j++)
		{
			// looping through coordinate cus vertices store as {x0,y0,z0,x1,y1,z1,..}
			File << vertices[dim * i + j] << " ";
		}
		File << std::endl;
	}
	// --------- Part to define cells ---------
	File << "CELLS " << elements_number << " " << _cell_size_() << std::endl; 
	// --------- Looping through global indices of vertices of the finite elements ---------
	for (size_t i = 0; i < elements_number; i++)
	{
		File << elements[i]->get_number_basis_func();
		for (auto elem : elements[i]->get_global_indices())
		{
			File << " " << elem;
		}
		File << std::endl;
	}
	// -------- Looping through elements to define their type ---------
	File << "CELL_TYPES " << elements_number << std::endl;
	for (IFiniteElement* elem : elements)
	{
		File << elem->get_element_type() << std::endl;
	}

	// -------- Part to store obtained numerical solution ---------
	File << "POINT_DATA " << solution.size() << std::endl;
	File << "SCALARS scalars double 1" << std::endl;
	File << "LOOKUP_TABLE default" << std::endl;
	for (auto& elem : solution)
	{
		File << elem << std::endl;
	}

	std::cout << "Data has been stored in " << filename << std::endl;
}


size_t FemGrid::_cell_size_() const
{
	size_t full_cellssize = 0;
	for (IFiniteElement* elem : elements)
	{
		full_cellssize += (elem->get_number_basis_func()) + 1;
	}
	return full_cellssize;
}

size_t FemGrid::get_elements_number() const
{
	return elements_number;
}

size_t FemGrid::get_vertices_number() const
{
	return vertices_number;
}

// std::vector<IFiniteElement*>* FemGrid::get_elements() const
// {
// 	return &elements;
// }

// std::vector<IBoundaryElement*>* FemGrid::get_boundary_elements() const
// {
// 	return &boundary_elements;
// }

IFiniteElement* FemGrid::get_element(size_t i) const
{
	return elements[i];
}

IBoundaryElement* FemGrid::get_boundary_element(size_t i) const
{
	return boundary_elements[i];
}

FemGrid::~FemGrid() { }


// #endif
