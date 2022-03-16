#ifndef __FINITELEMENTSGRID_CPP__
#define __FINITELEMENTSGRID_CPP__

#include "FemGrid.h"
#include <iostream>

FemGrid::FemGrid(size_t dim, std::vector<double> &vertices, std::vector<IFiniteElement> &elements, std::vector<IBoundaryElement> &boundary_element)
{
	this->dim = dim;
	this->boundary_elements = boundary_element;
	this->elements = elements;
	this->vertices = vertices;
	this->vertices_number = vertices.size();
	this->elements_number = elements.size();
}

void FemGrid::savevtk(const std::vector<double> &solution, const std::string &filename)
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
		for (size_t j = 0; j < dim; j++)
		{
			// looping through coordinate cus vertices store as {x0,y0,z0,x1,y1,z1,..}
			File << vertices[i + j] << " ";
		}
		File << std::endl;
	}
	// --------- Part to define cells ---------
	File << "CELLS " << elements_number << " " << _cell_size_ << std::endl; 
	// --------- Looping through global indices of vertices of the finite elements ---------
	for (size_t i = 0; i < elements_number; i++)
	{
		File << elements[i].Nbasis;
		for (auto elem : elements[i].GIndices)
		{
			File << " " << elem;
		}
		File << std::endl;
	}
	// -------- Looping through elements to define their type ---------
	File << "CELL_TYPES " << elements_number << std::endl;
	for (auto elem : elements)
	{
		File << elem.ElementType << std::endl;
	}

	// -------- Part to store obtained numerical solution ---------
	File << "POINT_DATA " << solution.size() << std::endl;
	File << "SCALARS scalars double 1" << std::endl;
	File << "LOOKUP_TABLE default" << std::endl;
	for (auto elem : solution)
	{
		File << elem << std::endl;
	}

	std::cout << "Data has been stored in " << filename << std::endl;
}


size_t FemGrid::_cell_size_()
{
	size_t full_cellssize = 0;
	for (auto elem : elements)
	{
		full_cellssize = elem.Nbasis + 1;
	}
	return full_cellssize;
}

size_t FemGrid::get_elements_number()
{
	return elements_number;
}

size_t FemGrid::get_vertices_number()
{
	return vertices_number;
}

#endif
