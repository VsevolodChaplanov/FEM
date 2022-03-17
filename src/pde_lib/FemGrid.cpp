#ifndef __FINITELEMENTSGRID_CPP__
#define __FINITELEMENTSGRID_CPP__

#include "../headers/FemGrid.h"

FemGrid::FemGrid(size_t dim, std::vector<double> &vertices, std::vector<IFiniteElement> &elements, std::vector<IBoundaryElement> &boundary_elements)
{
	this->dim = dim;
	this->boundary_elements = boundary_elements;
	this->elements = elements;
	this->vertices = vertices;
	this->vertices_number = vertices.size() / dim;
	this->elements_number = elements.size();
}

std::vector<size_t> FemGrid::boundary_element_indices(size_t boundary_element_type)
{
	std::vector<size_t> global_indices;
	for (IBoundaryElement& elem : boundary_elements)
	{
		std::vector<size_t> temp = elem.get_g_indices_for_belem_type(boundary_element_type);
		for (auto index : temp)
		{
			global_indices.push_back(index);
		}
	}
	return global_indices;
}

std::vector<double> FemGrid::approximate(double (*analytical_func)(double*))
{
	std::vector<double> appr_an_func;
	for (size_t i = 0; i < vertices_number; i++)
	{
		appr_an_func.push_back(analytical_func(get_vertex(i)));
	}
	return appr_an_func;
}

double FemGrid::norm2(const std::vector<double> &difference)
{
	double sum = 0;
	for (IFiniteElement& element : elements)
	{
		for (size_t i = 0; i < element.n_basis; i++)
		{
			sum += element.get_lumped(i) * difference[element.global_indices[i]];  // / elements[i].det_j
		}
	}
	// В формуле из pdf будто бы ошибка и там |D|_i -> если так то исправить добавлением в код коммента выше
	return 0;
}

double* FemGrid::get_coord_of_center_of_finelem(size_t i)
{
	return elements[i].get_center_coordinates();
}

double* FemGrid::get_vertex(size_t i)
{
	double* vertex = new double[dim];
	for (size_t j = 0; j < dim; j++)
	{
		vertex[j] = vertices[i * dim + j];
	}
	return vertex;
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
	File << "CELLS " << elements_number << " " << _cell_size_() << std::endl; 
	// --------- Looping through global indices of vertices of the finite elements ---------
	for (size_t i = 0; i < elements_number; i++)
	{
		File << elements[i].n_basis;
		for (auto elem : elements[i].global_indices)
		{
			File << " " << elem;
		}
		File << std::endl;
	}
	// -------- Looping through elements to define their type ---------
	File << "CELL_TYPES " << elements_number << std::endl;
	for (IFiniteElement& elem : elements)
	{
		File << elem.element_type << std::endl;
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


size_t FemGrid::_cell_size_()
{
	size_t full_cellssize = 0;
	for (IFiniteElement& elem : elements)
	{
		full_cellssize = elem.n_basis + 1;
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
