// #ifndef __FINITELEMENTSGRID_CPP__
// #define __FINITELEMENTSGRID_CPP__

#include "FemGrid.h"
#include "IBoundaryElement.h"
#include "LinearLineElement.h"

FemGrid::FemGrid(size_t dim, const std::vector<double> &vertices, const std::vector<IFiniteElement*> &elements, const std::vector<IBoundaryElement*> &boundary_elements) : 
	vertices(vertices),
	elements(elements),
	boundary_elements(boundary_elements),
	dim(dim),
	vertices_number(vertices.size() / dim),
	elements_number(elements.size()) 
{ 
	for (IFiniteElement* element : elements)
	{
		volume += element->get_volume();
	}
}

std::vector<size_t> FemGrid::boundary_element_indices(size_t boundary_element_type) const
{
	std::vector<size_t> global_indices;
	for (const IBoundaryElement* elem : boundary_elements)
	{
		if (elem->get_bound_type() == boundary_element_type)
		{
			std::vector<size_t> temp = elem->get_global_indices();
			for (auto index : temp)
			{
				global_indices.push_back(index);
			}
		}
	}
	return global_indices;
}

std::vector<size_t> FemGrid::boundary_element_indices(const std::function<bool(const double*)> seletor_func) const
{
	std::vector<size_t> global_indices;

	for (const IBoundaryElement* b_elem : boundary_elements)
	{
		std::vector<size_t> indices = b_elem->get_global_indices();

		for (const size_t index : indices)
		{
			if (seletor_func(&vertices[dim * index]))
			{
				global_indices.push_back(index);
			}
		}
	}
	return global_indices;
}

std::vector<double> FemGrid::approximate(std::function<double(const double*)> analytical_func) const
{
	std::vector<double> appr_an_func(vertices_number);
	for (size_t i = 0; i < vertices_number; i++)
	{
		appr_an_func[i] = analytical_func(get_vertex(i));
	}
	return appr_an_func;
}

double FemGrid::norm2(const std::vector<double> &difference) const
{
	double sum = 0;

	
	for (IFiniteElement* element : elements)
	{
		for (size_t i = 0; i < element->get_number_basis_func(); i++)
		{
			sum += element->get_lumped(i) * difference[element->get_global_indices()[i]] * difference[element->get_global_indices()[i]];  // / elements[i].det_j
		}
	}
	// В формуле из pdf будто бы ошибка и там |D|_i -> если так то исправить добавлением в код коммента выше
	return sqrt(sum / volume);
}

const double* FemGrid::get_vertex(size_t i) const
{
	return &vertices[i * dim];
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
		if (dim == 1)
		{
			File << vertices[dim * i] << " " << 0.0 << " " << 0.0 << std::endl;
			continue;
		}
		if (dim == 2)
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
	File << "CELLS " << elements_number + boundary_elements.size() << " " << _cell_size_() << std::endl; 
	// --------- Looping through global indices of vertices of the finite elements ---------
	for (const IBoundaryElement* b_elem : boundary_elements)
	{
		File << b_elem->get_number_basis_func();
		for(const size_t index : b_elem->get_global_indices())
		{
			File << " " << index;
		}
		File << std::endl;
	}
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
	File << "CELL_TYPES " << elements_number + boundary_elements.size() << std::endl;
	for (const IBoundaryElement* b_elem : boundary_elements)
	{
		if (b_elem->get_element_type() == 0)
		{
			File << 1 << std::endl;
		} else if (b_elem->get_element_type() == 1)
		{
			File << 3 << std::endl;
		}
	}
	
	for (IFiniteElement* elem : elements)
	{
		if (elem->get_element_type() == 0)
		{
			File << 1 << std::endl;
		} else if (elem->get_element_type() == 1)
		{
			File << 3 << std::endl;
		} else if (elem->get_element_type() == 2)
		{
			File << 5 << std::endl;
		}
		
	}

	// -------- Part to store obtained numerical solution ---------
	File << "POINT_DATA " << solution.size() << std::endl;
	File << "SCALARS scalars double 1" << std::endl;
	File << "LOOKUP_TABLE default" << std::endl;
	for (auto& elem : solution)
	{
		File << elem << std::endl;
	}

	File.close();

	std::cout << "Data has been stored in " << filename << std::endl;
}


size_t FemGrid::_cell_size_() const
{
	size_t full_cellssize = 0;
	for (const IFiniteElement* elem : elements)
	{
		full_cellssize += (elem->get_number_basis_func()) + 1;
	}
	for (const IBoundaryElement* b_elem : boundary_elements)
	{
		full_cellssize += (b_elem->get_number_basis_func()) + 1;
	}
	return full_cellssize;
}

size_t FemGrid::get_elements_number() const
{
	return elements_number;
}

size_t FemGrid::get_all_elements_number() const
{
	return elements_number + boundary_elements.size();
}

size_t FemGrid::get_vertices_number() const
{
	return vertices_number;
}

const IFiniteElement* FemGrid::get_element(size_t i) const
{
	return elements[i];
}

const IBoundaryElement* FemGrid::get_boundary_element(size_t i) const
{
	return boundary_elements[i];
}

FemGrid::~FemGrid() { }

