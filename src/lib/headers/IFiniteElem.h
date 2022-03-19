#ifndef __INTERFACE_FINITE_ELEMS__
#define __INTERFACE_FINITE_ELEMS__

#include <vector>
#include <functional>
#include <string>

// --------- Finite elements interface ---------
class IFiniteElement
{
public:

	// Enumerate vtk format types of finite elements
	enum ElementVTK_Type {
		VTK_VERTEX = 1,
		VTK_POLY_VERTEX,
		VTK_LINE,
		VTK_POLY_LINE,
		VTK_TRIANGLE,
		VTK_TRIANGLE_STRIP,
		VTK_POLYGON,
		VTK_PIXEL,
		VTK_QUAD,
		VTK_TETRA,
		VTK_VOXEL,
		VTK_HEXAHEDRON
	};
	
	static IFiniteElement* Factory(const std::vector<double> &vertices, const std::vector<std::size_t> &GIndices, ElementVTK_Type ElementType); 
	// Returns mass matrix element with local indices [i,j]
	virtual double get_mass(size_t i, size_t j) = 0; 
	// Returns stiffness matrix elements with local indices [i,j]
	virtual double get_stiffness(size_t i, size_t j) = 0;
	// Returns [i,i] element of lumped mass matrix
	virtual double get_lumped(size_t i) = 0;
	// Returns point in parametric space [x,y,z]
	virtual double* phys_to_param(double* point) = 0;
	// Returns point in
	virtual double* param_to_phys(double* p_point) = 0;
	// Center coordinates
	virtual double* get_center_coordinates() = 0;
	virtual ~IFiniteElement() = 0;
	

public:
	// Dimension of finite element
	size_t dim;
	// Number of basis functions
	size_t n_basis;
	// Element type according to vtk format
	size_t element_type;
	// Coordinates of center

	// Massive of basis functions
	std::vector<std::function<double(double*)>> basis_functions;
	// Global indines of local vertexes
	std::vector<size_t> global_indices;
	// Local mass matrix
	std::vector<double> mass_matrix;
	// Local stiffness matrux
	std::vector<double> stiffness_matrix;
	// Lumped mass matrix of element
	std::vector<double> lumped_mass_matrix;
	// Volume? of the finite element
	double det_j;
};


// Строит элемент на двух вершинах
// Хранит в себе
// 			* Базисные фукнции этих вершин
// 			* Глобальные индексы верших этого элемента
//			* Собирает в себе
//				** локлаьную матрицу жесткости
//				** локальную матрицу масс
//				** матрицу вз. масс
//
// Отрезок*
//
//			a ------------ b
//
// Хранится в параметроической плоскости как отрезок:
//
// 			0 ------------ 1
//
//
class LinElem : public IFiniteElement
{
public:

	// vertices a massive of two points [a,b] between which the finite element is constructed
	// GIndices - global indices of points "a" and "b", they will have local indices 0 and 1 
	LinElem(const std::vector<double> &vertices, const std::vector<std::size_t> &GIndices);
	double get_mass(size_t i, size_t j) override;														// Возвращает элемент локальной матрицы масс стояший на позиции [i,j]
	double get_stiffness(size_t i, size_t j) override;														// Возвращает элемент локальной матрицы жексткости стояший на позиции [i,j]
	double get_lumped(size_t i) override;																			// Возвращает элемент lumped mass matrix стояший на позиции [i,j]
	double* phys_to_param(double* point) override;
	double* param_to_phys(double* p_point) override;
	double* get_center_coordinates() override; 
	

private: // Methods

	double phi1(double* _param_point);	// Mathes to the left vertex
	double phi2(double* _param_point);	// Mathes to the right vertex

private: // Properties

	double length;
	double start_point;

public:
};
#endif
