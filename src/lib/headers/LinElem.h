#ifndef __LINEAR_FINITE_ELEMS__
#define __LINEAR_FINITE_ELEMS__

#include "IFiniteElem.h"

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
	// Returns mass matrix element [i][j]
	double get_mass(size_t i, size_t j) const override;
	// Returns stiffness matrix element [i][j]						
	double get_stiffness(size_t i, size_t j) const override;
	// Returns lumpred mass matrix element [i][i]
	double get_lumped(size_t i) const override;
	// Converts point in physical space to parametric space
	double* phys_to_param(double* point) const override;
	// Converts point from parametric space to physical space
	double* param_to_phys(double* p_point) const override;
	// Returns coordinates of center of the finite element
	double* get_center_coordinates() const override; 
	// Get volume of the element
	double get_volume() const override;
	// Get number of basis functions 
	size_t get_number_basis_func() const override;
	// Get element type
	size_t get_element_type() const override;
	// Get global indices
	std::vector<size_t> get_global_indices() const override;
	
private: // Methods

	// Basis function mathes to the left vertex
	double phi1(double* _param_point) const;
	// Basis function mathes to the right vertex
	double phi2(double* _param_point) const;

private: // Properties
	// Volume = length of the linear element
	const double length;
	// Determinant of Jacobi matrix
	const double det_j;
	// Volume of the linear elem
	const double Volume;
	// Reference point in physical space
	const double start_point;
	// Dimension of finite element
	const size_t dim;
	// Number of basis functions
	const size_t n_basis;
	// Element type according to vtk format
	const size_t element_type;
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
};

#endif
