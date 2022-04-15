#ifndef __LINEAR_FINITE_ELEMS__
#define __LINEAR_FINITE_ELEMS__

#include "IFiniteElement.h"

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
class LinearLineElement : public IFiniteElement
{
public:

	// vertices a massive of two points [a,b] between which the finite element is constructed
	// GIndices - global indices of points "a" and "b", they will have local indices 0 and 1 
	LinearLineElement(const std::vector<double> &vertices, const std::vector<size_t> &GIndices);
	// Returns mass matrix element [i][j]
	double get_mass(size_t i, size_t j) const override;
	// Returns stiffness matrix element [i][j]						
	double get_stiffness(size_t i, size_t j) const override;
	// Returns lumpred mass matrix element [i][i]
	double get_lumped(size_t i) const override;
	// Converts point in physical space to parametric space
	void phys_to_param(const double* phys_in, double* param_out) const override;
	// Converts point from parametric space to physical space
	void param_to_phys(const double* param_in, double *phys_outt) const override;
	// Returns coordinates of center of the finite element
	const double* get_center_coordinates() const override; 
	// Get volume of the element
	double get_volume() const override;
	// Get number of basis functions 
	size_t get_number_basis_func() const override;
	// Get element type
	size_t get_element_type() const override;
	// Get global indices
	const std::vector<size_t>& get_global_indices() const override;
	~LinearLineElement();

	// Returns mass matrix element
	std::vector<double> get_mass_matrix() const override;
	// Returns stiffness matrix element					
	std::vector<double> get_stiffness_matrix() const override;
	// Returns lumpred mass matrix element
	std::vector<double> get_lumped_matrix() const override;
	
private: // Methods

	// Basis function mathes to the left vertex
	double phi1(const double* _param_point) const;
	// Basis function mathes to the right vertex
	double phi2(const double* _param_point) const;

private: // Properties

	// Dimension
	size_t dim = 1; 
	// Volume of the linear elem
	const double volume;
	// Reference point in physical space
	const double* start_point;
	// Coordinate of the center of the finite element
	std::vector<double> center;
	// Massive of basis functions
	std::array<std::function<double(const double*)>, 2> basis_functions;
	// Global indines of local vertexes
	std::vector<size_t> global_indices; // Тут оставил вектор, чтобы вечно не гонять array -> vector -> array
	// Local mass matrix
	std::array<double, 4> mass_matrix;
	// Local stiffness matrux
	std::array<double, 4> stiffness_matrix;
	// Lumped mass matrix of element
	std::array<double, 2> lumped_mass_matrix; // В русской литературе "Вектор нагрузки"
};

#endif
