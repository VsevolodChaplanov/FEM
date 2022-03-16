#ifndef __INTERFACE_FINITE_ELEMS__
#define __INTERFACE_FINITE_ELEMS__

#include <cstddef>
#include <vector>
#include <functional>
#include <string>
#include <memory>

class IFiniteElement
{
public:
	enum ElementVTK_Type {
		VTK_VERTEX,
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
	virtual double GetMass(size_t i, size_t j) = 0; 
	// Returns stiffness matrix elements with local indices [i,j]
	virtual double GetStiffness(size_t i, size_t j) = 0;
	// Returns [i,i] element of lumped mass matrix
	virtual double GetLumped(size_t i) = 0;
	// Returns point [x,y,z]
	virtual double* PhysToParam(double* point) = 0;
	// Dummy
	size_t get_number_basis_finctions() { return Nbasis; }
	virtual double* get_coord_of_center_of_finelem() = 0;
	

public:
	size_t dim;
	size_t Nbasis; 													// Number of basis functions
	size_t ElementType;
	std::vector<std::function<std::vector<double>(std::vector<double>&)>> BasisFunctions;		// Массив базисных фукнций
	std::vector<std::size_t> GIndices; 										// Global indices of local vertices
	std::vector<std::vector<double>> MassMatrix;			// Локальная матрица масс
	std::vector<std::vector<double>> StiffnessMatrix; 	// Локальная матрица Жесткости
	std::vector<double> LumpedMassMatrix; 	// Локальная lumped mass матрица
};

#endif
