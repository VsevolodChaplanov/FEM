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

	static IFiniteElement* Factory(const std::vector<double> &vertices, const std::vector<std::size_t> &GIndices, const std::size_t ElementType); 
	// Пока возвращает двумерный массив нет смысла в виртуал -> переделать в одномерный
	virtual std::vector<std::vector<double>>* GetMass(const std::size_t i, const std::size_t j) = 0; 
	virtual std::vector<std::vector<double>>* GetStiffness(const std::size_t i, const std::size_t j) = 0;
	virtual std::vector<double>* GetLumped(const std::size_t i, const std::size_t j) = 0;
	virtual std::vector<double> PhysToParam(const double x, const double y = 0, const double z = 0) = 0;
	virtual std::vector<double> PhysToParam(const std::vector<double> &point) = 0;

protected:
	std::size_t Nbasis; 													// Number of basis functions
	std::vector<std::function<std::vector<double>(std::vector<double>&)>> BasisFunctions;		// Массив базисных фукнций
	std::vector<std::size_t> GIndices; 										// Global indices of local vertices
	std::vector<std::vector<double>> MassMatrix;			// Локальная матрица масс
	std::vector<std::vector<double>> StiffnessMatrix; 	// Локальная матрица Жесткости
	std::vector<double> LumpedMassMatrix; 	// Локальная lumped mass матрица
};

#endif
