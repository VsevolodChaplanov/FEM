#ifndef __PARTIAL_DEFF_EQUATION__
#define __PARTIAL_DEFF_EQUATION__

#include "../GridLib/Builder1D.h"
#include "../GridLib/Builder1D.cpp"
#include "../CompressedMatrix/CompressedM.h"
#include "../CompressedMatrix/CompressedM.cpp"
#include <functional>

class FEM_PDE
{
private:

	const std::size_t nn; // Number of elements
	CMatrix* M; // Mass Matrix
	CMatrix* S; // Stiffness Matrix
	CMatrix* Lhs;
	std::vector<double> Rhs;
	std::vector<double> Solution;

public:

	// Получает объект 
	FEM_PDE(Builder1D* Elements, const std::vector<double> &f_func_vec);
	void Solve(const CMatrix &Lhs, const std::vector<double> &Rhs, std::vector<double> Solution);
};

#endif

