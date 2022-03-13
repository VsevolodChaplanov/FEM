#ifndef __PARTIAL_DIFF_EQUATION__
#define __PARTIAL_DIFF_EQUATION__

#include "../GridLib/Builder1D.h"
#include "../GridLib/Builder1D.cpp"
#include "../CompressedMatrix/CompressedM.h"
#include "../CompressedMatrix/CompressedM.cpp"
#include "../Solvers/Solvers.hpp"
#include "../Solvers/Solvers.cpp"
#include <functional>
#include <string>


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
	FEM_PDE(const std::size_t NN);
	void AssembleSystem(Builder1D* Elements, const std::vector<double> &f_func_vec);
	void Solve(const std::string &Method);
};

#endif

