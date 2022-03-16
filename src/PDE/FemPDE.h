#ifndef __PARTIAL_DIFF_EQUATION__
#define __PARTIAL_DIFF_EQUATION__


#include "../GridLib/MeshBuilder.cpp"
#include "../CompressedMatrix/CompressedM.cpp"
#include "../Solvers/Solvers.cpp"
#include <functional>
#include <string>


// class FEM_PDE
// {
// private:

// 	const std::size_t nn; // Number of elements
// 	CMatrix* M; // Mass Matrix
// 	CMatrix* S; // Stiffness Matrix
// 	CMatrix* Lhs;
// 	std::vector<double> Rhs;
// 	std::vector<double> Solution;

// public:

// 	// Получает объект 
// 	FEM_PDE(const std::size_t NN);
// 	void AssembleSystem(ElemsMesh* Elements, const std::vector<double> &f_func_vec);
// 	void Solve(const std::string &Method);
// };

#endif

