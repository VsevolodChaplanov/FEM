#ifndef __PARTIAL_DIFF_EQUATION_CPP__
#define __PARTIAL_DIFF_EQUATION_CPP__

#include "FemPDE.h"

// FEM_PDE::FEM_PDE(const std::size_t NN) : nn(NN)
// {
// 	M = new CMatrix(nn+1);
// 	S = new CMatrix(nn+1);
// 	Lhs = new CMatrix(nn+1);
// 	Rhs.resize(nn+1);
// 	Solution.resize(nn+1);
// }

// void FEM_PDE::AssembleSystem(ElemsMesh* Elements, const std::vector<double> &f_func_vec)
// {
// 	for (size_t k = 0; k < nn; k++)
// 	{
// 		for (size_t i = 0; i < Elements->elems[k]->NBasisFunction; i++)
// 		{
// 			for (size_t j = 0; j < Elements->elems[k]->NBasisFunction; j++)
// 			{
// 				std::size_t GlobI = (Elements->elems[k])->GlobalIndices[i];
// 				std::size_t GlobJ = (Elements->elems[k])->GlobalIndices[j];
// 				M->SetValue(GlobI, GlobJ, M->GetValue(GlobI, GlobJ) + Elements->elems[k]->GetMass(i,j));
// 				S->SetValue(GlobI, GlobJ, S->GetValue(GlobI, GlobJ) + Elements->elems[k]->GetStiff(i,j)); // +Домножение на функцию k
// 			}
// 		}
// 	}
// 	Rhs = (*M) * f_func_vec;
// 	SummCM(M, S, Lhs);
// }

// void FEM_PDE::Solve(const std::string &Method)
// {
// 	IMatrixSolver* solver = IMatrixSolver::Fabric(Method);
// 	solver->solve((*Lhs), Rhs, Solution);
// }

#endif