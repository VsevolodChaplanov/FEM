#include "FEM_PDE.h"

FEM_PDE::FEM_PDE(Builder1D* Elements, const std::vector<double> &f_func_vec) : nn(Elements->GetNElem())
{
	M->resize(nn);
	S->resize(nn);
	Lhs->resize(nn);
	Rhs.resize(nn);
	for (size_t k = 0; k < nn; k++)
	{
		std::size_t ilvert = k;
		std::size_t irvert = k + 1;
		for (size_t i = 0; i < Elements->elems[k]->NBasisFunction; i++)
		{
			for (size_t j = 0; j < (*(Elements->elems[k])).NBasisFunction; j++)
			{
				std::size_t GlobI = (Elements->elems[k])->GlobalIndices[i];
				std::size_t GlobJ = (Elements->elems[k])->GlobalIndices[j];
				M->SetValue(GlobI, GlobJ, M->GetValue(GlobI, GlobJ) + Elements->elems[k]->GetMass(i,j));
				S->SetValue(GlobI, GlobJ, S->GetValue(GlobI, GlobJ) + Elements->elems[k]->GetStiff(i,j)); // +Домножение на функцию k
			}
		}
	}
	Rhs = (*M) * f_func_vec;
	Lhs = SummCM(M,S);
}
