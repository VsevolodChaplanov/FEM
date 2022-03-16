#ifndef __PARTIAL_DIFF_EQUATION_CPP__
#define __PARTIAL_DIFF_EQUATION_CPP__

#include "FemPDE.h"

FemPDE::FemPDE(FemGrid* finite_element_mesh, double (*f_analytical)(double*), double (*k_analytical)(double*))
{
	this->f_func = f_analytical;
	this->k_func = k_analytical;
	this->fin_elem_mesh = finite_element_mesh;
	this->Nvert = finite_element_mesh->get_vertices_number();
	this->Nelem = finite_element_mesh->get_elements_number(); 
	M_g.resize(Nvert);
	S_g.resize(Nvert);
	lhs_g.resize(Nvert);
	rhs_g.resize(Nvert);
}

void FemPDE::assemble()
{
	std::vector<double> k_vec = make_k_vec_center();
	for (auto ifinitelement : fin_elem_mesh->elements)
	{
		for (size_t i = 0; i < ifinitelement.Nbasis; i++)
		{
			for (size_t j = 0; j < ifinitelement.Nbasis; j++)
			{
				size_t GlobI = ifinitelement.GIndices[i];
				size_t GlobJ = ifinitelement.GIndices[j];
				M_g.SetValue(GlobI, GlobJ, M_g.GetValue(GlobI, GlobJ) + ifinitelement.GetMass(i,j));
				S_g.SetValue(GlobI, GlobJ, S_g.GetValue(GlobI, GlobJ) + ifinitelement.GetStiffness(i,j));
			}
		}
	}
	rhs_g = M_g * make_f_vec_center();
	SummCM(&M_g, &S_g, &lhs_g); // DummyFunc 
}

std::vector<double> FemPDE::solve(const std::string &Method, const double omega = 0)
{
	// Тут потом сделать чтобы фабрика возвращала std::unique_ptr
	IMatrixSolver* solver = IMatrixSolver::Fabric(Method, omega);
	std::vector<double> solution;
	solver->solve(lhs_g, rhs_g, solution);
	return solution;
}

std::vector<double> FemPDE::make_f_vec_center()
{
	std::vector<double> f_vec(Nvert);
	for (size_t i = 0; i < Nelem; i++)
	{
		f_vec[i] = f_func(fin_elem_mesh->get_vertex(i));
	}
	return f_vec;
}

std::vector<double> FemPDE::make_k_vec_center()
{
	std::vector<double> k_vec(Nvert);
	for (size_t i = 0; i < Nelem; i++)
	{
		k_vec[i] = k_func(fin_elem_mesh->elements[i].get_coord_of_center_of_finelem());
	}
	return k_vec;
}
#endif
