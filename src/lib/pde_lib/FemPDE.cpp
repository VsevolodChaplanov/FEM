#ifndef __PARTIAL_DIFF_EQUATION_CPP__
#define __PARTIAL_DIFF_EQUATION_CPP__

#include "../headers/FemPDE.h"
#include "BoundaryElements.cpp"
#include "Builder.cpp"
#include "CompressedM.cpp"
#include "FemGrid.cpp"
#include "LinElem.cpp"
#include "Preconditioners.cpp"
#include "Solvers.cpp"
#include "VectorOperations.cpp"

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
	// for (IFiniteElement* ifinitelement : fin_elem_mesh->elements)
	// {
	// 	for (size_t i = 0; i < ifinitelement->n_basis; i++)
	// 	{
	// 		for (size_t j = 0; j < ifinitelement->n_basis; j++)
	// 		{
	// 			size_t GlobI = ifinitelement->global_indices[i];
	// 			size_t GlobJ = ifinitelement->global_indices[j];
	// 			M_g.SetValue(GlobI, GlobJ, M_g.GetValue(GlobI, GlobJ) + ifinitelement->get_mass(i,j));
	// 			S_g.SetValue(GlobI, GlobJ, S_g.GetValue(GlobI, GlobJ) + ifinitelement->get_stiffness(i,j));
	// 		}
	// 	}
	// }
	for (size_t k = 0; k < Nelem; k++)
	{
		IFiniteElement* ifiniteelement = fin_elem_mesh->elements[k];
		for (size_t i = 0; i < ifiniteelement->n_basis; i++)
		{
			for (size_t j = 0; j < ifiniteelement->n_basis; j++)
			{
				size_t GlobI = ifiniteelement->global_indices[i];
				size_t GlobJ = ifiniteelement->global_indices[j];
				M_g.SetValue(GlobI, GlobJ, M_g.GetValue(GlobI, GlobJ) + ifiniteelement->get_mass(i,j));
				S_g.SetValue(GlobI, GlobJ, S_g.GetValue(GlobI, GlobJ) + k_vec[k] * ifiniteelement->get_stiffness(i,j));
			}
		}
	}
	SummCM(M_g, S_g, lhs_g);
	rhs_g = M_g * make_f_vec_vert();
}   

std::vector<double> FemPDE::solve(const std::string &Method, const double omega)
{
	std::vector<double> solution(this->Nvert);
	if (omega > 0)
	{
		// Тут потом сделать чтобы фабрика возвращала std::unique_ptr
		IMatrixSolver* solver = IMatrixSolver::Fabric(Method, omega);
		solver->solve(lhs_g, rhs_g, solution);
		return solution;
	} else
	{
		IMatrixSolver* solver = IMatrixSolver::Fabric(Method);
		solver->solve(lhs_g, rhs_g, solution);
		return solution;
	}
}

void FemPDE::apply_boundary_condition_dirichlet(double (*u_analytical)(double*), const std::vector<size_t> &boundary_element_indices)
{
	u_exact = u_analytical;
	for (auto index : boundary_element_indices)
	{
		lhs_g.SetZeroRow(index);
		lhs_g.SetValue(index, index, 1);
		rhs_g[index] = u_exact(fin_elem_mesh->get_vertex(index));
	}
}

std::vector<double> FemPDE::make_f_vec_vert()
{
	std::vector<double> f_vec(Nvert);
	for (size_t i = 0; i < Nvert; i++)
	{
		f_vec[i] = f_func(fin_elem_mesh->get_vertex(i));
	}
	return f_vec;
}

std::vector<double> FemPDE::make_k_vec_center()
{
	std::vector<double> k_vec(Nelem);
	for (size_t i = 0; i < Nelem; i++)
	{
		k_vec[i] = k_func(fin_elem_mesh->elements[i]->get_center_coordinates());
	}
	return k_vec;
}

#endif