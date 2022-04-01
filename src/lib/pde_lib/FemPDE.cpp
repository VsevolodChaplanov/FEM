// #ifndef __PARTIAL_DIFF_EQUATION_CPP__
// #define __PARTIAL_DIFF_EQUATION_CPP__

#include "FemPDE.h"
#include "CompressedM.h"
#include "FemGrid.h"
#include "Preconditioners.h"
#include "Solvers.h"
#include "VectorOperations.h"

FemPDE::FemPDE(FemGrid* finite_element_mesh, 
	double (*f_analytical)(const double*), 
	double (*k_analytical)(const double*), 
	const MatrixSolverParams* parameters) : Nvert(finite_element_mesh->get_vertices_number()),
		Nelem(finite_element_mesh->get_elements_number())
{
	this->f_func = f_analytical;
	this->k_func = k_analytical;
	this->fin_elem_mesh = finite_element_mesh;
	this->parameters = parameters;
}

void FemPDE::assemble()
{
	std::vector<double> k_vec = make_k_vec_center();

	CMatrix M_g(Nvert);
	CMatrix S_g(Nvert);

	lhs_g.resize(Nvert);
	rhs_g.resize(Nvert);

	for (size_t k = 0; k < Nelem; k++)
	{
		const IFiniteElement* ifiniteelement = fin_elem_mesh->get_element(k);
		for (size_t i = 0; i < ifiniteelement->get_number_basis_func(); i++)
		{
			for (size_t j = 0; j < ifiniteelement->get_number_basis_func(); j++)
			{
				size_t GlobI = ifiniteelement->get_global_indices()[i];
				size_t GlobJ = ifiniteelement->get_global_indices()[j];
				M_g.SetValue(GlobI, GlobJ, M_g.GetValue(GlobI, GlobJ) + ifiniteelement->get_mass(i,j));
				S_g.SetValue(GlobI, GlobJ, S_g.GetValue(GlobI, GlobJ) + k_vec[k] * ifiniteelement->get_stiffness(i,j));
			}
		}
	}
	summ_cm(M_g, S_g, lhs_g);
	rhs_g = M_g * fin_elem_mesh->approximate(f_func);
}   

std::vector<double> FemPDE::solve() const
{
	std::vector<double> solution(this->Nvert);

	IMatrixSolver* solver = IMatrixSolver::Factory(parameters);
	solver->solve(lhs_g, rhs_g, solution);
	delete solver;
	return solution;

}

void FemPDE::apply_boundary_condition_dirichlet(double (*u_analytical)(const double*), const std::vector<size_t> &boundary_element_indices)
{
	for (auto index : boundary_element_indices)
	{
		lhs_g.SetZeroRow(index);
		lhs_g.SetValue(index, index, 1);
		rhs_g[index] = u_analytical(fin_elem_mesh->get_vertex(index));
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
		k_vec[i] = k_func(fin_elem_mesh->get_element(i)->get_center_coordinates());
	}
	return k_vec;
}

void FemPDE::new_assembler()
{
	std::vector<double> k_vec = make_k_vec_center();

	GLobalMatrixAssembler* mass_assembler = new GLobalMatrixAssembler(Nvert);
	GLobalMatrixAssembler* stiffness_assembler = new GLobalMatrixAssembler(Nvert);


	lhs_g.resize(Nvert);
	rhs_g.resize(Nvert);

	for (size_t k = 0; k < Nelem; k++)
	{
		const IFiniteElement* ifiniteelement = fin_elem_mesh->get_element(k);
		mass_assembler->add_local_matrix(ifiniteelement->get_global_indices(), ifiniteelement->get_mass_matrix());
		stiffness_assembler->add_local_matrix(ifiniteelement->get_global_indices(), ifiniteelement->get_stiffness_matrix(), k_vec[k]);
	}
	summ_cm(mass_assembler->get_result(), stiffness_assembler->get_result(), lhs_g);
	rhs_g = mass_assembler->get_result() * fin_elem_mesh->approximate(f_func);
}

// #endif