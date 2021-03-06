#ifndef __PARTIAL_DIFF_EQUATION__
#define __PARTIAL_DIFF_EQUATION__

#include <functional>
#include <string>
#include "Builder.h"
#include "IBoundaryElem.h"
#include "IFiniteElem.h"
#include "CompressedM.h"
#include "Solvers.h"

//
// - ∇{k(r)∇{u(r)}} + u(r) = f(r)
//
class FemPDE
{
private:

	// 
	FemGrid* fin_elem_mesh;

	std::function<double(const double*)> f_func;
	std::function<double(const double*)> k_func;
	std::function<double(const double*)> u_exact;

	const MatrixSolverParams* parameters;

	const size_t Nelem;
	const size_t Nvert;

	// @kalininei Может быть их создать на этапе сборки?
	// Global mass matrix
	CMatrix M_g;
	// Global stiffness matrix
	CMatrix S_g;

	// lhs_g = (M_g + S_g)
	CMatrix lhs_g;
	// rhs_g = M * f_vector
	std::vector<double> rhs_g;
	
public:

	FemPDE(FemGrid* finite_element_mesh, double (*f_analytical)(const double*), double (*k_analytical)(const double*), const MatrixSolverParams* parameters);
	void assemble();
	std::vector<double> solve() const;
	void apply_boundary_condition_dirichlet(double (*u_analytical)(const double*), const std::vector<size_t> &boundary_element_indices);

private:

	std::vector<double> make_f_vec_vert();
	std::vector<double> make_k_vec_center();
};

#endif