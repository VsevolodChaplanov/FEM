#ifndef __PARTIAL_DIFF_EQUATION__
#define __PARTIAL_DIFF_EQUATION__

#include "Builder.h"
#include "IBoundaryElem.h"
#include "IFiniteElem.h"
#include <functional>
#include <string>
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

	std::function<double(double*)> f_func;
	std::function<double(double*)> k_func;
	std::function<double(double*)> u_exact;

	size_t Nelem;
	size_t Nvert;

	// Global mass matrix
	CMatrix M_g;
	// Global stiffness matrix
	CMatrix S_g;
	// lhs_g = (M_g + S_g)
	CMatrix lhs_g;
	// rhs_g = M * f_vector
	std::vector<double> rhs_g;
	
public:

	FemPDE(FemGrid* finite_element_mesh, double (*f_analytical)(double*), double (*k_analytical)(double*));
	void assemble();
	std::vector<double> solve(const std::string &Method, const double omega = 0);
	void apply_boundary_condition_dirichlet(double (*u_analytical)(double*), const std::vector<size_t> &boundary_element_indices);

private:

	std::vector<double> make_f_vec_center();
	std::vector<double> make_k_vec_center();
};

#endif