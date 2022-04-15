#ifndef __PARTIAL_DIFF_EQUATION__
#define __PARTIAL_DIFF_EQUATION__

#include <functional>
#include <string>
#include "Builder.h"
#include "IBoundaryElement.h"
#include "IFiniteElement.h"
#include "CompressedM.h"
#include "Solvers.h"
#include "GlobalAssemblers.h"

//
// - ∇{k(r)∇{u(r)}} + u(r) = f(r)
//
class FemPDE
{
private:

	// 
	const FemGrid* fin_elem_mesh;

	std::function<double(const double*)> f_func;
	std::function<double(const double*)> k_func;
	std::function<double(const double*)> u_exact;

	const MatrixSolverParams* parameters;

	const size_t Nelem;
	const size_t Nvert;

	// lhs_g = (M_g + S_g)
	CMatrix lhs_g;
	// rhs_g = M * f_vector
	std::vector<double> rhs_g;
	
public:

	FemPDE(const FemGrid* finite_element_mesh,
			double (*f_analytical)(const double*),
			double (*k_analytical)(const double*),
			const MatrixSolverParams* parameters);
	void assemble();
	std::vector<double> solve() const;
	void apply_boundary_condition_dirichlet(double (*u_analytical)(const double*), const std::vector<size_t> &boundary_element_indices);
	void apply_boundary_condition_dirichlet(double u, const std::vector<size_t> &boundary_element_indices);
	~FemPDE();

private:

	std::vector<double> make_k_vec_center();
};

#endif