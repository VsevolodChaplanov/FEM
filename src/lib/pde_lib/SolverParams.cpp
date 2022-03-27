// #ifndef __SOLVERS_PARAMETERS_CPP__
// #define __SOLVERS_PARAMETERS_CPP__

#include "SolverParams.h"


// Constuctor with selected params for SOR, SSOR methods
MatrixSolverParams::MatrixSolverParams(MatrixSolverParams::Methods method,
	MatrixSolverParams::Preconditioners precondition_method,
	size_t max_iterations,
	double eps,
	size_t save_steps,
	double omega_solve,
	double omega_precondition)
{
	solve_method = method;
	precondition_method = precondition_method;
	MAX_ITERATIONS = MAX_ITERATIONS;
	eps = eps;
	Save_steps = Save_steps;	
	omega_solve = omega_solve;
	omega_preconditioner = omega_precondition;
}

void MatrixSolverParams::set_params( MatrixSolverParams::Methods method,
	MatrixSolverParams::Preconditioners precondition_method,
	size_t MAX_ITERATIONS,
	double eps,
	size_t Save_steps,
	double omega_solve,
	double omega_precondition)
{
	solve_method = method;
	precondition_method = precondition_method;
	MAX_ITERATIONS = MAX_ITERATIONS;
	eps = eps;
	Save_steps = Save_steps;	
	omega_solve = omega_solve;
	omega_preconditioner = omega_precondition;
}

// #endif
