#ifndef __SOLVERS_PARAMETERS_CPP__
#define __SOLVERS_PARAMETERS_CPP__

#include "../headers/SolverParams.h"


// Constuctor with selected params for SOR, SSOR methods
SolversParams::SolversParams(SolversParams::Methods method = SolversParams::Methods::Thomas,
	SolversParams::Preconditioners precondition_method = SolversParams::Preconditioners::None,
	size_t max_iterations = 0,
	double eps = 0,
	size_t save_steps = 0,
	double omega_solve = 0,
	double omega_precondition = 0)
{
	solve_method = method;
	precondition_method = precondition_method;
	MAX_ITERATIONS = MAX_ITERATIONS;
	eps = eps;
	Save_steps = Save_steps;	
	omega_solve = omega_solve;
	omega_preconditioner = omega_precondition;
}

void SolversParams::set_params( SolversParams::Methods method = SolversParams::Methods::Thomas,
	SolversParams::Preconditioners precondition_method = SolversParams::Preconditioners::None,
	size_t MAX_ITERATIONS = 10000000,
	double eps = 1.e-5,
	size_t Save_steps = 10,
	double omega_solve = 1.95,
	double omega_precondition = 1.95)
{
	solve_method = method;
	precondition_method = precondition_method;
	MAX_ITERATIONS = MAX_ITERATIONS;
	eps = eps;
	Save_steps = Save_steps;	
	omega_solve = omega_solve;
	omega_preconditioner = omega_precondition;
}

#endif
