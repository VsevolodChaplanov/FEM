#ifndef __SOLVERS_PARAMETERS__
#define __SOLVERS_PARAMETERS__

#include <cmath>

struct SolversParams
{
public: 	// Params

	// Exist msolver methods
	enum Methods { Thomas,
		Jacobi,
		Seidel,
		SOR,
		SSOR,
		GD,
		MR,
		CG,
		LU,
		LDU
	};

	// Exist precondition methods
	// Only for GD and CG
	enum Preconditioners { None,
		Jacobi_P,
		SSOR_P,
		ILU01_P,
		ILU02_P
	};

public:

	SolversParams::Methods solve_method = SolversParams::Methods::Thomas;
	SolversParams::Preconditioners precondition_method = SolversParams::Preconditioners::None; 
	size_t MAX_ITERATIONS = 10000000;
	double eps = 1.e-5;
	size_t Save_steps = 10;
	double omega_solver = 1.95;	
	double omega_preconditioner = 1.95;

public: 	// Methods

	// Constructor with default params
	// 	Method - Thomas;
	// 	Without precondition
	// 	MAX_ITERATIONS = 10000000; 
	// 	eps = 1.e-5;
	// 	Save_steps = 10;
	// 	omega_solve = 1.95
	//	omega_preconditioner = 1.95
	SolversParams(SolversParams::Methods method,
		SolversParams::Preconditioners precondition_method,
		size_t max_iterations,
		double eps,
		size_t save_steps,
		double omega_solve,
		double omega_precondition);

	// Set params for SLE solver
	void set_params( SolversParams::Methods method,
		SolversParams::Preconditioners precondition_method,
		size_t max_iterations,
		double eps,
		size_t save_steps, 
		double omega_solve,
		double omega_precondition);	
};

#endif
