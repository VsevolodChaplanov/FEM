#ifndef __SOLVERS_PARAMETERS__
#define __SOLVERS_PARAMETERS__

#include <cmath>

struct MatrixSolverParams
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

	MatrixSolverParams::Methods solve_method = MatrixSolverParams::Methods::Thomas;
	MatrixSolverParams::Preconditioners precondition_method = MatrixSolverParams::Preconditioners::None; 
	size_t MAX_ITERATIONS = 10000;
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
	MatrixSolverParams(MatrixSolverParams::Methods method = MatrixSolverParams::Methods::Thomas,
		MatrixSolverParams::Preconditioners precondition_method = MatrixSolverParams::Preconditioners::None,
		size_t max_iterations = 10000,
		double eps = 1.e-5,
		size_t save_steps = 10,
		double omega_solve = 1.95,
		double omega_precondition = 1.95
	);

	// Set params for SLE solver
	void set_params(MatrixSolverParams::Methods method = MatrixSolverParams::Methods::Thomas,
		MatrixSolverParams::Preconditioners precondition_method = MatrixSolverParams::Preconditioners::None,
		size_t max_iterations = 10000,
		double eps = 1.e-5,
		size_t save_steps = 10,
		double omega_solve = 1.95,
		double omega_precondition = 1.95
	);
};

#endif
