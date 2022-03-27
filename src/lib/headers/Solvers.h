#ifndef __MATRIX_SOLVERS__
#define __MATRIX_SOLVERS__

#include <vector>
#include <iostream>
#include <memory>
#include "CompressedM.h"
#include "Preconditioners.h"
#include "SolverParams.h"

/*-----------------------------Solvers Interface-----------------------------*/
class IMatrixSolver
{
protected:

	const MatrixSolverParams* parameters;
	std::vector<double> R;

public:

	static IMatrixSolver* Factory(const MatrixSolverParams* parameters);
	IMatrixSolver(const MatrixSolverParams* parameters);
	virtual void solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) = 0;
};
/*-----------------------------Solvers Interface-----------------------------*/


/*-----------------------------Gradient decent solver-----------------------------*/
class GD_Solver : public IMatrixSolver
{
public:

	GD_Solver(const MatrixSolverParams* parameters) : IMatrixSolver(parameters) { }
	void solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Gradient decent solver-----------------------------*/


/*-----------------------------Minimal residuals solver-----------------------------*/
class MR_Solver : public IMatrixSolver
{
public:

	MR_Solver(const MatrixSolverParams* parameters) : IMatrixSolver(parameters) { }
	void solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Minimal residuals solver-----------------------------*/


/*-----------------------------Conjugate gradients solver-----------------------------*/
class CG_Solver : public IMatrixSolver
{
public:

	CG_Solver(const MatrixSolverParams* parameters) : IMatrixSolver(parameters) { }
	void solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Conjugate gradients solver-----------------------------*/


/*-----------------------------Symmetric Successive over relaxation solver-----------------------------*/
class SSOR_Solver : public IMatrixSolver
{
public:

	SSOR_Solver(const MatrixSolverParams* parameters) : IMatrixSolver(parameters) { }
	void solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Symmetric Successive over relaxation solver-----------------------------*/


/*-----------------------------Successive over relaxation solver-----------------------------*/
class SOR_Solver : public IMatrixSolver
{
public:

	SOR_Solver(const MatrixSolverParams* parameters) : IMatrixSolver(parameters) { }
	void solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Successive over relaxation solver-----------------------------*/


/*-----------------------------Seidel-Gauss solver-----------------------------*/
class Seidel_Solver : public IMatrixSolver
{
public:

	Seidel_Solver(const MatrixSolverParams* parameters) : IMatrixSolver(parameters) { }
	void solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Seidel-Gauss solver-----------------------------*/


/*-----------------------------Jacobi solver-----------------------------*/
class Jacobi_Solver : public IMatrixSolver
{
public:

	Jacobi_Solver(const MatrixSolverParams* parameters) : IMatrixSolver(parameters) { }
	void solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Jacobi solver-----------------------------*/


/*-----------------------------Thomas solver-----------------------------*/
class Thomas_Solver : public IMatrixSolver
{
protected:

	bool CheckTridiagonal(const CMatrix& Lhs);

public:

	Thomas_Solver(const MatrixSolverParams* parameters) : IMatrixSolver(parameters) { }
	void solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Thomas solver-----------------------------*/


/*-----------------------------Cholesky factorization solver-----------------------------*/
class LU_Solver : public IMatrixSolver
{
protected:

	void LU_decomposition(const CMatrix &Lhs, CMatrix &L);

public:

	LU_Solver(const MatrixSolverParams* parameters) : IMatrixSolver(parameters) { }
	void solve(const CMatrix &Lhs, const std::vector<double> &Rhs, std::vector<double> &u) override;

};
/*-----------------------------Cholesky factorization solver-----------------------------*/


/*-----------------------------LDU factorization solver-----------------------------*/
class LDU_Solver : public IMatrixSolver
{
protected:

	void LDU_decomposition(const CMatrix &Lhs, CMatrix &L, CMatrix &D);

public:

	LDU_Solver(const MatrixSolverParams* parameters) : IMatrixSolver(parameters) { }
	void solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------LDU factorization solver-----------------------------*/



// Возможно есть смысл изначально унаследовать отдельный интерфейс для методов предобсуловленных

/*-----------------------------Gradient decent solver (P)-----------------------------*/
class GD_Solver_P : public IMatrixSolver
{
protected:
	IPreconditioner* Preconditioner;
public:

	GD_Solver_P(const MatrixSolverParams* parameters, IPreconditioner* Preconditioner);
	void solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Gradient decent solver (P)-----------------------------*/


/*-----------------------------Conjugate gradients solver (P)-----------------------------*/
class CG_Solver_P : public IMatrixSolver
{
protected:
	IPreconditioner* Preconditioner;
public:

	CG_Solver_P(const MatrixSolverParams* parameters, IPreconditioner* Preconditioner);
	void solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Conjugate gradients solver (P)-----------------------------*/


// Оставил пока 1 структуру т.к. не знаю пока как хорошо сделать чтобы методы SOR SSOR и предобсулавливатель SSOR испольщовали правильные константы


// struct SolversParams_SOR_SSOR : public MatrixSolverParams
// {
// public:
// 	double omega_solver = 1.95;	
// 	double omega_preconditioner = 1.95;	

// public:

// 	// Constructor with default params
// 	// 	Method - Thomas;
// 	// 	Without precondition
// 	// 	MAX_ITERATIONS = 10000000; 
// 	// 	eps = 1.e-5;
// 	// 	Save_steps = 10;
// 	// 	omega = 1.95
// 	SolversParams_SOR_SSOR();
// 	// Constuctor with selected params
// 	SolversParams_SOR_SSOR(MatrixSolverParams::Methods method,
// 		MatrixSolverParams::Preconditioners precondition_method,
// 		size_t max_iterations,
// 		double eps,
// 		size_t save_steps, 
// 		double omega);

// 	// Set params for SLE solver
// 	void set_params( MatrixSolverParams::Methods method = Thomas,
// 	MatrixSolverParams::Preconditioners precondition_method = None,
// 	size_t max_iterations = 10000000,
// 	double eps = 1.e-5,
// 	size_t save_steps = 10, 
// 	double omega = 1.95 );
// };

#endif
