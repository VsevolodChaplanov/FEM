#ifndef __MATRIX_SOLVERS__
#define __MATRIX_SOLVERS__

#include <vector>
#include <iostream>
#include "../CompressedMatrix/CompressedM.cpp"
#include "Preconditioners.cpp"
#include "VectorOperations.cpp"

/*-----------------------------Solvers Interface-----------------------------*/
class IMatrixSolver
{
protected:

	const std::size_t MAX_ITERATIONS = 10000000;
	const double eps = 0.0001;
	std::size_t Save_steps = 10;
	std::vector<double> R;

public:
	virtual void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) = 0;
	static IMatrixSolver* Fabric(const std::string&);
	static IMatrixSolver* Fabric(const std::string&, const double omega);
	static IMatrixSolver* Fabric_P(const std::string&, const std::string&);
	static IMatrixSolver* Fabric_P(const std::string&, const std::string&, const double omega);
};
/*-----------------------------Solvers Interface-----------------------------*/


/*-----------------------------Gradient decent solver-----------------------------*/
class GD_Solver : public IMatrixSolver
{
public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Gradient decent solver-----------------------------*/


/*-----------------------------Minimal residuals solver-----------------------------*/
class MR_Solver : public IMatrixSolver
{
public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Minimal residuals solver-----------------------------*/


/*-----------------------------Conjugate gradients solver-----------------------------*/
class CG_Solver : public IMatrixSolver
{
public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Conjugate gradients solver-----------------------------*/


/*-----------------------------Symmetric Successive over relaxation solver-----------------------------*/
class SSOR_Solver : public IMatrixSolver
{
protected:

	double omega = 1.95;

public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
	void SetOmega(double omega);
};
/*-----------------------------Symmetric Successive over relaxation solver-----------------------------*/


/*-----------------------------Successive over relaxation solver-----------------------------*/
class SOR_Solver : public IMatrixSolver
{
protected:

	double omega = 1.95;

public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
	void SetOmega(double omega);
};
/*-----------------------------Successive over relaxation solver-----------------------------*/


/*-----------------------------Seidel-Gauss solver-----------------------------*/
class Seidel_Solver : public IMatrixSolver
{
public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Seidel-Gauss solver-----------------------------*/


/*-----------------------------Jacobi solver-----------------------------*/
class Jacobi_Solver : public IMatrixSolver
{
public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Jacobi solver-----------------------------*/


/*-----------------------------Thomas solver-----------------------------*/
class Thomas_Solver : public IMatrixSolver
{
protected:

	bool CheckTridiagonal(CMatrix& Lhs);

public:
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Thomas solver-----------------------------*/


/*-----------------------------Cholesky factorization solver-----------------------------*/
class LU_Solver : public IMatrixSolver
{
protected:

	void LU_decomposition(CMatrix &Lhs, CMatrix &L);

public:

	void solve(CMatrix &Lhs, const std::vector<double> &Rhs, std::vector<double> &u) override;

};
/*-----------------------------Cholesky factorization solver-----------------------------*/


/*-----------------------------LDU factorization solver-----------------------------*/
class LDU_Solver : public IMatrixSolver
{
protected:

	void LDU_decomposition(CMatrix &Lhs, CMatrix &L, CMatrix &D);

public:

	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------LDU factorization solver-----------------------------*/



// Возможно есть смысл изначально унаследовать отдельный интерфейс для методов предобсуловленных

/*-----------------------------Gradient decent solver (P)-----------------------------*/
class GD_Solver_P : public IMatrixSolver
{
protected:
	IPreconditioner* Preconditioner;
public:

	GD_Solver_P(IPreconditioner* Preconditioner);
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Gradient decent solver (P)-----------------------------*/


/*-----------------------------Conjugate gradients solver (P)-----------------------------*/
class CG_Solver_P : public IMatrixSolver
{
protected:
	IPreconditioner* Preconditioner;
public:

	CG_Solver_P(IPreconditioner* Preconditioner);
	void solve(CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u) override;
};
/*-----------------------------Conjugate gradients solver (P)-----------------------------*/


#endif
