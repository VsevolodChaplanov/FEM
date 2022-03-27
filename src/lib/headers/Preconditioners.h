#ifndef __PRECONDITION_METHODS__
#define __PRECONDITION_METHODS__

#include "CompressedM.h"
#include "SolverParams.h"
#include <string>
#include <vector>

class IPreconditioner
{
protected:

	const MatrixSolverParams* params;

public:

	IPreconditioner(const MatrixSolverParams* parameters);
	virtual std::vector<double> Precondition(const CMatrix& Lhs, const std::vector<double>& Rhs) = 0;
	static IPreconditioner* Factory(const MatrixSolverParams* params);
};

class Jacobi_P : public IPreconditioner
{
public:

	Jacobi_P(const MatrixSolverParams* parameters);
	std::vector<double> Precondition(const CMatrix& Lhs, const std::vector<double>& Rhs) override;
};

class SSOR_P : public IPreconditioner
{
public:

	SSOR_P(const MatrixSolverParams* parameters);
	std::vector<double> Precondition(const CMatrix& Lhs, const std::vector<double>& Rhs) override;
};

class ISO0_1_P : public IPreconditioner
{
public:

	ISO0_1_P(const MatrixSolverParams* parameters);
	std::vector<double> Precondition(const CMatrix& Lhs, const std::vector<double>& Rhs) override;
};

class ISO0_2_P : public IPreconditioner
{
public:

	ISO0_2_P(const MatrixSolverParams* parameters);
	std::vector<double> Precondition(const CMatrix& Lhs, const std::vector<double>& Rhs) override;
};

#endif