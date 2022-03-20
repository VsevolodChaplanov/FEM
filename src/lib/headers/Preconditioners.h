#ifndef __PRECONDITION_METHODS__
#define __PRECONDITION_METHODS__

#include "CompressedM.h"
#include "SolverParams.h"
#include <string>
#include <vector>

struct SolversParams;

class IPreconditioner
{
protected:

	const SolversParams* params;

public:

	IPreconditioner(const SolversParams* parameters);
	virtual std::vector<double> Precondition(CMatrix& Lhs, const std::vector<double>& Rhs) = 0;
	static IPreconditioner* Factory(const SolversParams* params);
};

class Jacobi_P : public IPreconditioner
{
public:

	Jacobi_P(const SolversParams* parameters);
	std::vector<double> Precondition(CMatrix& Lhs, const std::vector<double>& Rhs) override;
};

class SSOR_P : public IPreconditioner
{
public:

	SSOR_P(const SolversParams* parameters);
	std::vector<double> Precondition(CMatrix& Lhs, const std::vector<double>& Rhs) override;
};

class ISO0_1_P : public IPreconditioner
{
public:

	ISO0_1_P(const SolversParams* parameters);
	std::vector<double> Precondition(CMatrix& Lhs, const std::vector<double>& Rhs) override;
};

class ISO0_2_P : public IPreconditioner
{
public:

	ISO0_2_P(const SolversParams* parameters);
	std::vector<double> Precondition(CMatrix& Lhs, const std::vector<double>& Rhs) override;
};

#endif