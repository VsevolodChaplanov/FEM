#ifndef __PRECONDITION_METHODS__
#define __PRECONDITION_METHODS__

#include "../CompressedMatrix/CompressedM.cpp"
#include <string>
#include <vector>

class IPreconditioner
{
private:

protected:

	std::string P_method;

public:

	virtual std::vector<double> Precondition(CMatrix& Lhs, const std::vector<double>& Rhs) = 0;
	static IPreconditioner* Fabric(const std::string &Precondition_method);
	static IPreconditioner* Fabric(const std::string &Precondition_method, const double omega);
};

class Jacobi_P : public IPreconditioner
{
public:

	std::vector<double> Precondition(CMatrix& Lhs, const std::vector<double>& Rhs) override;
};

class SSOR_P : public IPreconditioner
{
protected:

	double omega = 1.5;

public:

	std::vector<double> Precondition(CMatrix& Lhs, const std::vector<double>& Rhs) override;
	void SetOmega(double omega);
};

class ISO0_1_P : public IPreconditioner
{
public:

	std::vector<double> Precondition(CMatrix& Lhs, const std::vector<double>& Rhs) override;
};

class ISO0_2_P : public IPreconditioner
{
public:

	std::vector<double> Precondition(CMatrix& Lhs, const std::vector<double>& Rhs) override;
};

#endif