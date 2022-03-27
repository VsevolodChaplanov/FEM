// #ifndef __PRECONDITION_METHODS_CPP__
// #define __PRECONDITION_METHODS_CPP__


#include "../headers/Preconditioners.h"
#include "../headers/CompressedM.h"



// Hadrcode
IPreconditioner* IPreconditioner::Factory(const MatrixSolverParams* parameters)
{
	if (parameters->precondition_method == MatrixSolverParams::Preconditioners::Jacobi_P)
	{
		//std::shared_ptr<IPreconditioner> preconditioner(new Jacobi_P(parameters));
		return new Jacobi_P(parameters);
	} else if (parameters->precondition_method == MatrixSolverParams::Preconditioners::SSOR_P)
	{
		//std::shared_ptr<IPreconditioner> preconditioner(new SSOR_P(parameters));
		return new SSOR_P(parameters);
	} else if (parameters->precondition_method == MatrixSolverParams::Preconditioners::ILU01_P)
	{
		//std::shared_ptr<IPreconditioner> preconditioner(new ISO0_1_P(parameters));
		return new ISO0_1_P(parameters);
	} else if (parameters->precondition_method == MatrixSolverParams::Preconditioners::ILU02_P)
	{
		//std::shared_ptr<IPreconditioner> preconditioner(new ISO0_2_P(parameters));
		return new ISO0_2_P(parameters);
	}
	return nullptr;
}

IPreconditioner::IPreconditioner(const MatrixSolverParams* parameters)
{
	params = parameters;
}

Jacobi_P::Jacobi_P(const MatrixSolverParams* parameters) : IPreconditioner(parameters) { }

SSOR_P::SSOR_P(const MatrixSolverParams* parameters) : IPreconditioner(parameters) { }

ISO0_1_P::ISO0_1_P(const MatrixSolverParams* parameters) : IPreconditioner(parameters) { }

ISO0_2_P::ISO0_2_P(const MatrixSolverParams* parameters) : IPreconditioner(parameters) { }

std::vector<double> Jacobi_P::Precondition(const CMatrix& Lhs, const std::vector<double>& Rhs)
{
	std::vector<double> Result(Rhs.size());
	for (size_t i = 0; i < Rhs.size(); i++)
	{
		Result[i] = Rhs[i] / Lhs.GetValue(i,i);
	}
	return Result;
}

std::vector<double> SSOR_P::Precondition(const CMatrix &Lhs, const std::vector<double> &Rhs)
{
	size_t N = Lhs.size();

	std::vector<double> z(N);
	std::vector<double> y(N);

	for (size_t i = 1; i <= N; i++)
	{
		double Sum = 0;
		for (auto elem : Lhs[i-1])
		{
			if (elem.first < i - 1)
			{
				Sum += elem.second * z[elem.first] / (2 - params->omega_preconditioner);
			}			
		}
		z[i-1] = params->omega_preconditioner * (2 - params->omega_preconditioner) / Lhs.GetValue(i-1, i-1) * (Rhs[i-1] - Sum);
	}
	

	for (size_t i = N; i > 0; i--)
	{
		double Sum = 0;
		for (auto elem : Lhs[i - 1])
		{
			if (elem.first > i - 1)
			{
				Sum += elem.second * y[elem.first] / (2 - params->omega_preconditioner);
			}				
		}
		y[i-1] = z[i-1] - params->omega_preconditioner * (2 - params->omega_preconditioner) / Lhs.GetValue(i-1, i-1) * Sum;
	}

	return y;
}

// Диагональ предобуславливателя совпадает с диагональню исходной матрицы
// Решается задача Matrix * y = Rhs 
// Как Bw = r
std::vector<double> ISO0_1_P::Precondition(const CMatrix &Lhs, const std::vector<double> &Rhs)
{

	size_t N = Lhs.size();
	std::vector<double> y(N);
	CMatrix D(N);

	for (size_t i = 1; i <= N; i++)
	{	
		double Value = 0;
		for (auto elem : Lhs[i-1])
		{
			if (elem.first < i - 1)
			{
				Value += elem.second * elem.second / D.GetValue(elem.first, elem.first);
				//Value += Lhs.GetValue(elem.first, i - 1) * elem.second / D.GetValue(elem.first, elem.first);
			}
		}
		D.SetValue(i - 1,i - 1 , Lhs.GetValue(i-1,i-1) - Value);
	}

	std::vector<double> z(N);
	for (size_t i = 1; i <= N; i++)
	{
		double Sum = 0;
		for (auto elem : Lhs[i-1])
		{
			if (elem.first < i - 1)
			{
				Sum += elem.second * z[elem.first];
			}			
		}
		z[i-1] = 1 / D.GetValue(i-1, i-1) * (Rhs[i-1] - Sum);
	}
	

	for (size_t i = N; i > 0; i--)
	{
		double Sum = 0;
		for (auto elem : Lhs[i - 1])
		{
			if (elem.first > i - 1)
			{
				Sum += elem.second * y[elem.first];
			}				
		}
		y[i-1] = z[i-1] - 1 / D.GetValue(i-1, i-1) * Sum;
	}
	return y;
}

std::vector<double> ISO0_2_P::Precondition(const CMatrix &Lhs, const std::vector<double> &Rhs)
{
	size_t N = Lhs.size();
	std::vector<double> y(N);
	CMatrix D(N);

	for (size_t i = 1; i <= N; i++)
	{	
		double Value = 0;
		for (auto elem : Lhs[i-1])
		{
			if (elem.first < i - 1)
			{
				double sum = 0;
				for (auto elem2 : Lhs[elem.first])
				{
					if (elem2.first > elem.first)
					{
						sum += elem2.second;
					}
				}
				
				Value += elem.second * elem.second / D.GetValue(elem.first, elem.first);
				//Value += Lhs.GetValue(elem.first, i - 1) * elem.second / D.GetValue(elem.first, elem.first);
			}
		}
		D.SetValue(i - 1,i - 1 , Lhs.GetValue(i-1,i-1) - Value);
	}

	std::vector<double> z(N);
	for (size_t i = 1; i <= N; i++)
	{
		double Sum = 0;
		for (auto elem : Lhs[i-1])
		{
			if (elem.first < i - 1)
			{
				Sum += elem.second * z[elem.first];
			}			
		}
		z[i-1] = 1 / D.GetValue(i-1, i-1) * (Rhs[i-1] - Sum);
	}
	

	for (size_t i = N; i > 0; i--)
	{
		double Sum = 0;
		for (auto elem : Lhs[i - 1])
		{
			if (elem.first > i - 1)
			{
				Sum += elem.second * y[elem.first];
			}				
		}
		y[i-1] = z[i-1] - 1 / D.GetValue(i-1, i-1) * Sum;
	}

	return y;	
}

// IPreconditioner* IPreconditioner::Fabric(const std::string &Precondition_method)
// {
// 	if (Precondition_method == "Jacobi")
// 	{
// 		return new Jacobi_P();
// 	}
// 	else if (Precondition_method == "SSOR")
// 	{
// 		return new SSOR_P();
// 	}
// 	else if (Precondition_method == "ISO0_1")
// 	{
// 		return new ISO0_1_P();
// 	}
// 	else if (Precondition_method == "ISO0_2")
// 	{
// 		return new ISO0_2_P();
// 	}
// 	return new Jacobi_P();
// }

// IPreconditioner* IPreconditioner::Fabric(const std::string &Precondition_method, const double omega)
// {
// 	SSOR_P* Preconditioner = new SSOR_P();
// 	Preconditioner->SetOmega(omega);
// 	return Preconditioner;
// }

// void SSOR_P::SetOmega(double omega)
// {
// 	this->omega = omega;
// }

// #endif