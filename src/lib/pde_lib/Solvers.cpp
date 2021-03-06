// #ifndef __MATRIX_SOLVERS_CPP__
// #define __MATRIX_SOLVERS_CPP__

#include "Solvers.h"
#include "VectorOperations.h"
#include "CompressedM.h"
#include "Preconditioners.h"
#include "SolverParams.h"


/*-----------------------------Realizations-----------------------------*/

// Hardcode
IMatrixSolver* IMatrixSolver::Factory(const MatrixSolverParams* parameters)
{
	if (parameters->solve_method == MatrixSolverParams::Thomas && parameters->precondition_method == MatrixSolverParams::Preconditioners::None)
	{	
		return new Thomas_Solver(parameters);
	} else if (parameters->solve_method == MatrixSolverParams::Methods::Jacobi && parameters->precondition_method == MatrixSolverParams::Preconditioners::None)
	{
		// std::shared_ptr<IMatrixSolver> solver(new Jacobi_Solver(parameters));
		return new Jacobi_Solver(parameters);
	} else if (parameters->solve_method == MatrixSolverParams::Methods::Seidel && parameters->precondition_method == MatrixSolverParams::Preconditioners::None)
	{
		// std::shared_ptr<IMatrixSolver> solver(new Seidel_Solver(parameters));
		return new Seidel_Solver(parameters);
	} else if (parameters->solve_method == MatrixSolverParams::Methods::SOR && parameters->precondition_method == MatrixSolverParams::Preconditioners::None)
	{
		// std::shared_ptr<IMatrixSolver> solver(new SOR_Solver(parameters));
		return new SOR_Solver(parameters);
	} else if (parameters->solve_method == MatrixSolverParams::Methods::SSOR && parameters->precondition_method == MatrixSolverParams::Preconditioners::None)
	{
		// std::shared_ptr<IMatrixSolver> solver(new SSOR_Solver(parameters));
		return new SSOR_Solver(parameters);
	} else if (parameters->solve_method == MatrixSolverParams::Methods::MR && parameters->precondition_method == MatrixSolverParams::Preconditioners::None)
	{
		// std::shared_ptr<IMatrixSolver> solver(new MR_Solver(parameters));
		return new MR_Solver(parameters);
	} else if (parameters->solve_method == MatrixSolverParams::Methods::GD && parameters->precondition_method == MatrixSolverParams::Preconditioners::None)
	{
		// std::shared_ptr<IMatrixSolver> solver(new GD_Solver(parameters));
		return new GD_Solver(parameters);
	} else if (parameters->solve_method == MatrixSolverParams::Methods::CG && parameters->precondition_method == MatrixSolverParams::Preconditioners::None)
	{
		// std::shared_ptr<IMatrixSolver> solver(new CG_Solver(parameters));
		return new CG_Solver(parameters);
	} else if (parameters->precondition_method != MatrixSolverParams::Preconditioners::None)
	{
		//std::shared_ptr<IPreconditioner> preconditioner = IPreconditioner::Factory(parameters);
		if (parameters->solve_method == MatrixSolverParams::Methods::GD)
		{
			//std::shared_ptr<IMatrixSolver> solver(new GD_Solver_P(parameters, preconditioner));
			return new GD_Solver_P(parameters, IPreconditioner::Factory(parameters));
		} else if (parameters->solve_method == MatrixSolverParams::Methods::CG) 
		{
			//std::shared_ptr<IMatrixSolver> solver(new CG_Solver_P(parameters, preconditioner));
			return new GD_Solver_P(parameters, IPreconditioner::Factory(parameters));
		}
	}
	return nullptr;
}

IMatrixSolver::IMatrixSolver(const MatrixSolverParams* parameters)
{
	this->parameters = parameters;
}

void CG_Solver::solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	double start_time = clock();
	//size_t N = Lhs.size();
	std::vector<double> r;
	std::vector<double> p;
	std::vector<double> R;

	r = vector_diff((Lhs * u), Rhs);
	p = r;
	double Norm = dot_product(r,p);
	std::vector<double> Ap = Lhs * p;
	double pAp = dot_product(p, Ap);
	double lambda = Norm / pAp;
	u = vector_diff(u, mult_n(p, lambda));

	for (size_t k = 1; k < parameters->MAX_ITERATIONS; k++)
	{
		r = vector_diff(r, mult_n(Ap, lambda));
		double NewNorm = dot_product(r,r);
		double alpha = NewNorm / Norm;
		p = vector_sum(r, mult_n(p, alpha));
		Ap = Lhs * p;
		pAp = dot_product(p, Ap);
		lambda = NewNorm / pAp;
		u = vector_diff(u, mult_n(p, lambda));
		Norm = NewNorm;
		if (k % parameters->Save_steps == 0)
		{
			R.push_back(sqrt(Norm));
		}
		if(sqrt(Norm) < parameters->eps * parameters->eps)
		{
			std::cout << "?????????????? ?????????????? ???? ????????????????: " << k <<std::endl;
			std::cout << "?????????? ?????????????????????? ???? ??????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
			std::cout << "?????????? ??????????????: " << sqrt(Norm) << std::endl;
			write_in_file(R, "Residuals");	
			break;
		}
	}
}

void GD_Solver::solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	double start_time = clock(); 				// ?????????? ??????????????
	std::vector<double> r;
	std::vector<double> R; 					// ?????? ???????????????????????? ?????????? ?????????????? ?? ??????????????????

	for (size_t k = 0; k < parameters->MAX_ITERATIONS; k++)
	{   
		r = vector_diff(Lhs * u, Rhs);
		double Norm = dot_product(r,r);
		double gamma = Norm / (dot_product(Lhs * r, r));

		u = vector_diff(u, mult_n(r, gamma));

		if(k % parameters->Save_steps == 0) 
		{
			R.push_back(sqrt(Norm));
		} 
		
		// ?????????????? ???????????? ???? ??????????
		if(sqrt(Norm) < (parameters->eps))
		{
			std::cout << "?????????????? ?????????????? ???? ????????????????: " << k <<std::endl;
			std::cout << "?????????? ?????????????????????? ???? ??????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << " ????????????" << std::endl;
			write_in_file(R, "Residuals");
			break;
		}
	} // k
}

void MR_Solver::solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
    double start_time = clock();
    std::vector<double> r(Rhs.size());
    // ?????? ???????????????????????? ?????????? ?????????????? ?? ??????????????????
    std::vector <double> R;

    for (size_t k = 0; k < parameters->MAX_ITERATIONS; k++)
    {
        r = vector_diff(Lhs * u, Rhs);

        std::vector<double> Ar = Lhs * r;
        double gamma = dot_product(r , Ar) / dot_product(Ar , Ar);
        u = vector_diff(u, mult_n(r, gamma));
        
        double Norm = sqrt(dot_product(r,r));

        if(k % parameters->Save_steps == 0)
        {
            R.push_back(Norm);
        } 

        if(Norm < parameters->eps)
        {
            std::cout << "?????????????? ?????????????? ???? ????????????????: " << k <<std::endl;
            std::cout << "?????????? ?????????????????????? ???? ??????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
			write_in_file(R, "Residuals");	
            break;
        }
    }
}

void SSOR_Solver::solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	size_t N = Lhs.size();
	std::vector<double> Temp_U(N, 0.);
	std::vector<double> r(N, 0.);

	double start_time = clock();
	for (size_t k = 0; k < parameters->MAX_ITERATIONS; k++)
	{
		// Forward Propogation
		for (size_t i = 0; i < N; i++)
		{
			double LowerSum = 0;
			double HigherSum = 0;

			for (auto el : Lhs[i])
			{
				if (el.first < i)
				{
					LowerSum += el.second * Temp_U[el.first];
				} // L * Temp_U
				if (el.first > i)
				{
					HigherSum += el.second * u[el.first];
				} // U * u
			} // End U * u & L * Temp_U // el
			
			Temp_U[i] = parameters->omega_solver / Lhs.GetValue(i,i) * ( - LowerSum - HigherSum + Rhs[i]) + (1 - parameters->omega_solver) * u[i];
		} // End i - filling solution vector // i

		// Back Propogation
		for (int i = N - 1; i >= 0; i--)
		{
			double LowerSum = 0;
			double HigherSum = 0;

			for (auto el : Lhs[i])
			{
				if (el.first < i)
				{
					LowerSum += el.second * Temp_U[el.first];
				} // L * Temp_U
				if (el.first > i)
				{
					HigherSum += el.second * u[el.first];
				} // U * u
			} // End U * u & L * Temp_U // el
			
			u[i] = parameters->omega_solver / Lhs.GetValue(i,i) * ( - LowerSum - HigherSum + Rhs[i]) + (1 - parameters->omega_solver) * Temp_U[i];
		} // End i - filling solution vector // i

		double res = norm_2(vector_diff(Lhs * u, Rhs));

		if (k % parameters->Save_steps == 0)
		{
			R.push_back(res);
		}

		if (res< parameters->eps)
		{
			std::cout << "?????????????? ?????????????? ???? ????????????????: " << k <<std::endl;
			std::cout << "?????????? ?????????????????????? ???? ??????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << " ????????????" << std::endl;
			write_in_file(R, "Residuals");
			break;
		}		
	} // Iterations cycle // k
}

void SOR_Solver::solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	size_t N = Lhs.size();
	std::vector<double> Temp_U(N, 0.);
	std::vector<double> r(N, 0.);

	double start_time = clock();
	for (size_t k = 0; k < parameters->MAX_ITERATIONS; k++)
	{
		for (size_t i = 0; i < N; i++)
		{
			double LowerSum = 0;
			double HigherSum = 0;

			for (auto el : Lhs[i])
			{
				if (el.first < i)
				{
					LowerSum += el.second * Temp_U[el.first];
				} // L * Temp_U
				if (el.first > i)
				{
					HigherSum += el.second * u[el.first];
				} // U * u
			} // End U * u & L * Temp_U // el
			
			Temp_U[i] = parameters->omega_solver / Lhs.GetValue(i,i) * ( - LowerSum - HigherSum + Rhs[i]) + (1 - parameters->omega_solver) * u[i];
		} // End i - filling solution vector // i


		double res = norm_2(vector_diff(Lhs * u, Rhs));

		if (k % parameters->Save_steps == 0)
		{
			R.push_back(res);
		}

		if (res < parameters->eps)
		{
			std::cout << "?????????????? ?????????????? ???? ????????????????: " << k <<std::endl;
			std::cout << "?????????? ?????????????????????? ???? ??????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << " ????????????" << std::endl;
			write_in_file(R, "Residuals");
			break;
		}

		u = Temp_U;		
	} // Iterations cycle // k
}

void Thomas_Solver::solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	if (CheckTridiagonal(Lhs) == false)
	{
		std::cout << "?????????????? ???? ???????????????? ????????????????????????????????" << std::endl;
	}
	double start_time = clock();
	size_t N = Lhs.size();
	// ?????????????? ???????????????????? ?????????? ?????????????? ???????? ???????????? ????????????????
	std::vector<double> P(N, 0.);
	// ?????????????? ???????????????????? ?????????? ?????????????? ???????? ???????????? ????????????????
	std::vector<double> Q(N, 0.);
	P[0] = -Lhs.GetValue(0,1)/(Lhs.GetValue(0,0));
	Q[0] = Rhs[0]/(Lhs.GetValue(0,0));
	for (size_t i = 1; i < N; i++)
	{
		P[i]=-Lhs.GetValue(i,i+1)/(Lhs.GetValue(i,i)+Lhs.GetValue(i,i-1)*P[i-1]);
		Q[i]=(-Lhs.GetValue(i,i-1)*Q[i-1]+Rhs[i])/(Lhs.GetValue(i,i)+Lhs.GetValue(i,i-1)*P[i-1]);
	}
	u[N-1] = (-Lhs.GetValue(N-1,N-2) * Q[N-2] + Rhs[N-1])/(Lhs.GetValue(N-1,N-1) + P[N-2] * Lhs.GetValue(N-1, N-2));
	for (size_t i = N-1; i > 0; i--)
	{
		u[i-1] = P[i-1] * (u[i]) + Q[i-1];
	}
	std::cout << "?????????? ?????????????????????? ???? ??????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << " ????????????" << std::endl;
}

bool Thomas_Solver::CheckTridiagonal(const CMatrix& Lhs)
{
	bool check_result = true;
	for (size_t i = 0; i < Lhs.size(); i++)
	{
		Lhs[i].size() <= 3 ? check_result = true : check_result = false;
	}
	return check_result;
}

void Jacobi_Solver::solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	double start_time = clock();
	std::size_t N = Lhs.size();
	std::vector<double> u_temp(N, 0.);

	for (size_t k = 0; k < parameters->MAX_ITERATIONS; k++)
	{
		for (size_t i = 0; i < N; i++)
		{
			double InnerSum = 0;
			for (auto elem : Lhs[i])
			{
				if (elem.first == i)
				{
					continue;
				}
				InnerSum += elem.second * u[elem.first];
			}

			u[i] = (-InnerSum + Rhs[i]) / Lhs.GetValue(i,i);
		}

		double res = norm_2(vector_diff(Lhs * u, Rhs));

		if (k % parameters->Save_steps == 0)
		{
			R.push_back(res);
		}
		

		if (res < parameters->eps)
		{
			std::cout << "?????????????? ?????????????? ???? ????????????????: " << k <<std::endl;
			std::cout << "?????????? ?????????????????????? ???? ??????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << " ????????????" << std::endl;
			write_in_file(R, "Residuals");
			break;
		}
	}
}

void Seidel_Solver::solve(const CMatrix& Lhs, const std::vector<double>& Rhs, std::vector<double>& u)
{
	size_t N = Lhs.size();
	std::vector <double> U(N, 0.);
	double start_time = clock();

	for (size_t k = 0; k < parameters->MAX_ITERATIONS; k++)
	{
		for (size_t i = 0; i < N; i++)
		{
			double LowerSum = 0;
			double HigherSum = 0;

			for (auto elem : Lhs[i])
			{
				if (elem.first < i)
				{
				  	LowerSum += elem.second * U[elem.first];
				}
				if (elem.first > i)
				{
					HigherSum += elem.second * u[elem.first];
				}
			}
			
			U[i] = (- HigherSum - LowerSum + Rhs[i]) / Lhs.GetValue(i,i);
		}


		double res = norm_2(vector_diff(Lhs * u, Rhs));

		if (k % parameters->Save_steps == 0)
		{
			R.push_back(res);
		}

		if (res < parameters->eps)
		{
			std::cout << "?????????????? ?????????????? ???? ????????????????: " << k <<std::endl;
			std::cout << "?????????? ?????????????????????? ???? ??????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << " ????????????" << std::endl;
			write_in_file(R, "Residuals");
			break;
		}
		u = U;
	}
	
}

void LU_Solver::LU_decomposition(const CMatrix &Lhs, CMatrix &L)
{
	std::size_t N = Lhs.size();
	double start_time = clock();

	// ?????? ?????????????????? ?????? ???????????????????? L(ower)
	for (size_t i = 0; i < N; i++) 
	{
		for (size_t j = 0; j <= i; j++) 
		{
			double sum = 0;
			for (size_t k = 0; k < j; k++)
			{
				sum += L.GetValue(i,k) * L.GetValue(j,k);
			}
			if (i == j)
			{
				L.SetValue(i,i, sqrt(Lhs.GetValue(i,i) - sum));
			}
			else
			{
				double Diff = Lhs.GetValue(i,j) - sum;
				if (Diff != 0)
				{
					L.SetValue(i,j, (1.0 / L.GetValue(j,j) * (Diff)));
				}
			}
		}
	}
	std::cout << "?????????????????? LU ???????????????????? ????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
}

void LU_Solver::solve(const CMatrix &Lhs, const std::vector<double> &Rhs, std::vector<double> &u)
{
	double start_time = clock();
	size_t N = Lhs.size();
	CMatrix L(N);
	LU_decomposition(Lhs, L);
	double start_clock_Solution = clock();
	std::vector<double> y(N);

	for (size_t i = 1; i <= N; i++)
	{
		double Sum = 0;
		for (size_t j = 1; j <= i - 1; j++)
		{
			Sum += L.GetValue(i-1, j-1) * y[j-1];
			//Sum += L[i-1][j-1] * y[j-1];
		}
		y[i-1] = 1 / L.GetValue(i-1,i-1) * (Rhs[i-1] - Sum);
		//y[i-1] = 1 / L[i-1][i-1] * (RHS[i-1] - Sum);
	}

	for (size_t i = N; i > 0; i--)
	{
		double Sum = 0;
		for (size_t j = i+1; j < N + 1; j++)
		{
			Sum += L.GetValue(j-1,i-1) * u[j-1];
			//Sum += L[j-1][i-1] * x[j-1];
		}
		u[i-1] = 1 / L.GetValue(i-1,i-1) * (y[i-1] - Sum);
		//x[i-1] = 1 / L[i-1][i-1] * (y[i-1] - Sum);
	}

	std::cout << "?????????????????? ?????????????? ???????? ?? ?????????????????????????? LU ???????????????????? ????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
	std::cout << "?????????? ?????????????? ?????? ?????????? ?????????????????? LU ????????????????????: " << (clock() - start_clock_Solution) / CLOCKS_PER_SEC << std::endl;
}

void LDU_Solver::LDU_decomposition(const CMatrix &Lhs, CMatrix &L, CMatrix &D)
{
	double start_time = clock();
	// ?????? ?????????????????? ?????? ???????????????????? D(iagonal) ?? L(ower)
	for (size_t i = 0; i < Lhs.size(); i++)
	{
		for (size_t j = 0; j <= i; j++)
		{
			double sum = 0;
			for (size_t k = 0; k < j; k++)
			{
				sum += L.GetValue(i,k) * L.GetValue(j,k) * D.GetValue(k,k);
			}
			D.SetValue(i,i, Lhs.GetValue(i,i) - sum);
			L.SetValue(i,i,1);
			double diff = (Lhs.GetValue(i,j) - sum);
			if (diff != 0)
			{
				L.SetValue(i,j, 1 / D.GetValue(j,j) * diff);
			}
		}
	}
	std::cout << "?????????????????? LDU ???????????????????? ????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
}

void LDU_Solver::solve(const CMatrix &Lhs, const std::vector<double> &Rhs, std::vector<double> &u)
{
	double start_time = clock();
	size_t N = Lhs.size();
	CMatrix L(N);
	CMatrix D(N);

	LDU_decomposition(Lhs, L, D);

	std::vector<double> y(N);
	std::vector<double> x(N);

	for (size_t i = 1; i <= N; i++)
	{
		double Sum = 0;
		for (size_t j = 1; j <= i - 1; j++)
		{
			Sum += D.GetValue(j-1,j-1) * L.GetValue(i-1,j-1) * y[j-1];
		}
		y[i-1] = 1 / D.GetValue(i-1,i-1) * (Rhs[i-1] - Sum);
	}

	for (size_t i = N; i > 0; i--)
	{
		double Sum = 0;
		for (size_t j = i+1; j < N + 1; j++)
		{
			Sum += L.GetValue(j-1,i-1) * u[j-1];
		}
		u[i-1] = 1 / L.GetValue(i-1, i-1) * (y[i-1] - Sum);
	}

	std::cout << "?????????????????? ?????????????? ???????? ?? ?????????????????????????? LDU ???????????????????? ????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
}

void CG_Solver_P::solve(const CMatrix &Lhs, const std::vector<double> &Rhs, std::vector<double> &u)
{
	double start_time = clock();
	size_t N = Lhs.size();
	std::vector<double> w;
	std::vector<double> r = vector_diff(Lhs * u, Rhs);
	std::vector<double> p = r;
	std::vector<double> Ap = Lhs * p;
	double rw = dot_product(r,r);
	double lambda = (dot_product(p,r)) / (dot_product(p, Ap));
	u = vector_diff(u, mult_n(p, lambda));

	for (size_t k = 0; k < parameters->MAX_ITERATIONS; k++)
	{	
		r = vector_diff(r, mult_n(Ap, lambda));
		w = Preconditioner->Precondition(Lhs, r);
		double RW = dot_product(r, w);
		double Norm2 = dot_product(r, r);
		double alpha = RW * (1 / rw);
		rw = RW;
		p = vector_sum(w, mult_n(p, alpha));
		Ap = Lhs * p;
		double pAp = dot_product(p, Ap);
		lambda = rw / pAp;
		u = vector_diff(u, mult_n(p, lambda));

        if(k % parameters->Save_steps == 0)
        {
            //Y_iter.push_back(Y);
            R.push_back(sqrt(Norm2));
        } 

		if(sqrt(Norm2) < parameters->eps * parameters->eps)
        {
            std::cout << "?????????????? ?????????????? ???? ????????????????: " << k <<std::endl;
            std::cout << "?????????? ?????????????????????? ???? ??????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
			write_in_file(R, "Residuals_P");
            break;
        }
	}
}

void GD_Solver_P::solve(const CMatrix &Lhs, const std::vector<double> &Rhs, std::vector<double> &u)
{
	//! ???????????? ?? ?????????????????? ????????????
	double start_time = clock();
    std::vector<double> r;
	std::vector<double> w;
    // ?????? ???????????????????????? ?????????? ?????????????? ?? ??????????????????
    std::vector <double> R;

    for (size_t k = 0; k < parameters->MAX_ITERATIONS; k++)
    {   
        // r = Ay-B
        r = vector_diff(Lhs * u, Rhs);
		w = Preconditioner->Precondition(Lhs, r);
        // gamma = (r,r) / (Ar,r)
		double Norm2 = dot_product(r,r);
        double gamma = dot_product(r, w) / dot_product(r, Lhs * w);

        u = vector_diff(u, mult_n(r, gamma));

        if(k % parameters->Save_steps == 0)
        {
            //Y_iter.push_back(Y);
            R.push_back(sqrt(Norm2));
        } 

        if(sqrt(Norm2) < parameters->eps)
        {
            std::cout << "?????????????? ?????????????? ???? ????????????????: " << k <<std::endl;
            std::cout << "?????????? ?????????????????????? ???? ??????????????: " << (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
			write_in_file(R, "Residuals_P");
            break;
        }
    }
}

CG_Solver_P::CG_Solver_P(const MatrixSolverParams* parameters, IPreconditioner* Preconditioner) : IMatrixSolver(parameters)
{
	this->Preconditioner = Preconditioner;
}

GD_Solver_P::GD_Solver_P(const MatrixSolverParams* parameters, IPreconditioner* Preconditioner) : IMatrixSolver(parameters)
{
	this->Preconditioner = Preconditioner;
}

// IMatrixSolver* IMatrixSolver::Fabric(const std::string& method_name)
// {
// 	if (method_name == "GD")
// 	{
// 		return new GD_Solver();
// 	} 
// 	else if (method_name == "MR")
// 	{
// 		return new MR_Solver();
// 	} 
// 	else if (method_name == "CG")
// 	{
// 		return new CG_Solver();
// 	} 
// 	else if (method_name == "SSOR")
// 	{
// 		return new SSOR_Solver();
// 	} 
// 	else if (method_name == "SOR")
// 	{
// 		return new SOR_Solver();
// 	} 
// 	else if (method_name == "Seidel")
// 	{
// 		return new Seidel_Solver();
// 	}
// 	else if (method_name == "Jacobi")
// 	{
// 		return new Jacobi_Solver();
// 	}
// 	else if (method_name == "LU")
// 	{
// 		return new LU_Solver();
// 	}
// 	else if (method_name == "LDU")
// 	{
// 		return new LDU_Solver();
// 	}
// 	return new Thomas_Solver();
// }

// IMatrixSolver* IMatrixSolver::Fabric(const std::string &method_name, const double omega)
// {
// 	if (method_name == "SSOR")
// 	{
// 		SSOR_Solver* solver = new SSOR_Solver();
// 		solver->SetOmega(omega);
// 		return solver;
// 	} else if (method_name == "SOR")
// 	{
// 		SOR_Solver* solver = new SOR_Solver();
// 		solver->SetOmega(omega);
// 		return solver;
// 	}
// 	SOR_Solver* solver = new SOR_Solver();
// 	solver->SetOmega(omega);
// 	return solver;
// }

// IMatrixSolver* IMatrixSolver::Fabric_P(const std::string &method_name, const std::string &Precondition_method)
// {
// 	IPreconditioner* Preconditioner = IPreconditioner::Fabric(Precondition_method);
// 	if (method_name == "CG")
// 	{
// 		return new CG_Solver_P(Preconditioner);
// 	}
// 	else if(method_name == "GD")
// 	{
// 		return new GD_Solver_P(Preconditioner);
// 	}
// 	return new GD_Solver_P(Preconditioner);
// }

// IMatrixSolver* IMatrixSolver::Fabric_P(const std::string &method_name, const std::string &Precondition_method, const double omega)
// {
// 	IPreconditioner* Preconditioner = IPreconditioner::Fabric(Precondition_method, omega);
// 	if (method_name == "CG")
// 	{
// 		return new CG_Solver_P(Preconditioner);
// 	}
// 	else if(method_name == "GD")
// 	{
// 		return new GD_Solver_P(Preconditioner);
// 	}
// 	return new GD_Solver_P(Preconditioner);
// }

// void SolversParams_SOR_SSOR::set_params( MatrixSolverParams::Methods method = Thomas,
// 	MatrixSolverParams::Preconditioners precondition_method = None,
// 	size_t MAX_ITERATIONS = 10000000,
// 	double eps = 1.e-5,
// 	size_t Save_steps = 10,
// 	double omega = 1.95)
// {
// 	solve_method = method;
// 	precondition_method = precondition_method;
// 	MAX_ITERATIONS = MAX_ITERATIONS;
// 	eps = eps;
// 	Save_steps = Save_steps;
// 	omega = omega;
// }

// #endif