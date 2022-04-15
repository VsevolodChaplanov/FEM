#include "PDE_lib_tests.h"

// Exist test cases
// Solve methods of SLAE 
// 		Thomas method
//			Tridiagonal check


// ------------------ Thomas solver test ------------------
TEST_CASE("Thomas solver check", "ThomasSolver")
{
	CMatrix lhs_ex(6);

	lhs_ex.SetValue(0, 0, 2);
	lhs_ex.SetValue(0, 1, 1);
	lhs_ex.SetValue(1, 0, 1);
	lhs_ex.SetValue(1, 1, 2);
	lhs_ex.SetValue(1, 2, 1);
	lhs_ex.SetValue(2, 1, 1);
	lhs_ex.SetValue(2, 2, 2);
	lhs_ex.SetValue(3, 2, 1);
	lhs_ex.SetValue(3, 3, 2);
	lhs_ex.SetValue(3, 4, 1);
	lhs_ex.SetValue(4, 3, 1);
	lhs_ex.SetValue(4, 4, 2);
	lhs_ex.SetValue(4, 5, 1);
	lhs_ex.SetValue(5, 4, 1);
	lhs_ex.SetValue(5, 5, 2);

	// ----------------- Tridiagonal check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		for (auto elem : lhs_ex[i])
		{
			INFO( "i = " << i);
			bool res = elem.first == i || elem.first == i - 1 || elem.first == i + 1;
			REQUIRE( res == true);
		}
	}
	// ----------------- Tridiagonal check ----------------- //


	std::vector<double> solution_ex(6, 1.);
	std::vector<double> rhs_ex = lhs_ex * solution_ex;

	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::Thomas, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10);
	IMatrixSolver* thomas_solver = IMatrixSolver::Factory(params);

	std::vector<double> test_num_solution(6);
	thomas_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		REQUIRE( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //

	delete params;
	delete thomas_solver;
}