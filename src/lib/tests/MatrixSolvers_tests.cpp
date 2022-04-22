#include "PDE_lib_tests.h"

// Exist test cases
// Solve methods of SLAE 
// 		Thomas method
//			Tridiagonal check
//		Jacobi method
//		SOR method
//		SSOR method
//		Gradient descent
//			Symm check
//		Conjugate gradients
//			Symm check
//		LU (Cholesky) factorization
//			Symm check
//		LDU factorization
//			Symm check


// ------------------ Thomas solver test ------------------
TEST_CASE("Thomas solver check", "[ThomasSolver][ExactMethods]")
{
	CMatrix lhs_ex(6);

	lhs_ex.SetValue(0, 0, 2);
	lhs_ex.SetValue(0, 1, 1);
	lhs_ex.SetValue(1, 0, 1);
	lhs_ex.SetValue(1, 1, 2);
	lhs_ex.SetValue(1, 2, 1);
	lhs_ex.SetValue(2, 1, 1);
	lhs_ex.SetValue(2, 2, 2);
	lhs_ex.SetValue(2, 3, 1);
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

// ------------------ Jacobi solver test ------------------
TEST_CASE( "Jacobi solver check", "[JacobiSolver]" )
{
	CMatrix lhs_ex(6);

	// -------------- Tridiagonal part -------------- //
	lhs_ex.SetValue(0, 0, 2);
	lhs_ex.SetValue(0, 1, 1);
	lhs_ex.SetValue(1, 0, 1);
	lhs_ex.SetValue(1, 1, 2);
	lhs_ex.SetValue(1, 2, 1);
	lhs_ex.SetValue(2, 1, 1);
	lhs_ex.SetValue(2, 2, 2);
	lhs_ex.SetValue(2, 3, 1);
	lhs_ex.SetValue(3, 2, 1);
	lhs_ex.SetValue(3, 3, 2);
	lhs_ex.SetValue(3, 4, 1);
	lhs_ex.SetValue(4, 3, 1);
	lhs_ex.SetValue(4, 4, 2);
	lhs_ex.SetValue(4, 5, 1);
	lhs_ex.SetValue(5, 4, 1);
	lhs_ex.SetValue(5, 5, 2);
	// -------------- Tridiagonal part -------------- //

	// -------------- Outer Diagonal elements -------------- //
	lhs_ex.SetValue(5,0, 1);
	lhs_ex.SetValue(0,5, 1);
	lhs_ex.SetValue(3, 5, 1);
	lhs_ex.SetValue(5, 3, 1);
	// -------------- Outer Diagonal elements -------------- //

	std::vector<double> solution_ex(6, 1.);
	std::vector<double> test_num_solution(6);
	std::vector<double> rhs_ex = lhs_ex * solution_ex;

	// ----------------------------------------------------------------------- //
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::Jacobi, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10);
	IMatrixSolver* jacobi_solver = IMatrixSolver::Factory(params);

	jacobi_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------------------------------------------------------------- //

	delete params;
	delete jacobi_solver;
}


// ------------------ Seidel solver test ------------------
TEST_CASE( "Seidel solver check", "[SeidelSolver]" )
{
	CMatrix lhs_ex(6);

	// -------------- Tridiagonal part -------------- //
	lhs_ex.SetValue(0, 0, 2);
	lhs_ex.SetValue(0, 1, 1);
	lhs_ex.SetValue(1, 0, 1);
	lhs_ex.SetValue(1, 1, 2);
	lhs_ex.SetValue(1, 2, 1);
	lhs_ex.SetValue(2, 1, 1);
	lhs_ex.SetValue(2, 2, 2);
	lhs_ex.SetValue(2, 3, 1);
	lhs_ex.SetValue(3, 2, 1);
	lhs_ex.SetValue(3, 3, 2);
	lhs_ex.SetValue(3, 4, 1);
	lhs_ex.SetValue(4, 3, 1);
	lhs_ex.SetValue(4, 4, 2);
	lhs_ex.SetValue(4, 5, 1);
	lhs_ex.SetValue(5, 4, 1);
	lhs_ex.SetValue(5, 5, 2);
	// -------------- Tridiagonal part -------------- //

	// -------------- Outer Diagonal elements -------------- //
	lhs_ex.SetValue(5,0, 1);
	lhs_ex.SetValue(0,5, 1);
	lhs_ex.SetValue(3, 5, 1);
	lhs_ex.SetValue(5, 3, 1);
	// -------------- Outer Diagonal elements -------------- //

	std::vector<double> solution_ex(6, 1.);
	std::vector<double> test_num_solution(6);
	std::vector<double> rhs_ex = lhs_ex * solution_ex;

	// ----------------------------------------------------------------------- //
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::Seidel, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10);
	IMatrixSolver* seidel_solver = IMatrixSolver::Factory(params);

	seidel_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------------------------------------------------------------- //

	delete params;
	delete seidel_solver;
}

// ------------------ SOR solver test ------------------
TEST_CASE( "SOR solver check", "[SORSolver]" )
{
	CMatrix lhs_ex(6);

	// -------------- Tridiagonal part -------------- //
	lhs_ex.SetValue(0, 0, 2);
	lhs_ex.SetValue(0, 1, 1);
	lhs_ex.SetValue(1, 0, 1);
	lhs_ex.SetValue(1, 1, 2);
	lhs_ex.SetValue(1, 2, 1);
	lhs_ex.SetValue(2, 1, 1);
	lhs_ex.SetValue(2, 2, 2);
	lhs_ex.SetValue(2, 3, 1);
	lhs_ex.SetValue(3, 2, 1);
	lhs_ex.SetValue(3, 3, 2);
	lhs_ex.SetValue(3, 4, 1);
	lhs_ex.SetValue(4, 3, 1);
	lhs_ex.SetValue(4, 4, 2);
	lhs_ex.SetValue(4, 5, 1);
	lhs_ex.SetValue(5, 4, 1);
	lhs_ex.SetValue(5, 5, 2);
	// -------------- Tridiagonal part -------------- //

	// -------------- Outer Diagonal elements -------------- //
	lhs_ex.SetValue(5,0, 1);
	lhs_ex.SetValue(0,5, 1);
	lhs_ex.SetValue(3, 5, 1);
	lhs_ex.SetValue(5, 3, 1);
	// -------------- Outer Diagonal elements -------------- //

	std::vector<double> solution_ex(6, 1.);
	std::vector<double> test_num_solution(6);
	std::vector<double> rhs_ex = lhs_ex * solution_ex;

	// ----------------------------------------------------------------------- //
	// ------------------ Relaxation parameter omega = 1.95 ------------------ //
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::SOR, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10, 1.95);
	IMatrixSolver* sor_solver = IMatrixSolver::Factory(params);

	sor_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------- Relaxation parameter omega = 1.95 ----------------- //
	// ----------------------------------------------------------------------- //


	// ----------------------------------------------------------------------- //
	// ------------------ Relaxation parameter omega = 1.5 ------------------ //
	params = new MatrixSolverParams(MatrixSolverParams::Methods::SOR, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10, 1.5);
	sor_solver = IMatrixSolver::Factory(params);

	sor_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------- Relaxation parameter omega = 1.5 ----------------- //
	// ----------------------------------------------------------------------- //

	// ----------------------------------------------------------------------- //
	// ------------------ Relaxation parameter omega = 1.99 ------------------ //
	params = new MatrixSolverParams(MatrixSolverParams::Methods::SSOR, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10, 1.99);
	sor_solver = IMatrixSolver::Factory(params);

	sor_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------- Relaxation parameter omega = 1.99 ----------------- //
	// ----------------------------------------------------------------------- //

	delete params;
	delete sor_solver;
}

// ------------------ SSOR solver test ------------------
TEST_CASE( "SSOR solver check", "[SSORSolver]" )
{
	CMatrix lhs_ex(6);

	// -------------- Tridiagonal part -------------- //
	lhs_ex.SetValue(0, 0, 2);
	lhs_ex.SetValue(0, 1, 1);
	lhs_ex.SetValue(1, 0, 1);
	lhs_ex.SetValue(1, 1, 2);
	lhs_ex.SetValue(1, 2, 1);
	lhs_ex.SetValue(2, 1, 1);
	lhs_ex.SetValue(2, 2, 2);
	lhs_ex.SetValue(2, 3, 1);
	lhs_ex.SetValue(3, 2, 1);
	lhs_ex.SetValue(3, 3, 2);
	lhs_ex.SetValue(3, 4, 1);
	lhs_ex.SetValue(4, 3, 1);
	lhs_ex.SetValue(4, 4, 2);
	lhs_ex.SetValue(4, 5, 1);
	lhs_ex.SetValue(5, 4, 1);
	lhs_ex.SetValue(5, 5, 2);
	// -------------- Tridiagonal part -------------- //

	// -------------- Outer Diagonal elements -------------- //
	lhs_ex.SetValue(5,0, 1);
	lhs_ex.SetValue(0,5, 1);
	lhs_ex.SetValue(3, 5, 1);
	lhs_ex.SetValue(5, 3, 1);
	// -------------- Outer Diagonal elements -------------- //

	std::vector<double> solution_ex(6, 1.);
	std::vector<double> test_num_solution(6);
	std::vector<double> rhs_ex = lhs_ex * solution_ex;

	// ----------------------------------------------------------------------- //
	// ------------------ Relaxation parameter omega = 1.95 ------------------ //
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::SSOR, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10, 1.95);
	IMatrixSolver* ssor_solver = IMatrixSolver::Factory(params);

	ssor_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------- Relaxation parameter omega = 1.95 ----------------- //
	// ----------------------------------------------------------------------- //


	// ----------------------------------------------------------------------- //
	// ------------------ Relaxation parameter omega = 1.5 ------------------ //
	params = new MatrixSolverParams(MatrixSolverParams::Methods::SSOR, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10, 1.5);
	ssor_solver = IMatrixSolver::Factory(params);

	ssor_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------- Relaxation parameter omega = 1.5 ----------------- //
	// ----------------------------------------------------------------------- //

	// ----------------------------------------------------------------------- //
	// ------------------ Relaxation parameter omega = 1.99 ------------------ //
	params = new MatrixSolverParams(MatrixSolverParams::Methods::SSOR, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10, 1.99);
	ssor_solver = IMatrixSolver::Factory(params);

	ssor_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------- Relaxation parameter omega = 1.99 ----------------- //
	// ----------------------------------------------------------------------- //

	delete params;
	delete ssor_solver;
}

TEST_CASE( "Gradient Descent SLAE solver check", "[GDSolver]" )
{
	CMatrix lhs_ex(6);

	// -------------- Tridiagonal part -------------- //
	lhs_ex.SetValue(0, 0, 2);
	lhs_ex.SetValue(0, 1, 1);
	lhs_ex.SetValue(1, 0, 1);
	lhs_ex.SetValue(1, 1, 2);
	lhs_ex.SetValue(1, 2, 1);
	lhs_ex.SetValue(2, 1, 1);
	lhs_ex.SetValue(2, 2, 2);
	lhs_ex.SetValue(2, 3, 1);
	lhs_ex.SetValue(3, 2, 1);
	lhs_ex.SetValue(3, 3, 2);
	lhs_ex.SetValue(3, 4, 1);
	lhs_ex.SetValue(4, 3, 1);
	lhs_ex.SetValue(4, 4, 2);
	lhs_ex.SetValue(4, 5, 1);
	lhs_ex.SetValue(5, 4, 1);
	lhs_ex.SetValue(5, 5, 2);
	// -------------- Tridiagonal part -------------- //

	// -------------- Outer Diagonal elements -------------- //
	lhs_ex.SetValue(5,0, 1);
	lhs_ex.SetValue(0,5, 1);
	lhs_ex.SetValue(3, 5, 1);
	lhs_ex.SetValue(5, 3, 1);
	// -------------- Outer Diagonal elements -------------- //

	std::vector<double> solution_ex(6, 1.);
	std::vector<double> test_num_solution(6);
	std::vector<double> rhs_ex = lhs_ex * solution_ex;


	// The method converges well enough for non-symmetric matrices as well
	// -------------- Symmetry check -------------- //
	for (size_t i = 0; i < 6; i++)
	{
		for (size_t j = 0; j < 6; j++)
		{
			CHECK( lhs_ex.GetValue(i, j) == lhs_ex.GetValue(j, i));
		}			
	}
	// -------------- Symmetry check -------------- //

	// ----------------------------------------------------------------------- //
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::GD, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10);
	IMatrixSolver* gd_solver = IMatrixSolver::Factory(params);

	gd_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------------------------------------------------------------- //

	delete params;
	delete gd_solver;
}

TEST_CASE( "Conjugate Gradients SLAE solver check", "[CGSolver][ExactMethods]" )
{
	CMatrix lhs_ex(6);

	// -------------- Tridiagonal part -------------- //
	lhs_ex.SetValue(0, 0, 2);
	lhs_ex.SetValue(0, 1, 1);
	lhs_ex.SetValue(1, 0, 1);
	lhs_ex.SetValue(1, 1, 2);
	lhs_ex.SetValue(1, 2, 1);
	lhs_ex.SetValue(2, 1, 1);
	lhs_ex.SetValue(2, 2, 2);
	lhs_ex.SetValue(2, 3, 1);
	lhs_ex.SetValue(3, 2, 1);
	lhs_ex.SetValue(3, 3, 2);
	lhs_ex.SetValue(3, 4, 1);
	lhs_ex.SetValue(4, 3, 1);
	lhs_ex.SetValue(4, 4, 2);
	lhs_ex.SetValue(4, 5, 1);
	lhs_ex.SetValue(5, 4, 1);
	lhs_ex.SetValue(5, 5, 2);
	// -------------- Tridiagonal part -------------- //

	// -------------- Outer Diagonal elements -------------- //
	lhs_ex.SetValue(5, 0, 1);
	lhs_ex.SetValue(0, 5, 1);
	lhs_ex.SetValue(3, 5, 1);
	lhs_ex.SetValue(5, 3, 1);
	// -------------- Outer Diagonal elements -------------- //

	std::vector<double> solution_ex(6, 1.);
	std::vector<double> test_num_solution(6);
	std::vector<double> rhs_ex = lhs_ex * solution_ex;


	// The method converges well enough for non-symmetric matrices as well
	// -------------- Symmetry check -------------- //
	for (size_t i = 0; i < 6; i++)
	{
		for (size_t j = 0; j < 6; j++)
		{
			CHECK( lhs_ex.GetValue(i, j) == lhs_ex.GetValue(j, i));
		}			
	}
	// -------------- Symmetry check -------------- //

	// ----------------------------------------------------------------------- //
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::CG, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10);
	IMatrixSolver* cg_solver = IMatrixSolver::Factory(params);

	cg_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// For symmtric matrices method has limit in iterations as:
	// iterations_reqiured <= N (matrix NxN)
	REQUIRE( cg_solver->get_iterations_number() <= 6 );

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------------------------------------------------------------- //

	delete params;
	delete cg_solver;
}

TEST_CASE( "Cholesky decomposition SLAE solver check", "[LUSolver][ExactMethods]" )
{
	CMatrix lhs_ex(6);

	// -------------- Tridiagonal part -------------- //
	lhs_ex.SetValue(0, 0, 2);
	lhs_ex.SetValue(0, 1, 1);
	lhs_ex.SetValue(1, 0, 1);
	lhs_ex.SetValue(1, 1, 2);
	lhs_ex.SetValue(1, 2, 1);
	lhs_ex.SetValue(2, 1, 1);
	lhs_ex.SetValue(2, 2, 2);
	lhs_ex.SetValue(2, 3, 1);
	lhs_ex.SetValue(3, 2, 1);
	lhs_ex.SetValue(3, 3, 2);
	lhs_ex.SetValue(3, 4, 1);
	lhs_ex.SetValue(4, 3, 1);
	lhs_ex.SetValue(4, 4, 2);
	lhs_ex.SetValue(4, 5, 1);
	lhs_ex.SetValue(5, 4, 1);
	lhs_ex.SetValue(5, 5, 2);
	// -------------- Tridiagonal part -------------- //

	// -------------- Outer Diagonal elements -------------- //
	lhs_ex.SetValue(5, 0, 1);
	lhs_ex.SetValue(0, 5, 1);
	lhs_ex.SetValue(3, 5, 1);
	lhs_ex.SetValue(5, 3, 1);
	// -------------- Outer Diagonal elements -------------- //

	std::vector<double> solution_ex(6, 1.);
	std::vector<double> test_num_solution(6);
	std::vector<double> rhs_ex = lhs_ex * solution_ex;


	// The method converges well enough for non-symmetric matrices as well
	// -------------- Symmetry check -------------- //
	for (size_t i = 0; i < 6; i++)
	{
		for (size_t j = 0; j < 6; j++)
		{
			CHECK( lhs_ex.GetValue(i, j) == lhs_ex.GetValue(j, i));
		}			
	}
	// -------------- Symmetry check -------------- //

	// ----------------------------------------------------------------------- //
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::LU, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10);
	IMatrixSolver* lu_solver = IMatrixSolver::Factory(params);

	lu_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------------------------------------------------------------- //

	delete params;
	delete lu_solver;
}

TEST_CASE( "LDU decomposition SLAE solver check", "[LDUSolver][ExactMethods]" )
{
	CMatrix lhs_ex(6);

	// -------------- Tridiagonal part -------------- //
	lhs_ex.SetValue(0, 0, 2);
	lhs_ex.SetValue(0, 1, 1);
	lhs_ex.SetValue(1, 0, 1);
	lhs_ex.SetValue(1, 1, 2);
	lhs_ex.SetValue(1, 2, 1);
	lhs_ex.SetValue(2, 1, 1);
	lhs_ex.SetValue(2, 2, 2);
	lhs_ex.SetValue(2, 3, 1);
	lhs_ex.SetValue(3, 2, 1);
	lhs_ex.SetValue(3, 3, 2);
	lhs_ex.SetValue(3, 4, 1);
	lhs_ex.SetValue(4, 3, 1);
	lhs_ex.SetValue(4, 4, 2);
	lhs_ex.SetValue(4, 5, 1);
	lhs_ex.SetValue(5, 4, 1);
	lhs_ex.SetValue(5, 5, 2);
	// -------------- Tridiagonal part -------------- //

	// -------------- Outer Diagonal elements -------------- //
	lhs_ex.SetValue(5, 0, 1);
	lhs_ex.SetValue(0, 5, 1);
	lhs_ex.SetValue(3, 5, 1);
	lhs_ex.SetValue(5, 3, 1);
	// -------------- Outer Diagonal elements -------------- //

	std::vector<double> solution_ex(6, 1.);
	std::vector<double> test_num_solution(6);
	std::vector<double> rhs_ex = lhs_ex * solution_ex;


	// The method converges well enough for non-symmetric matrices as well
	// -------------- Symmetry check -------------- //
	for (size_t i = 0; i < 6; i++)
	{
		for (size_t j = 0; j < 6; j++)
		{
			CHECK( lhs_ex.GetValue(i, j) == lhs_ex.GetValue(j, i));
		}			
	}
	// -------------- Symmetry check -------------- //

	// ----------------------------------------------------------------------- //
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::LDU, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10);
	IMatrixSolver* ldu_solver = IMatrixSolver::Factory(params);

	ldu_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------------------------------------------------------------- //

	delete params;
	delete ldu_solver;
}