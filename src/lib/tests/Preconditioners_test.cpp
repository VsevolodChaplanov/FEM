#include "PDE_lib_tests.h"
// Preconditioned methods
//		Conjugate gradients
//			Jacobi
//			SSOR
//			ILU0 (first impl)
// 			ILU0 (second impl)
//		Gradient Descent
//			Jacobi
//			SSOR
//			ILU0 (first impl)
// 			ILU0 (second impl)

TEST_CASE( "Conjugate Gradients SLAE solver check with Jacobi Precondition", "[CGSolverJacobi][ExactMethods][Preconditions]" )
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
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::CG, MatrixSolverParams::Preconditioners::Jacobi_P, 100, 1.0e-4, 10);
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

TEST_CASE( "Conjugate Gradients SLAE solver check with SSOR Precondition", "[CGSolverSSOR][ExactMethods][Preconditions]" )
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
	// ----------------------- Relaxation param 1.95 ------------------------- //
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::CG, MatrixSolverParams::Preconditioners::SSOR_P, 100, 1.0e-4, 10, 1.95, 1.95);
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
	// ----------------------- Relaxation param 1.95 ------------------------- //
	// ----------------------------------------------------------------------- //


	// ----------------------------------------------------------------------- //
	// ----------------------- Relaxation param 1.975 ------------------------ //
	params = new MatrixSolverParams(MatrixSolverParams::Methods::CG, MatrixSolverParams::Preconditioners::SSOR_P, 100, 1.0e-4, 10, 1.95, 1.975);
	cg_solver = IMatrixSolver::Factory(params);

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
	// ----------------------- Relaxation param 1.975 ------------------------ //
	// ----------------------------------------------------------------------- //


	// ----------------------------------------------------------------------- //
	// ----------------------- Relaxation param 1.99 ------------------------- //
	params = new MatrixSolverParams(MatrixSolverParams::Methods::CG, MatrixSolverParams::Preconditioners::SSOR_P, 100, 1.0e-4, 10, 1.95, 1.99);
	cg_solver = IMatrixSolver::Factory(params);

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
	// ----------------------- Relaxation param 1.99 ------------------------- //
	// ----------------------------------------------------------------------- //

	delete params;
	delete cg_solver;
}

TEST_CASE( "Conjugate Gradients SLAE solver check with ILU0 Precondition 1st", "[CGSolverILU01][ExactMethods][Preconditions]" )
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
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::CG, MatrixSolverParams::Preconditioners::ILU01_P, 100, 1.0e-4, 10);
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

TEST_CASE( "Conjugate Gradients SLAE solver check with ILU0 Precondition 2st", "[CGSolverILU02][ExactMethods][Preconditions]" )
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
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::CG, MatrixSolverParams::Preconditioners::ILU02_P, 100, 1.0e-4, 10);
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


TEST_CASE( "Gradient Descent SLAE solver check with Jacobi Precondition", "[CGSolverJacobi][Preconditions]" )
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
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::GD, MatrixSolverParams::Preconditioners::Jacobi_P, 100, 1.0e-4, 10);
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

TEST_CASE( "Gradient Descent SLAE solver check with SSOR Precondition", "[CGSolverSSOR][Preconditions]" )
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
	// ----------------------- Relaxation param 1.95 ------------------------- //
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::GD, MatrixSolverParams::Preconditioners::SSOR_P, 100, 1.0e-4, 10, 1.95, 1.95);
	IMatrixSolver* gd_solver = IMatrixSolver::Factory(params);

	gd_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------------- Relaxation param 1.95 ------------------------- //
	// ----------------------------------------------------------------------- //


	// ----------------------------------------------------------------------- //
	// ----------------------- Relaxation param 1.975 ------------------------ //
	params = new MatrixSolverParams(MatrixSolverParams::Methods::GD, MatrixSolverParams::Preconditioners::SSOR_P, 100, 1.0e-4, 10, 1.95, 1.975);
	gd_solver = IMatrixSolver::Factory(params);

	gd_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------------- Relaxation param 1.975 ------------------------ //
	// ----------------------------------------------------------------------- //


	// ----------------------------------------------------------------------- //
	// ----------------------- Relaxation param 1.99 ------------------------- //
	params = new MatrixSolverParams(MatrixSolverParams::Methods::GD, MatrixSolverParams::Preconditioners::SSOR_P, 100, 1.0e-4, 10, 1.95, 1.99);
	gd_solver = IMatrixSolver::Factory(params);

	gd_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	// ----------------- Numerical solution check ----------------- //
	for (size_t i = 0; i < 6; i++)
	{
		INFO( "i = " << i);
		CHECK( Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
	}
	// ----------------- Numerical solution check ----------------- //
	// ----------------------- Relaxation param 1.99 ------------------------- //
	// ----------------------------------------------------------------------- //

	delete params;
	delete gd_solver;
}

TEST_CASE( "Gradient Descent SLAE solver check with ILU0 Precondition 1st", "[CGSolverILU01][Preconditions]" )
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
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::GD, MatrixSolverParams::Preconditioners::ILU01_P, 100, 1.0e-4, 10);
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

TEST_CASE( "Gradient Descent SLAE solver check with ILU0 Precondition 2st", "[CGSolverILU02][Preconditions]" )
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
	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::GD, MatrixSolverParams::Preconditioners::ILU02_P, 100, 1.0e-4, 10);
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