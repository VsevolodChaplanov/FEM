#include "PDE_lib_tests.h"
#include <vector>
#include <limits>

// Tests of full solution
//		1D linear solution
// 		Approx degree

double u_ex(const double* point)
{
	return sin(9 * (point[0]+0.2) * (point[0]+0.2));
}

double f_fun(const double* point)
{
	return (point[0]+1.2) * u_ex(point);
}

double k_fun(const double* point)
{
	return 1.0 / (18 * 18 * (point[0] + 0.2));
}

// -------------------- Full solution check ------------------
TEST_CASE( "Check 1D linear solution", "[1D solution]" )
{
	std::vector<double> num_ex_sol {
		0.35227423327509,
		0.7206632643087959,
		0.97872263335775433,
		0.76006116867736739,
		-0.093379129996299423,
		-0.90236783697542888,
		-0.46810400653055212,
		0.75694560454116289,
		0.37528404123122855,
		-0.83598453120115468,
		0.38354275541260835
	};

	FemGrid femgridlinear = Builder::BuildLinear1DGrid(0, 1, 10);

	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::Thomas, MatrixSolverParams::Preconditioners::None, 1000, 1.e-5, 10);

	FemPDE fempde(&femgridlinear, f_fun, k_fun, params); 

	fempde.assemble();
	fempde.apply_boundary_condition_dirichlet(u_ex, femgridlinear.boundary_element_indices(1));
	fempde.apply_boundary_condition_dirichlet(u_ex, femgridlinear.boundary_element_indices(2));

	std::vector<double> sol = fempde.solve(); // вектор решения ин в самом методе

	// ------------------- Size check ------------------- //
	CHECK( sol.size() == num_ex_sol.size() );
	// ------------------- Size check ------------------- //

	// ------------------- Component-component comparing ------------------- //
	for (size_t i = 0; i < num_ex_sol.size(); i++)
	{
		CHECK( sol[i] == num_ex_sol[i] );
	}
	// ------------------- Component-component comparing ------------------- //

}

// ------------------- Approximation degree test ---------------------
TEST_CASE( "Degree of approximation", "[ApproxCheck]" )
{
	std::vector<double> norm2a_cont;

	for (const size_t N : {10, 100, 1000, 10000})
	{
		FemGrid femgridlinear = Builder::BuildLinear1DGrid(0, 1, N); // Для ГУ на границе стоит элемент порядка ниже
		MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::Thomas, MatrixSolverParams::Preconditioners::None, 1000, 1.e-5, 10);
		FemPDE fempde(&femgridlinear, f_fun, k_fun, params); 
		fempde.assemble();

		fempde.apply_boundary_condition_dirichlet(u_ex, femgridlinear.boundary_element_indices(1));
		fempde.apply_boundary_condition_dirichlet(u_ex, femgridlinear.boundary_element_indices(2));	
		std::vector<double> sol = fempde.solve();
		std::vector<double> u_ex_num = femgridlinear.approximate(u_ex);
		std::vector<double> difference = vector_diff(sol, u_ex_num);
		norm2a_cont.push_back(femgridlinear.norm2(difference));
	}

	std::array<double, 4> test_norms {0.0620393, 0.000684953, 6.85292e-06, 6.8539e-08};

	for (size_t i = 0; i < 4; i++)
	{
		CHECK( test_norms[i] == Approx(norm2a_cont[i]) );
	}
}
