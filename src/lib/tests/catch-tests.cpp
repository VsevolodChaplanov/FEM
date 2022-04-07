#define CATCH_CONFIG_RUNNER
#include "catch-tests.h"
#include "FemGrid.h"
#include "FemPDE.h"
#include "SolverParams.h"
#include "Builder.h"
#include "IFiniteElem.h"
#include "LinElem.h"
#include "GlobalAssemblers.h"
#include "VectorOperations.h"
#include "CompressedM.h"
#include <vector>
#include <limits>


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
TEST_CASE("Chech 1D linear solution", 
	"[1D solution]")
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

	// Builder mathes to left boundary 1 to right boundary 2
	FemGrid femgridlinear = Builder::BuildLinear1DGrid(0, 1, 10); // Для ГУ на границе стоит элемент порядка ниже
	// femgtidlinear->assign_boundary_type(selector , type)

	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::Thomas, MatrixSolverParams::Preconditioners::None, 1000, 1.e-5, 10);

	FemPDE fempde(&femgridlinear, f_fun, k_fun, params); 

	fempde.assemble();
	fempde.apply_boundary_condition_dirichlet(u_ex, femgridlinear.boundary_element_indices(1));
	fempde.apply_boundary_condition_dirichlet(u_ex, femgridlinear.boundary_element_indices(2));
	// 	fempde->apply_boundary_condition_dirichlet(u_exact, femgridlinear->boundary_element_indices()); All

	std::vector<double> sol = fempde.solve(); // вектор решения ин в самом методе

	SECTION( "Vector dimensions check" )
	{
		REQUIRE( sol.size() == num_ex_sol.size() );
	}

	SECTION( "Component-by-component comparing" )
	{
		// Сравнение покомпонентно
		for (size_t i = 0; i < num_ex_sol.size(); i++)
		{
			// Было бы наверно логично чтобы сравнение проводилось до 8 знака например
			REQUIRE( sol[i] == num_ex_sol[i] );
		}
	}

	// Approximate analytical function along calculation area
	// std::vector<double> u_ex_num = femgridlinear.approximate(u_ex);

	// SECTION("Component-by-component comparing with exact solution")
	// {// Сравнение покомпонентно
	// 	for (size_t i = 0; i < num_ex_sol.size(); i++)
	// 	{
	// 		// Было бы наверно логично чтобы сравнение проводилось до 8 знака например,
	// 		// Или вообще не стоит делать такое сравнение т.к. в 1 месте тест падает
	// 		REQUIRE( sol[i] == u_ex_num[i]);
	// 	}
	// }
}

// ----------------------- 1D linear element matrices check -----------------------
TEST_CASE( "1D linear element matrices check", "[LocalMatrices]" )
{
	IFiniteElement* lin_elem = IFiniteElement::Factory({0,1}, {0, 1}, 1);

	REQUIRE( lin_elem->get_number_basis_func() == 2 );
	REQUIRE( lin_elem->get_volume() == 1 );

	SECTION( "Mass matrix of linear element check" )
	{
		std::array<double, 4> test_mass {(double) 1 / 3, (double) 1 / 6, (double) 1 / 6, (double) 1 / 3};
		for (size_t i = 0; i < 2; i++)
		{
			for (size_t j = 0; j < 2; j++)
			{
				REQUIRE( lin_elem->get_mass(i, j) == test_mass[i * 2 + j] );
			}
		}
	}

	SECTION( "Siffness matrix of linear element check" )
	{
		std::array<double, 4> test_stiffness {1, - 1, - 1, 1};
		for (size_t i = 0; i < 2; i++)
		{
			for (size_t j = 0; j < 2; j++)
			{
				REQUIRE(lin_elem->get_stiffness(i, j) == test_stiffness[i * 2 + j]);
			}
		}
	}


	SECTION( "Lumped mass matrix of linear element check" )
	{
		std::array<double, 2> test_lumped_mass {0.5, 0.5};

		for (size_t i = 0; i < 2; i++)
		{
			REQUIRE( lin_elem->get_lumped(i) == test_lumped_mass[i] );
		}
	}

	SECTION( "Volume Check" )
	{
		REQUIRE( lin_elem->get_volume() == 1. );
	}

	SECTION( "Global indices check" )
	{
		const std::vector<double>& test_g_ind {0, 1};
		for (size_t i = 0; i < 2; i++)
		{
			REQUIRE( lin_elem->get_global_indices()[i] == test_g_ind[i] );
		}
		
	}

	SECTION( "Element types" )
	{
		REQUIRE( lin_elem->get_element_type() == 1 );
	}

	SECTION( "Centre coordinates" )
	{
		REQUIRE( *(lin_elem->get_center_coordinates()) == 0.5 );
	}

	delete lin_elem;
}

// -------------------- Point boundary element
TEST_CASE( "Point boundary element", "[PountElementCheck]" )
{
	IBoundaryElement* point_left_test = IBoundaryElement::Factory({0}, {0}, 0, 1);
	IBoundaryElement* point_right_test = IBoundaryElement::Factory({1}, {1}, 0, 2);



	SECTION( "Bound types" )
	{
		CHECK( point_left_test->get_bound_type() == 1 );
		CHECK( point_right_test->get_bound_type() == 2 );
	}

	SECTION( "Element types" )
	{
		CHECK( point_left_test->get_element_type() == 0 );
		CHECK( point_right_test->get_element_type() == 0 );
	}

	SECTION( "Centre coordinates" )
	{
		CHECK( *point_left_test->get_center_coordinates() == 0 );
		CHECK( *point_right_test->get_center_coordinates() == 1 );
	}

	SECTION( "Global indices" )
	{
		REQUIRE( point_left_test->get_global_indices().size() == 1 ); 
		REQUIRE( point_right_test->get_global_indices().size() == 1 );
		CHECK( point_left_test->get_global_indices()[0] == 0 ) ;
		CHECK( point_right_test->get_global_indices()[0] == 1 );
	}

	SECTION( "Volume" )
	{
		CHECK( point_left_test->get_volume() == 0);
		CHECK( point_right_test->get_volume() == 0);
	}

	SECTION( "Local Matrices" )
	{
		std::vector<double> one {1};
		std::vector<double> zer {0};
		CHECK( point_left_test->get_lumped_matrix() == one );
		CHECK( point_right_test->get_lumped_matrix() == one );
		CHECK( point_left_test->get_mass_matrix() == one );
		CHECK( point_right_test->get_mass_matrix() == one );
		CHECK( point_left_test->get_stiffness_matrix() == zer );
		CHECK( point_right_test->get_stiffness_matrix() == zer );
	}
}

// ---------------------- Global matrices check ----------------------
TEST_CASE( "Global matrices check", "[GlobalMatrices]" )
{
	FemGrid lin_ex = Builder::BuildLinear1DGrid(0, 1, 4);

	GLobalMatrixAssembler mass_g_matrix(5);
	GLobalMatrixAssembler stiffness_g_matrix(5);

	// Assemble test mass matrix
	for (size_t i = 0; i < 4; i++)
	{
		mass_g_matrix.add_local_matrix(lin_ex.get_element(i)->get_global_indices(), lin_ex.get_element(i)->get_mass_matrix());
	}

	// Assemble test stiffness matrix
	// k = 1
	for (size_t i = 0; i < 4; i++)
	{
		stiffness_g_matrix.add_local_matrix(lin_ex.get_element(i)->get_global_indices(), lin_ex.get_element(i)->get_stiffness_matrix(), 1);
	}

	CMatrix mass_ex(5);
	CMatrix stiff_ex(5);

	mass_ex.SetValue(0, 0, 0.083333333333333329);
	mass_ex.SetValue(0, 1, 0.041666666666666664);
	mass_ex.SetValue(1, 0, 0.041666666666666664);
	mass_ex.SetValue(1, 1, 0.16666666666666666);
	mass_ex.SetValue(1, 2, 0.041666666666666664);
	mass_ex.SetValue(2, 1, 0.041666666666666664);
	mass_ex.SetValue(2, 2, 0.1666666666666666);
	mass_ex.SetValue(2, 3, 0.041666666666666664);
	mass_ex.SetValue(3, 2, 0.041666666666666664);
	mass_ex.SetValue(3, 3, 0.16666666666666666);
	mass_ex.SetValue(3, 4, 0.041666666666666664);
	mass_ex.SetValue(4, 3, 0.041666666666666664);
	mass_ex.SetValue(4, 4, 0.083333333333333329);

	stiff_ex.SetValue(0, 0, 4);
	stiff_ex.SetValue(0, 1, -4);
	stiff_ex.SetValue(1, 0, -4);
	stiff_ex.SetValue(1, 1, 8);
	stiff_ex.SetValue(1, 2, -4);
	stiff_ex.SetValue(2, 1, -4);
	stiff_ex.SetValue(2, 2, 8);
	stiff_ex.SetValue(2, 3, -4);
	stiff_ex.SetValue(3, 2, -4);
	stiff_ex.SetValue(3, 3, 8);
	stiff_ex.SetValue(3, 4, -4);
	stiff_ex.SetValue(4, 3, -4);
	stiff_ex.SetValue(4, 4, 4);

	CMatrix mass_test = mass_g_matrix.get_result();
	CMatrix stiff_test = stiffness_g_matrix.get_result();


	SECTION( "Global matrices check" )
	{
		for (size_t i = 0; i < 5; i++)
		{
			for (size_t j = 0; j < 5; j++)
			{
				INFO( "i = " << i);
				INFO( "j = " << j);
				
				// Тест падает на величинах 1.666666667 == 1.666666667 
				REQUIRE( mass_test.GetValue(i, j) == Catch::Approx(mass_ex.GetValue(i, j)) );
				REQUIRE( stiff_test.GetValue(i, j) == Catch::Approx(stiff_ex.GetValue(i, j)) );
			}
		}
	}

	SECTION("Right hand side vector check", "[RhsCheck]")
	{
		std::vector<double> rhs_ex {0.125,
			0.24999999999999997,
			0.24999999999999997,
			0.24999999999999997,
			0.125
		};

		std::vector<double> I(5, 1.);
		std::vector<double> rhs_test = mass_test * I;

		for (size_t i = 0; i < 5; i++)
		{	
			INFO( "i = " << i);
			REQUIRE( rhs_test[i] == rhs_ex[i] );
		}
	}
}

// ------------------- Approximation degree test ---------------------
TEST_CASE( "Degree of approximation", "[ApproxCheck]" )
{
	std::vector<double> norm2a_cont;

	for (size_t N = 10; N < 10001; N *= 10)
	{
		FemGrid femgridlinear = Builder::BuildLinear1DGrid(0, 1, N); // Для ГУ на границе стоит элемент порядка ниже
		MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::Thomas, MatrixSolverParams::Preconditioners::None, 1000, 1.e-5, 10);
		FemPDE fempde(&femgridlinear, f_fun, k_fun, params); 
		fempde.new_assembler();

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
		REQUIRE( test_norms[i] == Catch::Approx(norm2a_cont[i]) );
	}
	
}

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

	SECTION( "Tridiagonal check" )
	{
		for (size_t i = 0; i < 6; i++)
		{
			for (auto elem : lhs_ex[i])
			{
				INFO( "i = " << i);
				bool res = elem.first == i || elem.first == i - 1 || elem.first == i + 1;
				REQUIRE( res == true);
			}
		}
	}

	std::vector<double> solution_ex(6, 1.);
	std::vector<double> rhs_ex = lhs_ex * solution_ex;

	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::Thomas, MatrixSolverParams::Preconditioners::None, 100, 1.0e-4, 10);
	IMatrixSolver* thomas_solver = IMatrixSolver::Factory(params);

	std::vector<double> test_num_solution(6);
	thomas_solver->solve(lhs_ex, rhs_ex, test_num_solution);

	SECTION("Solution check")
	{
		for (size_t i = 0; i < 6; i++)
		{
			INFO( "i = " << i);
			REQUIRE( Catch::Approx(test_num_solution[i]).epsilon(1.0e-4) == solution_ex[i] );
		}
	}
	

	delete params;
	delete thomas_solver;
}

// --------------- 1D linear Finite elements mesh builder test --------------
TEST_CASE("Builder checking", "[ClassBuilder]")
{
	FemGrid lin_test = Builder::BuildLinear1DGrid(0, 1, 4);

	LinElem Lin1({0, 0.25}, {0, 1});
	LinElem Lin2({0.25, 0.5}, {1, 2});
	LinElem Lin3({0.5, 0.75}, {2, 3});
	LinElem Lin4({0.75, 1.0}, {3, 4});

	std::vector<LinElem> lin_ex {Lin1, Lin2, Lin3, Lin4};

	for (size_t i = 0; i < 4; i++)
	{
		SECTION( "Volume Check" )
		{
			REQUIRE( lin_test.get_element(i)->get_volume() == lin_ex[i].get_volume() );
		}

		SECTION( "Global indices check" )
		{
			for (size_t i = 0; i < 2; i++)
			{
				REQUIRE( lin_test.get_element(i)->get_global_indices()[i] == lin_ex[i].get_global_indices()[i] );
			}			
		}

		SECTION( "Element types" )
		{
			REQUIRE( lin_test.get_element(i)->get_element_type() == lin_ex[i].get_element_type() );
		}

		SECTION( " Center coordinates" )
		{
			REQUIRE( *lin_test.get_element(i)->get_center_coordinates() == *lin_ex[i].get_center_coordinates() );
		}

		SECTION( " Lumped matrix " )
		{
			REQUIRE( lin_test.get_element(i)->get_lumped_matrix() == lin_ex[i].get_lumped_matrix() );
		}

		SECTION( " Mass matrix " )
		{
			REQUIRE( lin_test.get_element(i)->get_mass_matrix() == lin_ex[i].get_mass_matrix() );
		}

		SECTION( " Stiffness matrix " )
		{
			REQUIRE( lin_test.get_element(i)->get_stiffness_matrix() == lin_ex[i].get_stiffness_matrix() );
		}
	}
}

TEST_CASE( "Math operations", "[Math]" )
{
	std::vector<double> first_ex {1.5 , 0.5, 1.5, -1.5};
	std::vector<double> second_ex {0.5, -1.5, 0.5, -0.5};

	CMatrix mat_ex(4);

	mat_ex.SetValue(0, 0, 1);
	mat_ex.SetValue(0, 1, -2);
	mat_ex.SetValue(0, 2, 3);
	mat_ex.SetValue(0, 3, -1);
	mat_ex.SetValue(1, 0, 3);
	mat_ex.SetValue(1, 1, 6);	
	mat_ex.SetValue(1, 2, -3);
	mat_ex.SetValue(1, 3, -1);
	mat_ex.SetValue(2, 0, -1);
	mat_ex.SetValue(2, 1, 2);
	mat_ex.SetValue(2, 2, 4);
	mat_ex.SetValue(2, 3, -4);
	mat_ex.SetValue(3, 0, 2);
	mat_ex.SetValue(3, 1, -3);
	mat_ex.SetValue(3, 2, 4);
	mat_ex.SetValue(3, 3, 5);

	std::vector<double> test_result = mat_ex * first_ex;

	SECTION( "Matrix-vector multyply" )
	{
		CHECK( test_result[0] == 6.5 );
		CHECK( test_result[1] == 4.5);
		CHECK( test_result[2] == 11.5 );
		CHECK( test_result[3] == 0);		
	}

	SECTION( "Vector summ" )
	{
		std::vector<double> sum = vector_sum(first_ex, second_ex);
		CHECK( sum[0] == 2 );
		CHECK( sum[1] == -1 );
		CHECK( sum[2] == 2);
		CHECK( sum[3] == -2);
	}

	SECTION( "Scalar mult" )
	{
		double sc_test = dot_product(first_ex, second_ex);
		CHECK( sc_test == 1.5 );
	}

	CMatrix mat_ex_2(4);

	mat_ex_2.SetValue(0, 0, 3);
	mat_ex_2.SetValue(0, 1, 2);
	mat_ex_2.SetValue(0, 2, -1);
	mat_ex_2.SetValue(0, 3, 6);
	mat_ex_2.SetValue(1, 0, -2);
	mat_ex_2.SetValue(1, 1, 1);	
	mat_ex_2.SetValue(1, 2, -4);
	mat_ex_2.SetValue(1, 3, 7);
	mat_ex_2.SetValue(2, 0, 5);
	mat_ex_2.SetValue(2, 1, -4);
	mat_ex_2.SetValue(2, 2, 2);
	mat_ex_2.SetValue(2, 3, 4);
	mat_ex_2.SetValue(3, 0, 3);
	mat_ex_2.SetValue(3, 1, 4);
	mat_ex_2.SetValue(3, 2, 6);
	mat_ex_2.SetValue(3, 3, -2);

	SECTION( "Sum of sparse matrices" )
	{
		CMatrix res_test(4);
		summ_cm(mat_ex, mat_ex_2, res_test);
		REQUIRE( res_test.size() == 4);

		CHECK( res_test.GetValue(0, 0) == 4 );
		CHECK( res_test.GetValue(0, 1) == 0 );
		CHECK( res_test.GetValue(0, 2) == 2 );
		CHECK( res_test.GetValue(0, 3) == 5 );
		CHECK( res_test.GetValue(1, 0) == 1 );
		CHECK( res_test.GetValue(1, 1) == 7 );
		CHECK( res_test.GetValue(1, 2) == -7 );
		CHECK( res_test.GetValue(1, 3) == 6 );
		CHECK( res_test.GetValue(2, 0) == 4 );
		CHECK( res_test.GetValue(2, 1) == -2 );
		CHECK( res_test.GetValue(2, 2) == 6 );
		CHECK( res_test.GetValue(2, 3) == 0 );
		CHECK( res_test.GetValue(3, 0) == 5 );
		CHECK( res_test.GetValue(3, 1) == 1 );
		CHECK( res_test.GetValue(3, 2) == 10 );
		CHECK( res_test.GetValue(3, 3) == 3 );
	}

	SECTION( "Sparse matrix multyply" )
	{
		CMatrix res_test = mat_ex * mat_ex_2;
		REQUIRE( res_test.size() == 4);

		CHECK( res_test.GetValue(0, 0) == 19 );
		CHECK( res_test.GetValue(0, 1) == -16 );
		CHECK( res_test.GetValue(0, 2) == 7 );
		CHECK( res_test.GetValue(0, 3) == 6 );
		CHECK( res_test.GetValue(1, 0) == -21 );
		CHECK( res_test.GetValue(1, 1) == 20 );
		CHECK( res_test.GetValue(1, 2) == -39 );
		CHECK( res_test.GetValue(1, 3) == 50 );
		CHECK( res_test.GetValue(2, 0) == 1 );
		CHECK( res_test.GetValue(2, 1) == -32 );
		CHECK( res_test.GetValue(2, 2) == -23 );
		CHECK( res_test.GetValue(2, 3) == 32 );
		CHECK( res_test.GetValue(3, 0) == 47 );
		CHECK( res_test.GetValue(3, 1) == 5 );
		CHECK( res_test.GetValue(3, 2) == 48 );
		CHECK( res_test.GetValue(3, 3) == -3 );
	}
}

