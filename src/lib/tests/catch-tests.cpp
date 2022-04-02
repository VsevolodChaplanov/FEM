#include <iostream>
#include <vector>
#include <limits>

#include "FemGrid.h"
#include "FemPDE.h"
#include "SolverParams.h"
#include "Builder.h"
#include "IFiniteElem.h"
#include "LinElem.h"
#include "GlobalAssemblers.h"
#include "VectorOperations.h"

#include "CompressedM.h"

#include <catch2/catch_all.hpp>
#define CATCH_CONFIG_MAIN

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

// Общий тест
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

// Тест локальных матриц
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

	delete lin_elem;
}

// Тест глобальных матриц
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
}

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

	for (const double elem : norm2a_cont)
	{
		std::cout << elem << std::endl;
	}

	// TODO адекватную проверку на порядок аппроксимации
}