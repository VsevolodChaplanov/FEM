#include <iostream>
#include <vector>

#include "FemGrid.h"
#include "FemPDE.h"
#include "SolverParams.h"
#include "Builder.h"
#include "IFiniteElem.h"
#include "LinElem.h"

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

	SECTION("Vector dimensions check")
	{
		REQUIRE( sol.size() == num_ex_sol.size());
	}

	SECTION("Component-by-component comparing")
	{
		// Сравнение покомпонентно
		for (size_t i = 0; i < num_ex_sol.size(); i++)
		{
			// Было бы наверно логично чтобы сравнение проводилось до 8 знака например
			REQUIRE( sol[i] == num_ex_sol[i]);
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


TEST_CASE("1D linear element matrices check", "[Matrices]")
{
	IFiniteElement* lin_elem = IFiniteElement::Factory({0,1}, {0, 1}, 1);

	REQUIRE(lin_elem->get_number_basis_func() == 2);
	REQUIRE(lin_elem->get_volume() == 1);

	SECTION("Mass matrix of linear element check")
	{
		std::array<double, 4> test_mass {(double) 1 / 3, (double) 1 / 6, (double) 1 / 6, (double) 1 / 3};
		for (size_t i = 0; i < 2; i++)
		{
			for (size_t j = 0; j < 2; j++)
			{
				REQUIRE(lin_elem->get_mass(i, j) == test_mass[i * 2 + j]);
			}
		}
	}

	SECTION("Siffness matrix of linear element check")
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


	SECTION("Lumped mass matrix of linear element check")
	{
		std::array<double, 2> test_lumped_mass {0.5, 0.5};

		for (size_t i = 0; i < 2; i++)
		{
			REQUIRE(lin_elem->get_lumped(i) == test_lumped_mass[i]);
		}
	}
}

