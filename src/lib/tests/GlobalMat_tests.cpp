#include "PDE_lib_tests.h"

// Exist test cases
// 		Global matrices assemble 

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

	// ---------- Global matrices check ---------- // 
	for (size_t i = 0; i < 5; i++)
	{
		for (size_t j = 0; j < 5; j++)
		{
			INFO( "i = " << i);
			INFO( "j = " << j);
			 
			REQUIRE( mass_test.GetValue(i, j) == Approx(mass_ex.GetValue(i, j)) );
			REQUIRE( stiff_test.GetValue(i, j) == Approx(stiff_ex.GetValue(i, j)) );
		}
	}
	// ------------------------------------------- //

	// ---------- Right hand side vector check ---------- //
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
	// -------------------------------------------------- //
}