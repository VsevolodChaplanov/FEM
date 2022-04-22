#include "PDE_lib_tests.h"

// Exist test cases
// 		Selecting boundaries test
// 		Linear finite elements mesh builder test
// 		Finite elements mesh parser
// 		Finite elements mesh builder (from file)

bool select_func_left(const double* phys_p)
{
	if (phys_p[0] == 0) 
	{
		return true;
	}
	return false;
}

bool select_func_right(const double* phys_p)
{
	if (phys_p[0] == 1) 
	{
		return true;
	}
	return false;
}

bool select_func_up(const double* phys_p)
{
	if (phys_p[1] == 1) 
	{
		return true;
	}
	return false;
}

bool select_func_down(const double* phys_p)
{
	if (phys_p[1] == 0) 
	{
		return true;
	}
	return false;
}

TEST_CASE( "Selectors", "[Selectors]" )
{
	const std::string filename = "./test_resources/coarse_triangle.vtk";
	FemGrid grid = Builder::BuildFromFile(filename);

	std::vector<size_t> indbounds_left = grid.boundary_element_indices(select_func_left);
	std::vector<size_t> indbounds_right = grid.boundary_element_indices(select_func_right);
	std::vector<size_t> indbounds_up = grid.boundary_element_indices(select_func_up);
	std::vector<size_t> indbounds_down = grid.boundary_element_indices(select_func_down);

	//  ---------------- Left side check ---------------- //
	for (const size_t index : {0, 7, 3})
	{
		bool result = std::find(indbounds_left.begin(), indbounds_left.end(), index) == indbounds_left.end() ? false : true;
		CHECK( result );
		INFO( "Index : " << index );
	}
	//  ---------------- Left side check ---------------- //

	//  ---------------- Right side check ---------------- //
	for (const size_t index : {2, 5, 1})
	{
		bool result = std::find(indbounds_right.begin(), indbounds_right.end(), index) == indbounds_right.end() ? false : true;
		CHECK( result );
		INFO( "Index : " << index );
	}
	//  ---------------- Right side check ---------------- //

	//  ---------------- Up side check ---------------- //
	for (const size_t index : {3, 6, 2})
	{
		bool result = std::find(indbounds_up.begin(), indbounds_up.end(), index) == indbounds_up.end() ? false : true;
		CHECK( result );
		INFO( "Index : " << index );
	}
	//  ---------------- Up side check ---------------- //

	//  ---------------- Down side check ---------------- //
	for (const size_t index : {0, 4, 1})
	{
		bool result = std::find(indbounds_down.begin(), indbounds_down.end(), index) == indbounds_down.end() ? false : true;
		CHECK( result );
		INFO( "Index : " << index );
	}
	//  ---------------- Down side check ---------------- //
}


// --------------- 1D linear Finite elements mesh builder test --------------
TEST_CASE("Builder checking", "[ClassBuilder]")
{
	FemGrid lin_test = Builder::BuildLinear1DGrid(0, 1, 4);

	LinearLineElement Lin1({0, 0.25}, {0, 1});
	LinearLineElement Lin2({0.25, 0.5}, {1, 2});
	LinearLineElement Lin3({0.5, 0.75}, {2, 3});
	LinearLineElement Lin4({0.75, 1.0}, {3, 4});

	std::vector<LinearLineElement> lin_ex {Lin1, Lin2, Lin3, Lin4};

	for (size_t i = 0; i < 4; i++)
	{
		REQUIRE( lin_test.get_element(i)->get_volume() == lin_ex[i].get_volume() );
	
		for (size_t i = 0; i < 2; i++)
		{
			REQUIRE( lin_test.get_element(i)->get_global_indices()[i] == lin_ex[i].get_global_indices()[i] );
		}			

		REQUIRE( lin_test.get_element(i)->get_element_type() == lin_ex[i].get_element_type() );

		REQUIRE( *lin_test.get_element(i)->get_center_coordinates() == *lin_ex[i].get_center_coordinates() );

		REQUIRE( lin_test.get_element(i)->get_lumped_matrix() == lin_ex[i].get_lumped_matrix() );

		REQUIRE( lin_test.get_element(i)->get_mass_matrix() == lin_ex[i].get_mass_matrix() );
	
		REQUIRE( lin_test.get_element(i)->get_stiffness_matrix() == lin_ex[i].get_stiffness_matrix() );
	}
}

TEST_CASE( "Finite elements mesh parser (.vtk)", "[VtkParser]" )
{
	VtkFEMParser test_vtk("./test_resources/Rect.1.1mesh.1.vtk");

	test_vtk.load_mesh();

	REQUIRE( test_vtk.get_vertices_number() == 98 );
	REQUIRE( test_vtk.get_elements_number() == 198 );
	REQUIRE( test_vtk.get_cell_types().size() == 198);

	std::vector<double> vertices_test = test_vtk.get_vertices();
	std::vector<std::vector<size_t>> cells_test = test_vtk.get_cells();
	std::vector<size_t> cell_types_test = test_vtk.get_cell_types();

	// -------------- Vertices data check -------------- // 
			// First vertice
	for (size_t i = 0; i < 3; i++)
	{
		CHECK( vertices_test[i] == 0.0 );
	}
	
	// Random vertice from the centre
	CHECK( vertices_test[177] ==  0.31331768157997741);
	CHECK( vertices_test[178] ==  0.75237107335505982);
	CHECK( vertices_test[179] ==  0.0);

	// Last vertice
	CHECK( vertices_test[98 * 3 - 1] == 0.0);
	CHECK( vertices_test[98 * 3 - 2] == 0.8253550026152658);
	CHECK( vertices_test[98 * 3 - 3] == 0.7489644266062752);
	// -------------- Vertices data check -------------- // 

	// -------------- Cells data check -------------- // 
	// Nodes of the rectangle
	REQUIRE( cells_test[0].size() == 1);
	CHECK( cells_test[0][0] == 0 );
	// Edges of the rectangle
	REQUIRE( cells_test[4].size() == 2);
	CHECK( cells_test[4][0] == 0 );
	CHECK( cells_test[4][1] == 4 );
	// Internal triangles
	REQUIRE( cells_test[37].size() == 3);
	CHECK( cells_test[36][0] == 67);
	CHECK( cells_test[36][1] == 78);
	CHECK( cells_test[36][2] == 37);
	// Last triangle
	REQUIRE( cells_test[197].size() == 3);
	CHECK( cells_test[97][0] == 6);
	CHECK( cells_test[97][1] == 36);
	CHECK( cells_test[97][2] == 5);
	// Last triangle
	REQUIRE( cells_test[197].size() == 3);
	CHECK( cells_test[197][0] == 82);
	CHECK( cells_test[197][1] == 97);
	CHECK( cells_test[197][2] == 60);
	// -------------- Cells data check -------------- // 

	// -------------- Cell types data check -------------- // 

	// All cell types check
	for (std::vector<size_t>::iterator it = cell_types_test.begin(); it != cell_types_test.begin() + 4; it++)
	{
		CHECK( *it == 1 );
	}
	for (std::vector<size_t>::iterator it = cell_types_test.begin() + 4; it != cell_types_test.begin() + 36; it++)
	{
		CHECK( *it == 3 );
	}
	for (std::vector<size_t>::iterator it = cell_types_test.begin() + 36; it != cell_types_test.end(); it++)
	{
		CHECK( *it == 5 );
	}
	// -------------- Cell types data check -------------- // 

}

TEST_CASE( "Finite lements mesh builder from file .vtk", "[BuildAlgoTest]" )
{
	const std::string filename = "./test_resources/coarse_triangle.vtk";
	//FemGrid triang_test = Builder::BuildFromFile(filename);

	size_t p_dim = 2;

	FemGrid grid_test = Builder::BuildFromFile(filename);

	REQUIRE( grid_test.get_elements_number() == 14 );
	REQUIRE( grid_test.get_vertices_number() == 12 );
	REQUIRE( grid_test.get_all_elements_number() == 26);

	// --------------- Point bound check --------------- //
	const IBoundaryElement* belem_test_p = grid_test.get_boundary_element(0);
	CHECK(
		compare_vectors<size_t>(belem_test_p->get_global_indices(), {0})
	);
	CHECK(
		belem_test_p->get_element_type() == 0
	);
	delete belem_test_p;
	// --------------- Point bound check --------------- //

	// --------------- Line bound check --------------- //
	const IBoundaryElement* belem_test_l = grid_test.get_boundary_element(4);
	CHECK(
		compare_vectors<size_t>(belem_test_l->get_global_indices(), {0,4})
	);
	CHECK(
		belem_test_l->get_element_type() == 1
	);
	delete belem_test_l;
	// --------------- Line bound check --------------- //

	// --------------- Triangle element check --------------- //
	const IFiniteElement* elem_test = grid_test.get_element(0);
	CHECK(
		compare_vectors(elem_test->get_global_indices(), {0,9,7})
	);
	CHECK(
		elem_test->get_element_type() == 2
	);
	delete elem_test;
	// --------------- Triangle element check --------------- //
}

