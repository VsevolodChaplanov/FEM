#include "PDE_lib_tests.h"

// Exist test cases
// Finite elements propertices of 
// 		Linear line element 1D mesh
// 		Linear line element 2D mesh
// 		Linear line element 3D mesh
//		Linear triangle element 2D mesh
//		Linear triangle element 3D mesh
// 		Linear point boundary element 1D mesh
// 		Linear point boundary element 2D mesh
// 		Linear point boundary element 3D mesh
// 		Linear line boundary element 2D mesh
// 		Linear line boundary element 3D mesh 

TEST_CASE( "Linear line finite element 1D mesh", "[LineFE1D]" )
{
	LinearLineElement lin_elem_test({0, 1}, {0, 1});

	CHECK( lin_elem_test.get_number_basis_func() == 2 );
	CHECK( lin_elem_test.get_volume() == 1 );
	CHECK( compare_vectors(lin_elem_test.get_global_indices(), {0, 1}) );
	CHECK( lin_elem_test.get_element_type() == 1 );
	CHECK( *lin_elem_test.get_center_coordinates() == 0.5);	
}


TEST_CASE( "Linear line finite element 2D mesh", "[LineFE2D]" )
{
	LinearLineElement lin_elem_test({0,0, 1,1}, {0, 1});

	CHECK( lin_elem_test.get_number_basis_func() == 2 );
	CHECK( lin_elem_test.get_volume() == sqrt(2) );
	CHECK( compare_vectors(lin_elem_test.get_global_indices(), {0, 1}) );
	CHECK( lin_elem_test.get_element_type() == 1 );

	CHECK( lin_elem_test.get_center_coordinates()[0] == 0.5);
	CHECK( lin_elem_test.get_center_coordinates()[1] == 0.5);	
}


TEST_CASE( "Linear line finite element 3D mesh", "[LineFE3D]" )
{
	LinearLineElement lin_elem_test({0,0,0, 1,1,1}, {0, 1});

	CHECK( lin_elem_test.get_number_basis_func() == 2 );
	CHECK( lin_elem_test.get_volume() == sqrt(3) );
	CHECK( compare_vectors(lin_elem_test.get_global_indices(), {0, 1}) );
	CHECK( lin_elem_test.get_element_type() == 1 );
	CHECK( lin_elem_test.get_center_coordinates()[0] == 0.5);
	CHECK( lin_elem_test.get_center_coordinates()[1] == 0.5);
	CHECK( lin_elem_test.get_center_coordinates()[2] == 0.5);	
}

TEST_CASE( "Lenear triangle finite element 2D mesh", "[TriangleFE2D]" )
{
	LinearTriangleElement triang_elem_test({0,0, 1,0, 0,1}, {0, 1, 2});

	CHECK( triang_elem_test.get_center_coordinates()[0] == 1.0 / 3.0 );
	CHECK( triang_elem_test.get_center_coordinates()[1] == 1.0 / 3.0 );

	CHECK( triang_elem_test.get_volume() == 0.5);
	CHECK( triang_elem_test.get_element_type() == 2);
	CHECK( 
		compare_vectors(triang_elem_test.get_global_indices(), {0, 1, 2})
	);
}

TEST_CASE( "Lenear triangle finite element 3D mesh", "[TriangleFE3D]" )
{
	LinearTriangleElement triang_elem_test({0,0,0, 1,0,0, 0,1,0}, {0, 1, 2});

	CHECK( triang_elem_test.get_center_coordinates()[0] == 1.0 / 3.0 );
	CHECK( triang_elem_test.get_center_coordinates()[1] == 1.0 / 3.0 );
	CHECK( triang_elem_test.get_center_coordinates()[2] == 0.0 );

	CHECK( triang_elem_test.get_volume() == 0.5);
	CHECK( triang_elem_test.get_element_type() == 2);
	CHECK( 
		compare_vectors(triang_elem_test.get_global_indices(), {0, 1, 2})
	);
}

TEST_CASE( "Linear point boundary finite element 1D mesh", "[PointBFE1D]" )
{
	LinearPointBoundaryElement point_belem_test({0}, {0});

	CHECK( point_belem_test.get_element_type() == 0);
	CHECK( 
		compare_vectors(point_belem_test.get_global_indices(), {0})
	);
}

TEST_CASE( "Linear point boundary finite element 2D mesh", "[PointBFE2D]" )
{
	LinearPointBoundaryElement point_belem_test({0, 1}, {0});

	CHECK( point_belem_test.get_element_type() == 0);
	CHECK( 
		compare_vectors(point_belem_test.get_global_indices(), {0})
	);
}

TEST_CASE( "Linear point boundary finite element 3D mesh", "[PointBFE3D]" )
{
	LinearPointBoundaryElement point_belem_test({0, 1, 0}, {0});

	CHECK( point_belem_test.get_element_type() == 0);
	CHECK( 
		compare_vectors(point_belem_test.get_global_indices(), {0})
	);
}

TEST_CASE( "Linear line boundary finite element 3D mesh", "[LineBFE3D]" )
{
	LinearLineBoundaryElement line_belem_test({0,0,0, 1,1,1}, {0, 1});

	CHECK( line_belem_test.get_center_coordinates()[0] == 0.5 );
	CHECK( line_belem_test.get_center_coordinates()[1] == 0.5 );
	CHECK( line_belem_test.get_center_coordinates()[2] == 0.5 );

	CHECK( line_belem_test.get_element_type() == 1 );
	CHECK( compare_vectors(line_belem_test.get_global_indices(), {0, 1}) );
	CHECK( line_belem_test.get_volume() == sqrt(3) );
}

TEST_CASE( "Linear line boundary finite element 2D mesh", "[LineBFE2D]" )
{
	LinearLineBoundaryElement line_belem_test({0,0, 1,1}, {0, 1});

	CHECK( line_belem_test.get_center_coordinates()[0] == 0.5 );
	CHECK( line_belem_test.get_center_coordinates()[1] == 0.5 );

	CHECK( line_belem_test.get_element_type() == 1 );
	CHECK( compare_vectors(line_belem_test.get_global_indices(), {0, 1}) );
	CHECK( line_belem_test.get_volume() == sqrt(2) );
}