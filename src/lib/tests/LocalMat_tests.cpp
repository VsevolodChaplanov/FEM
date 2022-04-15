#include "PDE_lib_tests.h"

// Exist test cases
// Local matrices of 
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

TEST_CASE( "Local matrices of linear line element 1D mesh", "[LocMatricesLine1D]" )
{
	LinearLineElement lin_elem_test({0, 1}, {0, 1});

	std::vector<double> test_mass {(double) 1 / 3, (double) 1 / 6, (double) 1 / 6, (double) 1 / 3};
	std::vector<double> test_stiffness {1, - 1, - 1, 1};
	std::vector<double> test_lumped_mass {0.5, 0.5};

	CHECK( compare_vectors(lin_elem_test.get_mass_matrix(), test_mass) );
	CHECK( compare_vectors(lin_elem_test.get_stiffness_matrix(), test_stiffness) );
	CHECK( compare_vectors(lin_elem_test.get_lumped_matrix(), test_lumped_mass) );
}

TEST_CASE( "Local matrices of linear line element 2D mesh", "[LocMatricesLine2D]" )
{
	LinearLineElement lin_elem_test({0,0, 1,1}, {0, 1});
	double length_test_elem = sqrt(2);

	std::vector<double> test_mass {(double) length_test_elem / 3, (double) length_test_elem / 6, (double) length_test_elem / 6, (double) length_test_elem / 3};
	std::vector<double> test_stiffness {1.0 / length_test_elem, - 1.0 / length_test_elem, - 1.0 / length_test_elem, 1.0 / length_test_elem};
	std::vector<double> test_lumped_mass {length_test_elem / 2.0, length_test_elem / 2.0};

	CHECK( compare_vectors(lin_elem_test.get_mass_matrix(), test_mass) );
	CHECK( compare_vectors(lin_elem_test.get_stiffness_matrix(), test_stiffness) );
	CHECK( compare_vectors(lin_elem_test.get_lumped_matrix(), test_lumped_mass) );
}

TEST_CASE( "Local matrices of linear line element 3D mesh", "[LocMatricesLine3D]" )
{
	LinearLineElement lin_elem_test({0,0,0, 1,1,1}, {0, 1});
	double length_test_elem = sqrt(3);

	std::vector<double> test_mass {(double) length_test_elem / 3, (double) length_test_elem / 6, (double) length_test_elem / 6, (double) length_test_elem / 3};
	std::vector<double> test_stiffness {1.0 / length_test_elem, - 1.0 / length_test_elem, - 1.0 / length_test_elem, 1.0 / length_test_elem};
	std::vector<double> test_lumped_mass {length_test_elem / 2.0, length_test_elem / 2.0};

	CHECK( compare_vectors(lin_elem_test.get_mass_matrix(), test_mass) );
	CHECK( compare_vectors(lin_elem_test.get_stiffness_matrix(), test_stiffness) );
	CHECK( compare_vectors(lin_elem_test.get_lumped_matrix(), test_lumped_mass) );
}

TEST_CASE( "Local matricees of linear triangle element 3D mesh", "[Local3DTriangleMatrices]" )
{
	// ------- Элемент задается на трёхмерной сетке
	// ------- 1 ------ {x0, y0, z0, x1, y1, z1, x2, y2, z2}
	LinearTriangleElement triang_elem_test( {0,0,0, 1,0,0, 0,1,0}, {0, 1, 2});

	CHECK( compare_vectors(triang_elem_test.get_mass_matrix(), 
	{1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0,
	 1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0,
	 1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0})
	);

	CHECK( compare_vectors(triang_elem_test.get_stiffness_matrix(),
	{1.0, -1.0 / 2.0, -1.0 / 2.0,
	 -1.0 / 2.0, 1.0 / 2.0, 0,
	 -1.0 / 2.0, 0, 1.0 / 2.0})
	);
}

TEST_CASE( "Local matrices of linear triangle elements 2D mesh", "[Local2DTriangleMatrices]" )
{
	// ------- Элемент задается на трёхмерной сетке
	// ------- 1 ------ {x0, y0, z0, x1, y1, z1, x2, y2, z2}
	LinearTriangleElement triang_elem_test( {0,0, 1,0, 0,1}, {0, 1, 2});

	CHECK( compare_vectors(triang_elem_test.get_mass_matrix(), 
	{1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0,
	 1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0,
	 1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0})
	);

	CHECK( compare_vectors(triang_elem_test.get_stiffness_matrix(),
	{1.0, -1.0 / 2.0, -1.0 / 2.0,
	 -1.0 / 2.0, 1.0 / 2.0, 0,
	 -1.0 / 2.0, 0, 1.0 / 2.0})
	);
}

TEST_CASE( "Local matrices of the linear point boundary element 1D mesh", "[PointBound1D]" )
{
	LinearPointBoundaryElement point_elem_test({0}, {0}, 0);

	CHECK( 
		compare_vectors(point_elem_test.get_mass_matrix(), {1})
	);

	CHECK(
		compare_vectors(point_elem_test.get_stiffness_matrix(), {0})
	);

	CHECK(
		compare_vectors(point_elem_test.get_lumped_matrix(), {1})
	);
}

TEST_CASE( "Local matrices of the linear point boundary element 2D mesh", "[PointBound2D]" )
{
	LinearPointBoundaryElement point_elem_test({0, 0}, {0}, 0);

	CHECK( 
		compare_vectors(point_elem_test.get_mass_matrix(), {1})
	);

	CHECK(
		compare_vectors(point_elem_test.get_stiffness_matrix(), {0})
	);

	CHECK(
		compare_vectors(point_elem_test.get_lumped_matrix(), {1})
	);
}

TEST_CASE( "Local matrices of the linear point boundary element 3D mesh", "[PointBound3D]" )
{
	LinearPointBoundaryElement point_elem_test({0, 0, 0}, {0}, 0);

	CHECK( 
		compare_vectors(point_elem_test.get_mass_matrix(), {1})
	);

	CHECK(
		compare_vectors(point_elem_test.get_stiffness_matrix(), {0})
	);

	CHECK(
		compare_vectors(point_elem_test.get_lumped_matrix(), {1})
	);
}

TEST_CASE( "Local matrices of the linear line boundary element 3D mesh", "[LineBound3D]" )
{
	LinearLineBoundaryElement line_elem_test({0,0,0, 1,1,1}, {0, 1});

	double length_test_elem = sqrt(3);

	std::vector<double> test_mass {(double) length_test_elem / 3, (double) length_test_elem / 6, (double) length_test_elem / 6, (double) length_test_elem / 3};
	std::vector<double> test_stiffness {1.0 / length_test_elem, - 1.0 / length_test_elem, - 1.0 / length_test_elem, 1.0 / length_test_elem};
	std::vector<double> test_lumped_mass {length_test_elem / 2.0, length_test_elem / 2.0};


	CHECK(
		compare_vectors(line_elem_test.get_stiffness_matrix(), test_stiffness)
	);

	CHECK(
		compare_vectors(line_elem_test.get_lumped_matrix(), test_lumped_mass)
	);

	CHECK(
		compare_vectors(line_elem_test.get_mass_matrix(), test_mass)
	);
}

TEST_CASE( "Local matrices of the linear line boundary element 2D", "[LineBound2D]" )
{
	LinearLineBoundaryElement line_elem_test({0,0, 1,1}, {0, 1});

	double length_test_elem = sqrt(2);

	std::vector<double> test_mass {(double) length_test_elem / 3, (double) length_test_elem / 6, (double) length_test_elem / 6, (double) length_test_elem / 3};
	std::vector<double> test_stiffness {1.0 / length_test_elem, - 1.0 / length_test_elem, - 1.0 / length_test_elem, 1.0 / length_test_elem};
	std::vector<double> test_lumped_mass {length_test_elem / 2.0, length_test_elem / 2.0};


	CHECK(
		compare_vectors(line_elem_test.get_stiffness_matrix(), test_stiffness)
	);

	CHECK(
		compare_vectors(line_elem_test.get_lumped_matrix(), test_lumped_mass)
	);

	CHECK(
		compare_vectors(line_elem_test.get_mass_matrix(), test_mass)
	);
}
