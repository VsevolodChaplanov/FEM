#include "PDE_lib_tests.h"

// Exist test of the following methods
// 		Vectors multyply (dot product)
//		Matrix multyply
// 		Matrix-vector multyply
//		Vector sum
// 		Matrix sum
// 		Vectors equality
// 		Vector length
// 		Vector centre

TEST_CASE( "Matrix-vector multyply", "[MatVecMult]" )
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

	CHECK( test_result[0] == 6.5 );
	CHECK( test_result[1] == 4.5);
	CHECK( test_result[2] == 11.5 );
	CHECK( test_result[3] == 0 );	
}

TEST_CASE( "Sum of two vectors", "[VectorsSum]")
{
	std::vector<double> first_ex {1.5 , 0.5, 1.5, -1.5};
	std::vector<double> second_ex {0.5, -1.5, 0.5, -0.5};
	std::vector<double> sum = vector_sum(first_ex, second_ex);
	
	CHECK( sum[0] == 2 );
	CHECK( sum[1] == -1 );
	CHECK( sum[2] == 2 );
	CHECK( sum[3] == -2 );
}

TEST_CASE( "Dot product of two vectors", "[VectorsDotProduct]")
{
	std::vector<double> first_ex {1.5 , 0.5, 1.5, -1.5};
	std::vector<double> second_ex {0.5, -1.5, 0.5, -0.5};
	double sc_test = dot_product(first_ex, second_ex);
	
	CHECK( sc_test == 1.5 );
}

TEST_CASE( "Sum of the two sparse matrices", "[SparseMatSum]" )
{
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

TEST_CASE( "Sparse matrix multyply", "[SparseMatMult]" )
{
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

TEST_CASE( "Equality of two vectors", "[VectorEquality]" )
{
	std::vector<double> first_d_test {0.5, -0.1, 0.5};
	std::vector<double> second_d_test {0.5, -0.1, 0.5};
	std::vector<double> thrid_d_test {0.5, -0.1, 0.5, -0.8};

	CHECK_FALSE( 
		compare_vectors(first_d_test, thrid_d_test) 
	);
	CHECK(
		compare_vectors(first_d_test, second_d_test)
	);

	std::vector<size_t> first_ui_test {5, 1, 5};
	std::vector<size_t> second_ui_test {5, 1, 5};
	std::vector<size_t> thrid_ui_test {5, 1, 5, 8};

	CHECK_FALSE( 
		compare_vectors(first_ui_test, thrid_ui_test) 
	);
	CHECK(
		compare_vectors(first_ui_test, second_ui_test)
	);
}

TEST_CASE( "Lenght of the vector", "[VectorLength]" )
{
	std::vector<double> first_test {0.0, 0.5};
	std::vector<double> second_test {0.0, 0.0, 0.5, 0.5};
	std::vector<double> thrid_test {0.0, 0.0, 0.0, 0.5, 0.5, 0.5};

	CHECK( 
		vector_lenght(first_test) == 0.5
	);
	CHECK( 
		vector_lenght(second_test) == Approx(sqrt(0.5*0.5+0.5*0.5))
	);
	CHECK( 
		vector_lenght(thrid_test) == Approx(sqrt(0.5*0.5+0.5*0.5+0.5*0.5))
	);
}