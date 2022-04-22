set(SRC PDE_lib_tests.cpp
	MatVec_Math_tests.cpp
	FiniteElements_tests.cpp
	LocalMat_tests.cpp
	GlobalMat_tests.cpp
	MatrixSolvers_tests.cpp
	FiniteElemMesh_tests.cpp
	Preconditioners_test.cpp
	main.cpp)
set(HEADERS PDE_lib_tests.h)

set(LIB_TESTS fem_pde_lib-test)

add_executable(${LIB_TESTS} ${SRC} ${HEADERS})

target_link_libraries(${LIB_TESTS} PRIVATE ${FEM_PDE_lib} ${CATCH2_FRAMEWORK})
target_include_directories(${LIB_TESTS} PRIVATE ${CMAKE_SOURCE_DIR}/src/lib/include)

