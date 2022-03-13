include_directories(GridLib)
include_directories(PDE)
include_directories(Math)
include_directories(Solvers)
include_directories(CompressedMatrix)

set(SRC main.cpp)

add_executable(${TEST_EXE} ${SRC})		  
add_subdirectory("GridLib/")
add_subdirectory("PDE/")
add_subdirectory("Math/")
add_subdirectory("Solvers/")
add_subdirectory("CompressedMatrix/")


target_link_libraries(${TEST_EXE} 	FINELEM
									VECT_LIB
									COMPRESSEDMAT_LIB
									PDE_LIB
									VECT_LIB
									SOLVERS_LIB)

