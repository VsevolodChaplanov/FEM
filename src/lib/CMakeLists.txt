set(LIB_CPPS pde_lib/Builder.cpp
	pde_lib/CompressedM.cpp
	pde_lib/FemGrid.cpp
	pde_lib/FemPDE.cpp
	pde_lib/IFiniteElement.cpp
	pde_lib/IBoundaryElement.cpp
	pde_lib/LinearLineElement.cpp
	pde_lib/LinearPointBoundaryElement.cpp
	pde_lib/SolverParams.cpp
	pde_lib/Preconditioners.cpp
	pde_lib/Solvers.cpp
	pde_lib/VectorOperations.cpp
	pde_lib/GlobalAssemblers.cpp
	pde_lib/FiniteElemMeshParser.cpp
	pde_lib/LinearTriangleElement.cpp
	pde_lib/LinearLineBoundaryElement.cpp)

set(HEADERS include/Builder.h
	include/CompressedM.h
	include/FemGrid.h
	include/FemPDE.h
	include/IFiniteElement.h
	include/IBoundaryElement.h
	include/LinearLineElement.h
	include/LinearPointBoundaryElement.h
	include/SolverParams.h
	include/Preconditioners.h
	include/Solvers.h
	include/VectorOperations.h
	include/GlobalAssemblers.h
	include/FiniteElemMeshParser.h
	include/LinearTriangleElement.h
	include/LinearLineBoundaryElement.h
)

set(FEM_PDE_lib fem_lib)

find_package(Boost REQUIRED)

add_library(${FEM_PDE_lib} SHARED ${LIB_CPPS} ${HEADERS})
target_include_directories(${FEM_PDE_lib} PRIVATE include)
target_include_directories(${FEM_PDE_lib} PRIVATE ${Boost_INCLUDE_DIR})
target_link_libraries(${FEM_PDE_lib} PRIVATE ${Boost_LIBARARIES})
target_compile_options(${FEM_PDE_lib} PRIVATE -Wall)

add_subdirectory(tests)