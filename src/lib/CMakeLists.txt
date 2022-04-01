set(LIB_CPPS pde_lib/Builder.cpp
	pde_lib/CompressedM.cpp
	pde_lib/FemGrid.cpp
	pde_lib/FemPDE.cpp
	pde_lib/IFiniteElem.cpp
	pde_lib/IBoundaryElem.cpp
	pde_lib/LinElem.cpp
	pde_lib/PointBoundElem.cpp
	pde_lib/SolverParams.cpp
	pde_lib/Preconditioners.cpp
	pde_lib/Solvers.cpp
	pde_lib/VectorOperations.cpp
	pde_lib/GlobalAssemblers.cpp
)

set(HEADERS include/Builder.h
	include/CompressedM.h
	include/FemGrid.h
	include/FemPDE.h
	include/IFiniteElem.h
	include/IBoundaryElem.h
	include/LinElem.h
	include/PointBoundElem.h
	include/SolverParams.h
	include/Preconditioners.h
	include/Solvers.h
	include/VectorOperations.h
	include/GlobalAssemblers.h
)


add_library(${FEM_PDE_lib} STATIC ${LIB_CPPS} ${HEADERS})
target_include_directories(${FEM_PDE_lib} PRIVATE include)

add_subdirectory(tests)