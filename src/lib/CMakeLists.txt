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
)


add_library(${FEM_PDE_lib} STATIC ${LIB_CPPS})
target_include_directories(${FEM_PDE_lib} PRIVATE headers)
