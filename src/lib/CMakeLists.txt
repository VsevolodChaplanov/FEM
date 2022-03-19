set(LIB_CPPS pde_lib/VectorOperations.cpp
	pde_lib/Builder.cpp
	pde_lib/CompressedM.cpp
	pde_lib/FemGrid.cpp
	pde_lib/FemPDE.cpp
	pde_lib/BoundaryElements.cpp
	pde_lib/LinElem.cpp
	pde_lib/Preconditioners.cpp
	pde_lib/Solvers.cpp
)

set(HEADERS headers/Builder.h
	headers/CompressedM.h
	headers/FemGrid.h
	headers/FemPDE.h
	headers/IBoundaryElem.h
	headers/IFiniteElem.h
	headers/Preconditioners.h
	headers/Solvers.h
)

add_library(FEM_PDE_lib STATIC ${LIB_CPPS} ${HEADERS})
# target_include_directories(${FEMLIB} PRIVATE "${PROJECT_SOURCE_DIR}/lib/pde_lib"
# 	"${PROJECT_SOURCE_DIR}/lib/headers")

