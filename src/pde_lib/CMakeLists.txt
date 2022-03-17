set(LIB_CPPS Builder.cpp
	CompressedM.cpp
	#CubeElem.cpp
	FemGrid.cpp
	FemPDE.cpp
	IBoundaryElem.cpp
	LinElem.cpp
	Preconditioners.cpp
	#RectElem.cpp
	Solvers.cpp
	VecMath.cpp
	VectorOperations.cpp
)

add_library(FEM_lib STATIC ${LIB_CPPS})