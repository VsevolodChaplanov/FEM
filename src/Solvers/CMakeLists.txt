set(SOLVERS Solvers.cpp
				Preconditioners.cpp
				VectorOperations.cpp)

add_library(SOLVERS_LIB STATIC ${SOLVERS})