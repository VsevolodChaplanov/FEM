#ifndef __FEM_PDE_LIB_TESTS__
#define __FEM_PDE_LIB_TESTS__

#include <iostream>
#include <catch2/catch_all.hpp>

#include "FemGrid.h"
#include "FemPDE.h"
#include "SolverParams.h"
#include "Builder.h"
#include "IFiniteElem.h"
#include "LinElem.h"
#include "GlobalAssemblers.h"
#include "VectorOperations.h"
#include "CompressedM.h"
#include "FiniteElemMeshParser.h"

double u_ex(const double* point);
double f_fun(const double* point);
double k_fun(const double* point);

#endif 