#ifndef __FEM_PDE_LIB_TESTS__
#define __FEM_PDE_LIB_TESTS__

#include <iostream>
#include "FemGrid.h"
#include "FemPDE.h"
#include "SolverParams.h"
#include "Builder.h"
#include "IFiniteElement.h"
#include "LinearLineElement.h"
#include "LinearTriangleElement.h"
#include "LinearLineBoundaryElement.h"
#include "GlobalAssemblers.h"
#include "VectorOperations.h"
#include "CompressedM.h"
#include "FiniteElemMeshParser.h"
#include <catch2/catch.hpp>

double u_ex(const double* point);
double f_fun(const double* point);
double k_fun(const double* point);

#endif 