#include <iostream>
#include "GridLib/Builder1D.cpp"
#include "PDE/FEM_PDE.h"
#include "PDE/FEM_PDE.cpp"

double f(double x)
{
	return 1;
}

int main(int argc, char const *argv[])
{
	std::vector<double> f_vec(4, 1.);
	
	Builder1D* FirstElem = new Builder1D(1,5, 4);
	FEM_PDE* My_Pde = new FEM_PDE(4);
	My_Pde->AssembleSystem(FirstElem, f_vec);
	My_Pde->Solve("GD");
	return 0;
}
