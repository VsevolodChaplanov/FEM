#include <iostream>
#include "lib/Builder.h"
#include "lib/FemGrid.h"
#include "lib/FemPDE.h"

double f_fun(double* point)
{
	return 0;
}

double k_fun(double* point)
{
	return 0;
}

double u_ex(double * point)
{

}

// Selector variant
bool left_boundaries(double* center_face) 
{
	if(center_face[0] <= 1e-8)
	{
		return true;
	}
	return false;
}

int main(int argc, char const *argv[])
{
	// Builder mathes to left boundary 1 to right boundary 2
	FemGrid femgridlinear = Builder::BuildLinear1DGrid(0, 1, 10); // Для ГУ на границе стоит элемент порядка ниже
	// femgtidlinear->assign_boundary_type(selector , type)

	FemPDE fempde(&femgridlinear, f_fun, k_fun); 
	// fempde->set_bc() 
	fempde.assemble();
	fempde.apply_boundary_condition_dirichlet(u_exact, femgridlinear.boundary_element_indices(1));
	fempde.apply_boundary_condition_dirichlet(u_exact, femgridlinear.boundary_element_indices(2));
	// 	fempde->apply_boundary_condition_dirichlet(u_exact, femgridlinear->boundary_element_indices()); All

	std::vector<double> sol = fempde.solve("GD"); // вектор решения ин в самом методе

	std::vector<double> u_ex_num = femgridlinear.approximate(u_ex);
	
	std::vector<double> difference = u_ex_num - sol;

	double norm2 = femgridlinear.norm2(difference);
	double normmax = max_abs(difference);


	femgridlinear.savevtk(sol, "Test_solution.vtk");

	return 0;
}

//
//
// 
//
// boundary_element_indices - есть селектор