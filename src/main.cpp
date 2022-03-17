#include <iostream>
#include "pde_lib/Builder.cpp"
#include "pde_lib/FemGrid.cpp"
#include "pde_lib/FemPDE.cpp"
#include "pde_lib/VecMath.cpp"

double f_fun(double* point)
{
	return 1;
}

double k_fun(double* point)
{
	return 1;
}

double u_ex(double * point)
{
	return 1;
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
	fempde.apply_boundary_condition_dirichlet(u_ex, femgridlinear.boundary_element_indices(1));
	fempde.apply_boundary_condition_dirichlet(u_ex, femgridlinear.boundary_element_indices(2));
	// 	fempde->apply_boundary_condition_dirichlet(u_exact, femgridlinear->boundary_element_indices()); All

	std::vector<double> sol = fempde.solve("GD"); // вектор решения ин в самом методе

	std::vector<double> u_ex_num = femgridlinear.approximate(u_ex);
	
	std::vector<double> difference = vector_difference(sol, u_ex_num);

	// Было double norm2 = femgridlinear.norm2(difference);
	// Но я сделал отдельную функцию для такой второй нормы по соображениям, что вторая норма может быть вычислена вроде для любого такого вектора
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