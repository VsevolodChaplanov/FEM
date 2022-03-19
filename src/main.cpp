#include <iostream>

#include "lib/pde_lib/Builder.cpp"

double u_ex(double * point)
{
	return sin(9 * (point[0]+0.2) * (point[0]+0.2));
}

double f_fun(double* point)
{
	return (point[0]+1.2) * u_ex(point);
}

double k_fun(double* point)
{
	return 1.0 / (18 * 18 * (point[0] + 0.2));
}

int main(int argc, char const *argv[])
{
	// Builder mathes to left boundary 1 to right boundary 2
	FemGrid femgridlinear = Builder::BuildLinear1DGrid(0, 1, 100); // Для ГУ на границе стоит элемент порядка ниже
	// femgtidlinear->assign_boundary_type(selector , type)

	FemPDE fempde(&femgridlinear, f_fun, k_fun); 
	// fempde->set_bc() 
	fempde.assemble();
	fempde.apply_boundary_condition_dirichlet(u_ex, femgridlinear.boundary_element_indices(1));
	fempde.apply_boundary_condition_dirichlet(u_ex, femgridlinear.boundary_element_indices(2));
	// 	fempde->apply_boundary_condition_dirichlet(u_exact, femgridlinear->boundary_element_indices()); All

	std::vector<double> sol = fempde.solve("Thomas"); // вектор решения ин в самом методе

	// Approximate analytical function along calculation area
	std::vector<double> u_ex_num = femgridlinear.approximate(u_ex);
	
	// Obtain difference between vectors of numerical solution and analytical function approximation
	std::vector<double> difference = VDiff(sol, u_ex_num);

	// Calculate the maximum deviation 
	double normmax = max_abs(difference);
	// Calculate second norm
	double norm2 = norm_2(difference);

	// Calculate the second norm within the approximation
	double norm2a = femgridlinear.norm2(difference);


	// Save obtained data as .vtk format
	femgridlinear.savevtk(sol, "Test_solution.vtk");

	// Get params
	std::cout << "Maximun absolute diviation: " << normmax << std::endl;
	std::cout << "Second norm: " << norm2 << std::endl;
	std::cout << "Second norm within approximation: " << norm2a << std::endl;

	return 0;
}
