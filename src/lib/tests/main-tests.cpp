#include <vector>
#include <limits>
#include "catch-tests.h"

int main(int argc, char const *argv[])
{
	int result = Catch::Session().run(argc, argv);

	// Builder mathes to left boundary 1 to right boundary 2
	FemGrid femgridlinear = Builder::BuildLinear1DGrid(0, 1, 10); // Для ГУ на границе стоит элемент порядка ниже
	// femgtidlinear->assign_boundary_type(selector , type)

	MatrixSolverParams* params = new MatrixSolverParams(MatrixSolverParams::Methods::Thomas, MatrixSolverParams::Preconditioners::None, 1000, 1.e-5, 10);

	FemPDE fempde(&femgridlinear, f_fun, k_fun, params); 
	FemPDE fempde_test(&femgridlinear, f_fun, k_fun, params); 

	fempde.assemble();
	fempde_test.new_assembler();


	fempde.apply_boundary_condition_dirichlet(u_ex, femgridlinear.boundary_element_indices(1));
	fempde.apply_boundary_condition_dirichlet(u_ex, femgridlinear.boundary_element_indices(2));
	// 	fempde->apply_boundary_condition_dirichlet(u_exact, femgridlinear->boundary_element_indices()); All

	std::vector<double> sol = fempde.solve(); // вектор решения ин в самом методе

	// Approximate analytical function along calculation area
	std::vector<double> u_ex_num = femgridlinear.approximate(u_ex);
	
	// Obtain difference between vectors of numerical solution and analytical function approximation
	std::vector<double> difference = vector_diff(sol, u_ex_num);

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

	delete params;
	return 0;
}