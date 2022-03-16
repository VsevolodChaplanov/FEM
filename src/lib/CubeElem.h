#ifndef __CUBE_ELEMENT__
#define __CUBE_ELEMENT__


#include <vector>
#include <cmath>
#include "IFiniteElem.h"

class CubeElem : public IFiniteElement
{
public:

	CubeElem(const std::vector<double> &vertices, const std::vector<std::size_t> &GIndices);
	std::vector<std::vector<double>>* GetMass(const std::size_t , const std::size_t ) override;
	std::vector<std::vector<double>>* GetStiffness(const std::size_t , const std::size_t ) override;
	std::vector<double>* GetLumped(const std::size_t, const std::size_t ) override;
	std::vector<double> PhysToParam(const double x, const double y = 0, const double z = 0) override;
	std::vector<double> PhysToParam(const std::vector<double> &) override;

private:

	double phi1(const double _xi, const double _eta, const double _zeta); // Matches to point 0
	double phi2(const double _xi, const double _eta, const double _zeta); // Matches to point 1
	double phi3(const double _xi, const double _eta, const double _zeta); // Matches to point 3
	double phi4(const double _xi, const double _eta, const double _zeta); // Matches to point 2
	double phi5(const double _xi, const double _eta, const double _zeta); // Matches to point 4
	double phi6(const double _xi, const double _eta, const double _zeta); // Matches to point 5
	double phi7(const double _xi, const double _eta, const double _zeta); // Matches to point 7
	double phi8(const double _xi, const double _eta, const double _zeta); // Matches to point 6

private:

	const double Lx;
	const double Ly;
	const double Lz;
	const double DetJ = Lx * Ly * Lz;

public:

	//
	// TODO 
	// 		Либо стоит хранить эти матрицы в статическом массиве
	//		Либо хранить только нижний или верхний треугольные подматрицы т.к. симметричные 


	std::vector<std::size_t> GlobalIndices;
	const std::vector<std::vector<double>> MassMatrix { {DetJ / 216 * 8, DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 1, DetJ / 216 * 2},
														{DetJ / 216 * 4, DetJ / 216 * 8, DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 1},
														{DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 8, DetJ / 216 * 4, DetJ / 216 * 1, DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 2},
														{DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 8, DetJ / 216 * 2, DetJ / 216 * 1, DetJ / 216 * 2, DetJ / 216 * 4},
														{DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 1, DetJ / 216 * 2, DetJ / 216 * 8, DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 4},
														{DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 1, DetJ / 216 * 4, DetJ / 216 * 8, DetJ / 216 * 4, DetJ / 216 * 2},
														{DetJ / 216 * 1, DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 8, DetJ / 216 * 4},
														{DetJ / 216 * 2, DetJ / 216 * 1, DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 8} };

	const std::vector<std::vector<double>> StiffnessMatrix {	{	(Ly*Lz)/(9*Lx)+(Lx*(Ly*Ly+Lz*Lx))/(9*Ly*Lz),
																	-((Ly*Lz)/(9*Lx))+(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	(Lx*Ly)/(36*Lz)-(Lx*Lz)/(18*Ly)-(Ly*Lz)/(18*Lx),
																	(Lx*Ly)/(18*Lz)-(Lx*Lz)/(9*Ly)+(Ly*Lz)/(18*Lx),
																	-((Lx*Ly)/(9*Lz))+(Lx*Lz)/(18*Ly)+(Ly*Lz)/(18*Lx),
																	-((Lx*Ly)/(18*Lz))+(Lx*Lz)/(36*Ly)-(Ly*Lz)/(18*Lx),
																	-((Ly*Lz)/(36*Lx))-(Lx*(Ly*Ly+Lz*Lz))/(36*Ly*Lz), 
																	(Ly*Lz)/(36*Lx)-(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz)
																},
																{	-((Ly*Lz)/(9*Lx))+(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	(Ly*Lz)/(9*Lx)+(Lx*(Ly*Ly+Lz*Lz))/(9*Ly*Lz),
																	(Lx*Ly)/(18*Lz)-(Lx*Lz)/(9*Ly)+(Ly*Lz)/(18*Lx),
																	(Lx*Ly)/(36*Lz)-(Lx*Lz)/(18*Ly)-(Ly*Lz)/(18*Lx),
																	-((Lx*Ly)/(18*Lz))+(Lx*Lz)/(36*Ly)-(Ly*Lz)/(18*Lx),
																	-((Lx*Ly)/(9*Lz))+(Lx*Lz)/(18*Ly)+(Ly*Lz)/(18*Lx),
																	(Ly*Lz)/(36*Lx)-(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	-((Ly*Lz)/(36*Lx))-(Lx*(Ly*Ly+Lz*Lz))/(36*Ly*Lz)
																},
																{	(Lx*Ly)/(36*Lz)-(Lx*Lz)/(18*Ly)-(Ly*Lz)/(18*Lx),
																	(Lx*Ly)/(18*Lz)-(Lx*Lz)/(9*Ly)+(Ly*Lz)/(18*Lx),
																	(Ly*Lz)/(9*Lx)+(Lx*(Ly*Ly+Lz*Lz))/(9*Ly*Lz),
																	-((Ly*Lz)/(9*Lx))+(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	-((Ly*Lz)/(36*Lx))-(Lx*(Ly*Ly+Lz*Lz))/(36*Ly*Lz),
																	(Ly*Lz)/(36*Lx)-(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	-((Lx*Ly)/(9*Lz))+(Lx*Lz)/(18*Ly)+(Ly*Lz)/(18*Lx),
																	-((Lx*Ly)/(18*Lz))+(Lx*Lz)/(36*Ly)-(Ly*Lz)/(18*Lx)
																},
																	{(Lx*Ly)/(18*Lz)-(Lx*Lz)/(9*Ly)+(Ly*Lz)/(18*Lx),
																	(Lx*Ly)/(36*Lz)-(Lx*Lz)/(18*Ly)-(Ly*Lz)/(18*Lx),
																	-((Ly*Lz)/(9*Lx))+(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	(Ly*Lz)/(9*Lx)+(Lx*(Ly*Ly+Lz*Lz))/(9*Ly*Lz),
																	(Ly*Lz)/(36*Lx)-(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	-((Ly*Lz)/(36*Lx))-(Lx*(Ly*Ly+Lz*Lz))/(36*Ly*Lz),
																	-((Lx*Ly)/(18*Lz))+(Lx*Lz)/(36*Ly)-(Ly*Lz)/(18*Lx),
																	-((Lx*Ly)/(9*Lz))+(Lx*Lz)/(18*Ly)+(Ly*Lz)/(18*Lx)
																},
																{	-((Lx*Ly)/(9*Lz))+(Lx*Lz)/(18*Ly)+(Ly*Lz)/(18*Lx),
																	-((Lx*Ly)/(18*Lz))+(Lx*Lz)/(36*Ly)-(Ly*Lz)/(18*Lx),
																	-((Ly*Lz)/(36*Lx))-(Lx*(Ly*Ly+Lz*Lz))/(36*Ly*Lz),
																	(Ly*Lz)/(36*Lx)-(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	(Ly*Lz)/(9*Lx)+(Lx*(Ly*Ly+Lz*Lz))/(9*Ly*Lz),
																	-((Ly*Lz)/(9*Lx))+(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	(Lx*Ly)/(36*Lz)-(Lx*Lz)/(18*Ly)-(Ly*Lz)/(18*Lx),
																	(Lx*Ly)/(18*Lz)-(Lx*Lz)/(9*Ly)+(Ly*Lz)/(18*Lx)
																},
																{	-((Lx*Ly)/(18*Lz))+(Lx*Lz)/(36*Ly)-(Ly*Lz)/(18*Lx),
																	-((Lx*Ly)/(9*Lz))+(Lx*Lz)/(18*Ly)+(Ly*Lz)/(18*Lx),
																	(Ly*Lz)/(36*Lx)-(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	-((Ly*Lz)/(36*Lx))-(Lx*(Ly*Ly+Lz*Lz))/(36*Ly*Lz),
																	-((Ly*Lz)/(9*Lx))+(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	(Ly*Lz)/(9*Lx)+(Lx*(Ly*Ly+Lz*Lz))/(9*Ly*Lz),
																	(Lx*Ly)/(18*Lz)-(Lx*Lz)/(9*Ly)+(Ly*Lz)/(18*Lx),
																	(Lx*Ly)/(36*Lz)-(Lx*Lz)/(18*Ly)-(Ly*Lz)/(18*Lx)
																},
																{	-((Ly*Lz)/(36*Lx))-(Lx*(Ly*Ly+Lz*Lz))/(36*Ly*Lz),
																	(Ly*Lz)/(36*Lx)-(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	-((Lx*Ly)/(9*Lz))+(Lx*Lz)/(18*Ly)+(Ly*Lz)/(18*Lx),
																	-((Lx*Ly)/(18*Lz))+(Lx*Lz)/(36*Ly)-(Ly*Lz)/(18*Lx),
																	(Lx*Ly)/(36*Lz)-(Lx*Lz)/(18*Ly)-(Ly*Lz)/(18*Lx),
																	(Lx*Ly)/(18*Lz)-(Lx*Lz)/(9*Ly)+(Ly*Lz)/(18*Lx),
																	(Ly*Lz)/(9*Lx)+(Lx*(Ly*Ly+Lz*Lz))/(9*Ly*Lz),
																	-((Ly*Lz)/(9*Lx))+(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz)
																},
																{	(Ly*Lz)/(36*Lx)-(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	-((Ly*Lz)/(36*Lx))-(Lx*(Ly*Ly+Lz*Lz))/(36*Ly*Lz),
																	-((Lx*Ly)/(18*Lz))+(Lx*Lz)/(36*Ly)-(Ly*Lz)/(18*Lx),
																	-((Lx*Ly)/(9*Lz))+(Lx*Lz)/(18*Ly)+(Ly*Lz)/(18*Lx),
																	(Lx*Ly)/(18*Lz)-(Lx*Lz)/(9*Ly)+(Ly*Lz)/(18*Lx),
																	(Lx*Ly)/(36*Lz)-(Lx*Lz)/(18*Ly)-(Ly*Lz)/(18*Lx),
																	-((Ly*Lz)/(9*Lx))+(Lx*(Ly*Ly+Lz*Lz))/(18*Ly*Lz),
																	(Ly*Lz)/(9*Lx)+(Lx*(Ly*Ly+Lz*Lz))/(9*Ly*Lz)
																}
															};
	const std::vector<double> LumpedMatrix {DetJ / 8, DetJ / 8,  DetJ / 8, DetJ / 8, DetJ / 8, DetJ / 8, DetJ / 8, DetJ / 8};
};


#endif