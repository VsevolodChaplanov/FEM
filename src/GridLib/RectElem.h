#ifndef __RECTANGULAR_ELEMENT__
#define __RECTANGULAR_ELEMENT__

#include <cmath>
#include <vector>
#include "IFiniteElem.h"

class RectElem : public IFiniteElement
{
public:

	RectElem(const std::vector<double> &vertices, const std::vector<std::size_t> &GIndices);
	std::vector<std::vector<double>>* GetMass(const std::size_t , const std::size_t ) override;
	std::vector<std::vector<double>>* GetStiffness(const std::size_t , const std::size_t ) override;
	std::vector<double>* GetLumped(const std::size_t, const std::size_t ) override;
	std::vector<double> PhysToParam(const double x, const double y = 0, const double z = 0) override;
	std::vector<double> PhysToParam(const std::vector<double> &) override;

private:

	double phi1(const double &, const double &);
	double phi2(const double &, const double &);
	double phi3(const double &, const double &);
	double phi4(const double &, const double &);

private:

	const double Lx;
	const double Ly;
	const double Volume = Lx * Ly;
	const std::vector<std::size_t> GlobalIndices;

public:

	std::vector<std::vector<double>> MassMatrix {
														{Lx * Ly / 36 * 4, Lx * Ly / 36 * 2, Lx * Ly / 36 * 1, Lx * Ly / 36 * 2},
													    {Lx * Ly / 36 * 2, Lx * Ly / 36 * 4, Lx * Ly / 36 * 2, Lx * Ly / 36 * 1},
													    {Lx * Ly / 36 * 1, Lx * Ly / 36 * 2, Lx * Ly / 36 * 4, Lx * Ly / 36 * 2},
													    {Lx * Ly / 36 * 2, Lx * Ly / 36 * 1, Lx * Ly / 36 * 2, Lx * Ly / 36 * 4} 
													  };
	std::vector<std::vector<double>> StiffnessMatrix {
															{
															 1 / (6 * Lx * Ly) * 2 * (Ly * Ly + Lx * Lx), 
															 1 / (6 * Lx * Ly) * (Lx * Lx - 2 * Ly * Ly), 
															 1 / (6 * Lx * Ly) * (- Ly * Ly - Lx * Lx), 
															 1 / (6 * Lx * Ly) * (Ly * Ly - 2 * Lx * Lx)
															}, 
															{
															 1 / (6 * Lx * Ly) * (Lx * Lx - 2 * Ly * Ly),
															 1 / (6 * Lx * Ly) * 2 * (Ly * Ly + Lx * Lx),
															 1 / (6 * Lx * Ly) * (Ly * Ly - 2 * Lx * Lx),
															 1 / (6 * Lx * Ly) * (- Ly * Ly - Lx * Lx)
															},
															{
															 1 / (6 * Lx * Ly) * (- Ly * Ly - Lx * Lx),
															 1 / (6 * Lx * Ly) * (Ly * Ly - 2 * Lx * Lx),
															 1 / (6 * Lx * Ly) * 2 * (Ly * Ly + Lx * Lx),
															 1 / (6 * Lx * Ly) * (Lx * Lx - 2 * Ly * Ly)
															},
															{
															 1 / (6 * Lx * Ly) * (Ly * Ly - 2 * Lx * Lx),
															 1 / (6 * Lx * Ly) * (- Ly * Ly - Lx * Lx),
															 1 / (6 * Lx * Ly) * (Lx * Lx - 2 * Ly * Ly),
														 	 1 / (6 * Lx * Ly) * 2 * (Ly * Ly + Lx * Lx)
															} 
														   };
	std::vector<double> LumpedMatrix { 1/4, 1/4, 1/4, 1/4 };
};

#endif
