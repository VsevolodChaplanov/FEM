#ifndef __RECTANGULAR_ELEMENT__
#define __RECTANGULAR_ELEMENT__

#include <cmath>
#include <vector>
#include "../Math/VecLib.h"

class RectElem
{
public:

	RectElem(const std::vector<std::vector<double>> &vertices, const std::vector<std::size_t> &GIndices);
	double GetMass(const std::size_t &, const std::size_t &);
	double GetSteffness(const std::size_t &, const std::size_t &);
	double GetLumped(const std::size_t &, const std::size_t &);
	void P_to_Param(const double &, const double &);

private:

	double phi1(const double &, const double &);
	double phi2(const double &, const double &);
	double phi3(const double &, const double &);
	double phi4(const double &, const double &);

public:

	const std::vector<std::vector<double>> MassMatrix {
														{Lx * Ly / 36 * 4, Lx * Ly / 36 * 2, Lx * Ly / 36 * 1, Lx * Ly / 36 * 2},
													    {Lx * Ly / 36 * 2, Lx * Ly / 36 * 4, Lx * Ly / 36 * 2, Lx * Ly / 36 * 1},
													    {Lx * Ly / 36 * 1, Lx * Ly / 36 * 2, Lx * Ly / 36 * 4, Lx * Ly / 36 * 2},
													    {Lx * Ly / 36 * 2, Lx * Ly / 36 * 1, Lx * Ly / 36 * 2, Lx * Ly / 36 * 4} 
													  };
	const std::vector<std::vector<double>> StiffnessMatrix {
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
	const std::vector<double> LumpedMatrix { 1/4, 1/4, 1/4, 1/4 };

private:

	const double Lx;
	const double Ly;
	const double Volume = Lx * Ly;
	const std::vector<std::size_t> GlobalIndices;
};

#endif
