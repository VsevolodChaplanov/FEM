#ifndef __CUBE_ELEMENT__
#define __CUBE_ELEMENT__


#include <vector>
#include <cmath>

struct Point
{
public:

	double x = 0;
	double y = 0;
	double z = 0;

public:

	Point GetPoint();
	double GetX();
	double GetY();
	double GetZ();
	Point(const double &, const double &, const double &);
	Point(const std::vector<double> &);
};


class CubeElem
{
public:

	CubeElem(const std::vector<std::vector<double>> &vertices, const std::vector<std::size_t> &GIndices);
	double GetMass(const std::size_t &, const std::size_t &, const std::size_t &);
	double GetStiff(const std::size_t &, const std::size_t &, const std::size_t &);
	double GetLumpedMass(const std::size_t &, const std::size_t &, const std::size_t &);
	double P_to_Param(const double &, const double &, const double &);

private:

	double phi1();
	double phi2();
	double phi3();
	double phi4();
	double phi5();
	double phi6();
	double phi7();
	double phi8();

public:

	//
	// TODO 
	// 		Либо стоит хранить эти матрицы в статическом массиве
	//		Либо хранить только нижний или верхний треугольные подматрицы т.к. симметричные 


	std::vector<std::size_t> GlobalIndices;
	const std::vector<std::vector<double>> MassMatrix { {DetJ / 216 * 8, DetJ / 216 * 4, DetJ / 216 * 4, DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 2, DetJ / 216 * 2, DetJ / 216 * 2},
														{DetJ / 216 * 4, DetJ / 216 * 8, DetJ / 216 * 2, DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 1, DetJ / 216 * 4, DetJ / 216 * 4},
														{DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 8, DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 4, DetJ / 216 * 1, DetJ / 216 * 4},
														{DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 2, DetJ / 216 * 8, DetJ / 216 * 1, DetJ / 216 * 4, DetJ / 216 * 4, DetJ / 216 * 1},
														{DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 4, DetJ / 216 * 1, DetJ / 216 * 8, DetJ / 216 * 2, DetJ / 216 * 2, DetJ / 216 * 8},
														{DetJ / 216 * 2, DetJ / 216 * 1, DetJ / 216 * 4, DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 8, DetJ / 216 * 2, DetJ / 216 * 2},
														{DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 1, DetJ / 216 * 4, DetJ / 216 * 2, DetJ / 216 * 2, DetJ / 216 * 8, DetJ / 216 * 2},
														{DetJ / 216 * 2, DetJ / 216 * 4, DetJ / 216 * 4, DetJ / 216 * 1, DetJ / 216 * 8, DetJ / 216 * 2, DetJ / 216 * 2, DetJ / 216 * 8} };
	const std::vector<std::vector<double>> StiffnessMatrix { 	{	1 / (36 * DetJ) * 4 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly + Lz + Lz)), 
																	1 / (36 * DetJ) * (-4 * Ly * Ly * Lz * Lz + 2 * Lx * Lx * (Ly * Ly + Lz * Lz)), 
																	1 / (36 * DetJ) * (-4 * Lx * Lx * Ly * Ly + 2 * (Lx * Lx + Ly * Ly) * Lz * Lz), 
																	1 / (36 * DetJ) * 2 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly - 2 * Lz * Lz)), 
																	1 / (36 * DetJ) * (-2 * Ly * Ly * Lz * Lz + Lx * Lx * (-2 * Ly * Ly + Lz * Lz)), 
																	1 / (36 * DetJ) * (Ly * Ly * Lz * Lz - 2 * Lx * Lx * (Ly * Ly + Lz * Lz)), 
																	1 / (36 * DetJ) * (Lx * Lx * Ly * Ly - 2 * (Lx * Lx + Ly * Ly) * Lz * Lz), 
																	1 / (36 * DetJ) * (-2 * Ly * Ly * Lz * Lz + Lx * Lx * (-2 * Ly * Ly + Lz * Lz))	},
																{	1 / (36 * DetJ) * (-4 * Ly * Ly * Lz * Lz + 2 * Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * 4 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly + Lz + Lz)),
																	1 / (36 * DetJ) * (-2 * Ly * Ly * Lz * Lz + Lx * Lx * (-2 * Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (Lx * Lx * Ly * Ly - 2 * (Lx * Lx + Ly * Ly) * Lz * Lz),
																	1 / (36 * DetJ) * (-4 * Lx * Lx * Ly * Ly + 2 * (Lx * Lx + Ly * Ly) * Lz * Lz),
																	1 / (36 * DetJ) * (- Ly * Ly * Lz * Lz - Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * 2 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly - 2 * Lz * Lz)),
																	1 / (36 * DetJ) * (-4 * Lx * Lx * Ly * Ly + 2 * (Lx * Lx + Ly * Ly) * Lz * Lz)},
																{	1 / (36 * DetJ) * (-4 * Lx * Lx * Ly * Ly + 2 * (Lx * Lx + Ly * Ly) * Lz * Lz),
																	1 / (36 * DetJ) * (-2 * Ly * Ly * Lz * Lz + Lx * Lx * (-2 * Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * 4 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly + Lz + Lz)),
																	1 / (36 * DetJ) * (Ly * Ly * Lz * Lz - 2 * Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (-4 * Ly * Ly * Lz * Lz + 2 * Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * 2 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly - 2 * Lz * Lz)),
																	1 / (36 * DetJ) * (- Ly * Ly * Lz * Lz - Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (-4 * Ly * Ly * Lz * Lz + 2 * Lx * Lx * (Ly * Ly + Lz * Lz))},
																{	1 / (36 * DetJ) * 2 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly - 2 * Lz * Lz)),
																	1 / (36 * DetJ) * (Lx * Lx * Ly * Ly - 2 * (Lx * Lx + Ly * Ly) * Lz * Lz),
																	1 / (36 * DetJ) * (Ly * Ly * Lz * Lz - 2 * Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * 4 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly + Lz + Lz)),
																	1 / (36 * DetJ) * (- Ly * Ly * Lz * Lz - Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (-4 * Lx * Lx * Ly * Ly + 2 * (Lx * Lx + Ly * Ly) * Lz * Lz),
																	1 / (36 * DetJ) * (-4 * Ly * Ly * Lz * Lz + 2 * Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (- Ly * Ly * Lz * Lz - Lx * Lx * (Ly * Ly + Lz * Lz))},
																{	1 / (36 * DetJ) * (-2 * Ly * Ly * Lz * Lz + Lx * Lx * (-2 * Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (-4 * Lx * Lx * Ly * Ly + 2 * (Lx * Lx + Ly * Ly) * Lz * Lz),
																	1 / (36 * DetJ) * (-4 * Ly * Ly * Lz * Lz + 2 * Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (- Ly * Ly * Lz * Lz - Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * 4 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly + Lz + Lz)),
																	1 / (36 * DetJ) * (Lx * Lx * Ly * Ly - 2 * (Lx * Lx + Ly * Ly) * Lz * Lz),
																	1 / (36 * DetJ) * (Ly * Ly * Lz * Lz - 2 * Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * 4 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly + Lz + Lz))},
																{	1 / (36 * DetJ) * (Ly * Ly * Lz * Lz - 2 * Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (- Ly * Ly * Lz * Lz - Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * 2 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly - 2 * Lz * Lz)),
																	1 / (36 * DetJ) * (-4 * Lx * Lx * Ly * Ly + 2 * (Lx * Lx + Ly * Ly) * Lz * Lz),
																	1 / (36 * DetJ) * (Lx * Lx * Ly * Ly - 2 * (Lx * Lx + Ly * Ly) * Lz * Lz),
																	1 / (36 * DetJ) * 4 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly + Lz + Lz)),
																	1 / (36 * DetJ) * (-2 * Ly * Ly * Lz * Lz + Lx * Lx * (-2 * Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (Lx * Lx * Ly * Ly - 2 * (Lx * Lx + Ly * Ly) * Lz * Lz)},
																{	1 / (36 * DetJ) * (Lx * Lx * Ly * Ly - 2 * (Lx * Lx + Ly * Ly) * Lz * Lz),
																	1 / (36 * DetJ) * 2 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly - 2 * Lz * Lz)),
																	1 / (36 * DetJ) * (- Ly * Ly * Lz * Lz - Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (-4 * Ly * Ly * Lz * Lz + 2 * Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (Ly * Ly * Lz * Lz - 2 * Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (-2 * Ly * Ly * Lz * Lz + Lx * Lx * (-2 * Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * 4 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly + Lz + Lz)),
																	1 / (36 * DetJ) * (Ly * Ly * Lz * Lz - 2 * Lx * Lx * (Ly * Ly + Lz * Lz))},
																{	1 / (36 * DetJ) * (-2 * Ly * Ly * Lz * Lz + Lx * Lx * (-2 * Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (-4 * Ly * Ly * Lz * Lz + 2 * Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (-4 * Ly * Ly * Lz * Lz + 2 * Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * (- Ly * Ly * Lz * Lz - Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * 4 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly + Lz + Lz)),
																	1 / (36 * DetJ) * (Lx * Lx * Ly * Ly - 2 * (Lx * Lx + Ly * Ly) * Lz * Lz),
																	1 / (36 * DetJ) * (Ly * Ly * Lz * Lz - 2 * Lx * Lx * (Ly * Ly + Lz * Lz)),
																	1 / (36 * DetJ) * 4 * (Ly * Ly * Lz * Lz + Lx * Lx * (Ly * Ly + Lz + Lz))}
															};
	const std::vector<double> LumpedMatrix {1 / 8, 1 / 8,  1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8};

private:

	const double Lx;
	const double Ly;
	const double Lz;
	const double DetJ = Lx * Ly * Lz;

};

#endif