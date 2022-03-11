#include <iostream>
#include "boost/geometry/geometry.hpp"
#include "GridLib/LinElem.h"

int main(int argc, char const *argv[])
{
	LinElem* MyfirstFE = new LinElem({0,1}, {0,1});
	std::cout << "В матрице масс на месте 0,0 стоит: " << MyfirstFE->GetMass(0,0) << std::endl;
	return 0;
}
