#ifndef __INTERFACE_FINITE_ELEMS__
#define __INTERFACE_FINITE_ELEMS__

#include <cstddef>
#include <vector>
#include <functional>
#include <string>

class Element
{
public:

private:

	std::vector<std::function<double(const double&)>> BasisFunctions;		// Массив базисных фукнций

	std::size_t Nbasis; 													// Number of basis functions

	std::vector<std::size_t> GIndices; 										// Global indices of local vertices

};

#endif
