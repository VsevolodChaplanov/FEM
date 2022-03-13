#ifndef __BUILDER_1D__
#define __BUILDER_1D__

#include "LinElem.cpp"

#include <cmath>
#include <vector>

class Builder1D
{
public: // Properties

	// Хранение элементов
	std::vector<LinElem*> elems;	

private:

	// Число элементов
	const std::size_t nn;

public: // Methods

	// Construct elements on a uniform grid
	Builder1D(const double &, const double &, const std::size_t &);
	// Construct elements on a unstructured grid
	Builder1D(const std::vector<double> &Mesh);
	// Get number of elements
	std::size_t GetNElem(void);

private:

};



#endif