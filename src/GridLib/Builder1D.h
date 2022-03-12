#ifndef __BUILDER_1D__
#define __BUILDER_1D__

#include "LinElem.cpp"

class Builder1D
{
public: // Properties

private:

	// Хранение элементов
	std::vector<LinElem*> Elements;

public: // Methods

	// Construct elements on a uniform grid
	Builder1D(const double &, const double &, const std::size_t &);
	// Construct elements on a unstructured grid
	Builder1D(const std::vector<double> &Mesh);

private:

};



#endif