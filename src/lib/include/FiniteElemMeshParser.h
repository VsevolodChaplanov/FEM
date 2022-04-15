#ifndef __FINITE_ELEMENTS_MESH_PARSER__
#define __FINITE_ELEMENTS_MESH_PARSER__

#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <boost/token_iterator.hpp>
#include <boost/algorithm/algorithm.hpp>

// Parses finite elements mesh generated by gmsh
class IFEMParser
{
protected:

	std::string filename;
	// Continous storage of the vertices as:
	// 	{x0,y0,z0,x1,y1,z1, ... }
	size_t Nvert;
	size_t Nelem;
	std::vector<double> vertices;
	std::vector<std::vector<size_t>> cells;
	std::vector<size_t> cell_types;

	enum DataType {
		Null,
		POINTS,
		CELLS,
		CELL_TYPES
	};

	DataType state = DataType::Null;

public:

	IFEMParser(const std::string &filename);
	static IFEMParser* Factory(const std::string &filename);
	virtual void load_mesh() = 0;
	const std::vector<double>& get_vertices() const;
	const std::vector<size_t>& get_cell_types() const;
	const std::vector<std::vector<size_t>>& get_cells() const;
	size_t get_vertices_number() const;
	size_t get_elements_number() const;
	virtual ~IFEMParser();
};

// Parses finite elements mesh generated by gmsh in vtk format
class VtkFEMParser : public IFEMParser
{
private:

public:

	VtkFEMParser(const std::string &filename);
	void load_mesh() override;
	~VtkFEMParser();
};


#endif