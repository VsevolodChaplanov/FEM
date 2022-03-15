#ifndef __MESHBUILDER__
#define __MESHBUILDER__

#include <cmath>
#include <vector>
#include "IFiniteElem.h"

// Used Builder pattern


// Abstract builder
class IBuilder
{
public:

	std::vector<double> _temp_vertex_;

protected:

	std::vector<double> bounds;
	std::vector<std::size_t> partition_num;
	std::size_t elnum;

public:

	std::vector<double>* GetTempVerts() {return &_temp_vertex_;}
	void SetElementsNumber(const std::size_t NumOfElems) { this->elnum = NumOfElems; }
	virtual void SetElementType(const std::size_t);
	virtual void buildElement() = 0;
	virtual IFiniteElement* CreateElement(const std::vector<double> &localverts, const std::vector<std::size_t> &GIndices) = 0;
	virtual void SetBounds(const std::vector<double> &) = 0;
	virtual void SetPartition(const std::vector<std::size_t> &) = 0;
	virtual void MakeTempVerts() = 0;

public:

};


// Строитель сетки из одномерных элементов
class BuilderLinear : public IBuilder
{
private:

	double hx;
	std::size_t ElemType;

public:

	BuilderLinear();
	void SetElementType(const std::size_t vtx_N) override;
	void buildElement() override;
	IFiniteElement* CreateElement(const std::vector<double> &localverts, const std::vector<std::size_t> &GIndices) override;
	void SetBounds(const std::vector<double> &) override;
	void SetPartition(const std::vector<std::size_t> &) override;
	void MakeTempVerts();

protected:
};

class BuilderRect : public IBuilder
{
private:

	double hx;
	double hy;

public:

	BuilderRect();
	void buildElement() override;
	IFiniteElement* CreateElement(const std::vector<double> &localverts, const std::vector<std::size_t> &GIndices) override;
	void SetBounds(const std::vector<double> &) override;
	void SetPartition(const std::vector<std::size_t> &) override;

private:
};

class BuilderCube : public IBuilder
{
private:

	double hx;
	double hy;
	double hz;

public:

	BuilderCube();
	void buildElement() override;
	virtual void CreateElement(IFiniteElement* element);
	void SetBounds(const std::vector<double> &);
	void SetPartition(const std::vector<std::size_t> &);
};

// Director
class ElemsMesh
{
protected:

	std::vector<IFiniteElement*> mesh;
	IBuilder* builder;
	std::vector<double> vertices;
	std::vector<std::size_t> partition_nums;
	std::size_t vert_num = 1;
public:

public:

	// vertices - Вершины определяют границы области, partition_nums - определяют число элементов вдоль осей
	// 
	ElemsMesh(const std::vector<double> &vertices, const std::vector<std::size_t> &partition_nums, IBuilder* builder); 
	void SetElemBuilder(IBuilder* builder);
	void SetVertNum() 
	{
		for (auto elem : partition_nums)
		{
			vert_num *= (elem + 1);
		}
	}
	void ConstructMesh();
};

ElemsMesh::ElemsMesh(const std::vector<double> &vertices, const std::vector<std::size_t> &partition_nums, IBuilder* builder) : builder(builder) , vertices(vertices), partition_nums(partition_nums) 
{  
	SetVertNum();
}

void ElemsMesh::SetElemBuilder(IBuilder* builder)
{
	this->builder = builder;
}

void ElemsMesh::ConstructMesh()
{
	builder->SetPartition(partition_nums);
	builder->SetBounds(vertices);
	builder->MakeTempVerts();
	for (size_t i = 0; i < vert_num; i++)
	{
		mesh.push_back(builder->CreateElement( {builder->_temp_vertex_[i], builder->_temp_vertex_[i+1]}, {i, i+1}));
	}
}

#endif