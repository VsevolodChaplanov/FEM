#ifndef __GLOBAL_MATRIX_ASSEMBLERS__
#define __GLOBAL_MATRIX_ASSEMBLERS__

#include "CompressedM.h"
#include "FemGrid.h"

class GLobalMatrixAssembler
{
private:

	CMatrix Matrix;

public:

	GLobalMatrixAssembler(size_t N);
	void add_local_matrix(const std::vector<size_t> &global_indices, const std::vector<double> &local_matrix, double coeff = 1);
	CMatrix get_result() const;
};

#endif