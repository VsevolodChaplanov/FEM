#include "GlobalAssemblers.h"

GLobalMatrixAssembler::GLobalMatrixAssembler(size_t N)
{
	Matrix.resize(N);
}

CMatrix GLobalMatrixAssembler::get_result() const
{
	return Matrix;
}

void GLobalMatrixAssembler::add_local_matrix(const std::vector<size_t> &global_indices, const std::vector<double> &local_matrix, double coeff)
{
	size_t N = global_indices.size();
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			size_t GlobI = global_indices[i];
			size_t GlobJ = global_indices[j];
			Matrix.SetValue(GlobI, GlobJ, Matrix.GetValue(GlobI, GlobJ) + coeff * local_matrix[i * N + j]);
		}
	}
}

