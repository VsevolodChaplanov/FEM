// #ifndef __VECTOROPERATIONS_CPP__
// #define __VECTOROPERATIONS_CPP__

#include "VectorOperations.h"

// Процедуры основных векторных и матричных операций

double norm_2(const std::vector<double> &vec)
{
	double n2 = 0;
	for (auto elem : vec)
	{
		n2 += elem * elem;
	}
	return sqrt(n2 / (double) vec.size());
}

double max_abs(const std::vector<double> &vec)
{
	double max_elem = fabs(vec[0]);
	for (auto elem : vec)
	{
		if (fabs(elem) > max_elem)
		{
			max_elem = fabs(elem);
		}
		
	}
	return max_elem;
}

// Скалярное произведение векторов
double dot_product(const std::vector<double> &a, const std::vector<double> &b)
{
    double Result = 0;

    for (size_t i = 0; i < a.size(); i++)
    {
        Result += a[i] * b[i];
    } // i

    return Result;
}

// Умножение вектора на число
std::vector<double> mult_n(const std::vector<double> &a, const double &b)
{
    std::vector<double> Result(a.size(), 0.);
    for (size_t i = 0; i < a.size(); i++)
    {
        Result[i] = a[i] * b;
    } // i
    return Result;
}

// Сложение векторов
std::vector<double> vector_sum(const std::vector <double> &a, const std::vector <double> &b)
{
    std::vector<double> Result(a.size(), 0.);
    for (size_t i = 0; i < a.size(); i++)
    {
        Result[i] = a[i] + b[i];
    }
    return Result;
}

// Вычитание векторов
std::vector<double> vector_diff(const std::vector<double> &a, const std::vector<double> &b)
{
    std::vector<double> Result(a.size(), 0.);
    for (size_t i = 0; i < a.size(); i++)
    {
        Result[i] = a[i] - b[i];
    } // i
    return Result;
}

// Запись вектора a в файл с названием Filename.csv
void write_in_file(const std::vector<double> &a, const std::string &Filename)
{
	std::ofstream file;
    file.open(Filename + ".csv");
    if (file.is_open())
    {
		for (auto elem : a)
		{
			file << elem << std::endl;
		} // elem
    } //endif
    file.close();
}

bool check_matrix_sym(CMatrix &Matrix)
{
    bool result = true;
    for (size_t i = 0; i < Matrix.size(); i++)
    {
        for (auto elem : Matrix[i])
        {
            if (elem.second != Matrix.GetValue(elem.first, i))
            {
                result = false;
            }
        }
    }
    return result;
}

void summ_cm(const CMatrix &A, const CMatrix &B, CMatrix &Lhs)
{
	std::size_t N = A.size();
	for (size_t i = 0; i < N; i++)
	{
		for (auto elem : A[i])
		{
			Lhs.SetValue(i, elem.first, Lhs.GetValue(i, elem.first) + elem.second);
		}
		for (auto elem : B[i])
		{
			Lhs.SetValue(i, elem.first, Lhs.GetValue(i, elem.first) + elem.second);
		}
	}
}

// #endif