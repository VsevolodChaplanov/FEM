#ifndef __VECTOROPERATIONS__
#define __VECTOROPERATIONS__

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "CompressedM.h"

// Процедуры основных векторных и матричных операций

double norm_2(const std::vector<double> &vec);

double max_abs(const std::vector<double> &vec);

// Скалярное произведение векторов
double dot_product(const std::vector<double> &a, const std::vector<double> &b);

// Умножение вектора на число
std::vector<double> mult_n(const std::vector<double> &a, const double &b);

// Сложение векторов
std::vector<double> vector_sum(const std::vector <double> &a, const std::vector <double> &b);

// Вычитание векторов
std::vector<double> vector_diff(const std::vector<double> &a, const std::vector<double> &b);

// Запись вектора a в файл с названием Filename.csv
void write_in_file(const std::vector<double> &a, const std::string &Filename);

bool check_matrix_sym(CMatrix &Matrix);

void summ_cm(const CMatrix &A, const CMatrix &B, CMatrix &Lhs);

#endif