#ifndef __COMPRESSEDFORMAT__
#define __COMPRESSEDFORMAT__

#include <vector>
#include <map>
#include <string>
#include <random>
#include <iostream>

// Класс для харнения разряженных матриц в виде:
// { {{Coloumn, Value}, {Coloumn, Value}, ...}, 
//	 {{Coloumn, Value}, ... },
// 	 {{Coloumn, Value}, {Coloumn, Value}, ...} }
class CMatrix
{
private:
	
	// Число строк в матрице 
	int N;

	// Формал хранения отличных от 0 элементов матрицы в виде:
	// { {{Coloumn, Value}, {Coloumn, Value}, ...}, 
	//	 {{Coloumn, Value}, ... },
	// 	 {{Coloumn, Value}, {Coloumn, Value}, ...} }
	// Первым элементом в map хранится индекс столбца в котором находится значение Value
	std::vector<std::map<size_t, double>> CompressedMatrix;

public:

	CMatrix();

	// Конструктор создания пустой матрица в N строк
	CMatrix(size_t N);

	// Конструктор создания заполненной матрицы случайными числами
	// 3 Случайных числа в строку + добавление чисел для симметрии -> может быть больше 7 чисел
	CMatrix(size_t N, const std::string &Random);

	// Получить элемент на месте [i,j]
	double GetValue(size_t i, size_t j) const;

	// Добавить элемент (а) в матрицу на место [i,j] 
	void SetValue(size_t i, size_t j, double a);

	// Обнулить строку матрицы
	void SetZeroRow(size_t i);

	// Умонжение разреженной матрицы на вектор
	std::vector<double> operator*(const std::vector<double> &Vector) const;

	// Умножение матрицы на число
	CMatrix operator*(double a) const;

	// Умножение матриц разряженного типа
	CMatrix operator*(CMatrix &Matrix) const;

	// Перегрузка скобок, т.к. часто необходимо обращаться к строкам матрицы
	// Так будет быстрее чем GetValue и удобней
	std::map<size_t, double> operator[](size_t i) const;

	// Перегрузка оператора = 
	// Используется для конструкций CMatrix C = A * B
	CMatrix operator=(CMatrix &Matrix_A);

	// Возвращает число строк матрицы
	size_t size() const;

	void resize(size_t N);

	// Деструктор по умолчанию
	~CMatrix();
};

#endif