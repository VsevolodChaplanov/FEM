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
	std::vector<std::map<int , double>> CompressedMatrix;

public:

	CMatrix();

	// Конструктор создания пустой матрица в N строк
	CMatrix(const int &N);

	// Конструктор создания заполненной матрицы случайными числами
	// 3 Случайных числа в строку + добавление чисел для симметрии -> может быть больше 7 чисел
	CMatrix(const int &N, const std::string &Random);

	// Получить элемент на месте [i,j]
	double GetValue(const int i, const int j) const;

	// Добавить элемент (а) в матрицу на место [i,j] 
	void SetValue(int i, int j, double a);

	// Обнулить строку матрицы
	void SetZeroRow(int i);

	// Умонжение разреженной матрицы на вектор
	std::vector<double> operator*(const std::vector<double> &Vector);

	// Умножение матрицы на число
	CMatrix operator*(const double &a);

	// Умножение матриц разряженного типа
	CMatrix operator*(CMatrix &Matrix);

	// Перегрузка скобок, т.к. часто необходимо обращаться к строкам матрицы
	// Так будет быстрее чем GetValue и удобней
	std::map<int, double> operator[](int i);

	// Перегрузка оператора = 
	// Используется для конструкций CMatrix C = A * B
	CMatrix operator=(CMatrix &Matrix_A);

	// Возвращает число строк матрицы
	int size() const;

	void resize(const std::size_t& N);

	// Деструктор по умолчанию
	~CMatrix();
};

#endif