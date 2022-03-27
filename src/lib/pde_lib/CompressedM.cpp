// #ifndef __COMPRESSEDFORMAT_CPP__
// #define __COMPRESSEDFORMAT_CPP__

#include "CompressedM.h"

// ============= C++ ============= //
// ========= Реализации ========== //

CMatrix::CMatrix() { }

CMatrix::CMatrix(size_t N)
{
	this->N = N;
	CompressedMatrix.resize(N);
}

CMatrix::CMatrix(size_t N, const std::string &Random)
{
	this->N = N;
	std::random_device rd;  // Не детерминированный генератор
							// Не высокая производительность генератора
							// Высокое качество генерируемых случайных чисел
	std::mt19937 gen(rd()); // 32-разрядный механизм типа "Вихрь Мерсенна"

	// Генератор равномерного распределения чисел от a до b 	// ++
	std::uniform_int_distribution<> Number(-5,-1); 				// Получить распределение от -5 до -1

	CompressedMatrix.resize(N);

	for (size_t i = 1; i < N; i++)
	{	
		// До i-1 чтобы не генерировались числа на диагональ
		std::uniform_int_distribution<> Coloumn(0,  i - 1); // Получить распределение от -1 до 1
		// Здесь можно записать всё в несколько строк, но таким образом удобнее пользоваться отладкой
		int first = Coloumn(gen);
		int first_number = Number(gen);
		//int second = Coloumn(gen);
		//int second_number = Number(gen);
		//int third = Coloumn(gen);
		//int third_number = Number(gen);
		CompressedMatrix[i].insert(std::make_pair(first, first_number));
		//CompressedMatrix[i].insert(std::make_pair(second, second_number));
		//CompressedMatrix[i].insert(std::make_pair(third, third_number));

		CompressedMatrix[first].insert(std::make_pair(i, first_number));
		//CompressedMatrix[second].insert(std::make_pair(i, second_number));
		//CompressedMatrix[third].insert(std::make_pair(i, third_number));
	} // i

	// На "диагональ" поставт сумму всех элементов в строке + 1 для диагонального преобладания
	for (size_t i = 0; i < N; i++)
	{
		double Sum = 0;
		for (auto elem : CompressedMatrix[i])
		{
			Sum += elem.second;
		} // elem
		CompressedMatrix[i].insert(std::make_pair(i, fabs(Sum) + 1));
	} // i
}

// Получить элемент на месте [i,j]
double CMatrix::GetValue(size_t i, size_t j) const 
{
	// Такой алгоритм быстрее чем 
	// return CompressedMatrix[i][j]; Не понятно почему
	auto fnd = CompressedMatrix[i].find(j);
	if (fnd == CompressedMatrix[i].end())
	{
		return 0;
	} 
	else
	{
		return fnd->second;
	}
}

// Добавить элемент (а) в матрицу на место [i,j] 
void CMatrix::SetValue(size_t i, size_t j, double a)
{
	CompressedMatrix[i][j] = a;
}

// Умонжение разреженной матрицы на вектор
std::vector<double> CMatrix::operator*(const std::vector<double> &Vector) const
{
	if (Vector.size() == ! N)
	{
		throw std::runtime_error("Dimensions of vector and matrix doesn't match");
	}
	
	// Вектор в который записывается результат
	std::vector<double> Result(N);
	for (size_t i = 0; i < N; i++)
	{	
		// Внутренняя сумма (*Скалярное произведение строки на вектор*)
		double Sum = 0;
		for (auto elem : CompressedMatrix[i])
		{
			Sum += Vector[elem.first] * elem.second;
		} // elem  	// Наиболее быстрый найденный вариант
					// Хоть и в геттере есть возвращение 0 на любой элемент матрицы если его нет в @CompressedMatrix
					// Пробег лишь по не нулевым элементам быстрее
		Result[i] = Sum;
	} // i
	return Result;
}

// Умножение матрицы на число
CMatrix CMatrix::operator*(const double a) const
{
	// Возвращается новая матрица
	CMatrix Result(N);
	for (size_t i = 0; i < N; i++)
	{
		for (auto elem : CompressedMatrix[i])
		{
			Result.SetValue(i, elem.first, elem.second * a);
		} // elem 	// Наискорейший вариант
	} // i
	return Result;
}

// Умножение матриц разряженного типа
CMatrix CMatrix::operator*(CMatrix &Matrix) const
{
	// Возвращается новая матрица
	CMatrix Result(N);
	for (size_t i = 0; i < N; i++) 
	{			
		for (size_t j = 0; j < N; j++) 
		{
			// Внутренняя сумма
			double Sum = 0;
			for (auto elem : CompressedMatrix[i]) 
			{
				// ! Тут поменял
				Sum += elem.second * Matrix.GetValue(elem.first, j);
			} // elem 
			if (Sum != 0)
			{
				Result.SetValue(i,j, Sum);
			} //endif
		} // j
	} // i
	return Result;
}

// Перегрузка скобок, т.к. часто необходимо обращаться к строкам матрицы
// Так будет быстрее чем GetValue и удобней
std::map<size_t, double> CMatrix::operator[](size_t i) const
{
	return CompressedMatrix[i];
}

// Перегрузка оператора = 
// Используется для конструкций CMatrix C = A * B
CMatrix CMatrix::operator=(CMatrix &Matrix_A)
{
	for (size_t i = 0; i < Matrix_A.size(); i++)
	{
		this->CompressedMatrix[i] = Matrix_A[i];
	}
	return *this;
}
// CMatrix CMatrix::operator=(const CMatrix &Matrix_A)
// {
// 	// Результат новая матрица B
// 	CMatrix B(Matrix_A.size());

// 	for (size_t i = 0; i < Matrix_A.size(); i++)
// 	{
// 		B[i] = Matrix_A[i];
// 	} // i
	
// 	return B;
// }


// Возвращает число строк матрицы
size_t CMatrix::size() const
{
	return this->N;
}

void CMatrix::resize(size_t N)
{
	CompressedMatrix.resize(N);
	this->N = N;
}

void CMatrix::SetZeroRow(size_t i)
{
	CompressedMatrix[i].clear();
}

CMatrix::~CMatrix() { }

// #endif