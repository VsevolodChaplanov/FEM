#ifndef __LINELEM__
#define __LINELEM__

#include <vector>

// Строит элемент на двух вершинах
// Хранит в себе
// 			* Базисные фукнции этиъ вершин
// 			* Глобальные индексы верших этого элемента
//			* Собирает в себе
//				** локлаьную матрицу жесткости
//				** локальную матрицу масс
//				** матрицу вз. масс
//
// Отрезок*
//
//			a ------------ b
//
// Хранится в параметроической плоскости как отрезок:
//
// 			0 ------------ 1
//
//
class LinElem
{
public:

	// Принимает (массив из координат двух вершин [0.9, 1],  глобальные индексы этих вершин [9, 10])
	LinElem(const std::vector<double> &vertices, const std::vector<std::size_t> &GIndices); 						// Конструктор 
																													// По 2-ум вершинам строит определяет линейный элемент

	double GetMass(const std::size_t &, const std::size_t &);														// Возвращает элемент локальной матрицы масс стояший на позиции [i,j]
	double GetStiff(const std::size_t &, const std::size_t &);														// Возвращает элемент локальной матрицы жексткости стояший на позиции [i,j]
	double GetLumped(const std::size_t &);																			// Возвращает элемент lumped mass matrix стояший на позиции [i,j]
	double P_to_Param(const double &);																				// Переводит точку в физическом пространстве, в точку в параметрическом пространстве

private: // Methods

	double phi1(const double &_xi);	// Mathes to the left vertex																			// Базисная функция линейного элемента phi(x) = 1 -x
	double phi2(const double &_xi);	// Mathes to the right vertex														// Базнсная фукнция линейного элемента phi(x) = x

private: // Properties

	double Length; 																									// Длинна элемента
	double StartPoint;		                                  														// Точка начала отрезка

public:

	std::vector<std::size_t> GlobalIndices;	
	const std::size_t NBasisFunction = GlobalIndices.size(); // Тут поменять когда будет массив базисных функций																		// Массив глобальных идексов объекта
	const std::vector<std::vector<double>> MassMatrix {{Length / 3, Length / 6}, {Length / 6, Length / 3}};			// Локальная матрица масс
	const std::vector<std::vector<double>> StiffnessMatrix {{Length / 3, Length / 6}, {Length / 6, Length / 3}}; 	// Локальная матрица Жесткости
	const std::vector<double> LumpedMassMatrix {Length / 2, Length / 2}; 	// Локальная lumped mass матрица
};

// TODO
// 		* Можно изменить динамикческие массивы веткоров на статические
//				т.к. они имеют всегда заданный размер, должно быть быстрее и эффективнее?
//		* В интерфейсе написана обертка для базиса 

#endif