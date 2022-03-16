#ifndef __LINELEM__
#define __LINELEM__

#include <vector>
#include "IFiniteElem.h"

// Строит элемент на двух вершинах
// Хранит в себе
// 			* Базисные фукнции этих вершин
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
class LinElem : public IFiniteElement
{
public:

	// Принимает (массив из координат двух вершин [0.9, 1],  глобальные индексы этих вершин [9, 10])
	LinElem(const std::vector<double> &vertices, const std::vector<std::size_t> &GIndices); 						// Конструктор 
																													// По 2-ум вершинам строит определяет линейный элемент
	std::vector<std::vector<double>>* GetMass(const std::size_t , const std::size_t ) override;														// Возвращает элемент локальной матрицы масс стояший на позиции [i,j]
	std::vector<std::vector<double>>* GetStiffness(const std::size_t , const std::size_t ) override;														// Возвращает элемент локальной матрицы жексткости стояший на позиции [i,j]
	std::vector<double>* GetLumped(const std::size_t i, const std::size_t j) override;																			// Возвращает элемент lumped mass matrix стояший на позиции [i,j]
	std::vector<double> PhysToParam(const double x, const double y = 0, const double z = 0) override;																				// Переводит точку в физическом пространстве, в точку в параметрическом пространстве
	std::vector<double> PhysToParam(const std::vector<double> &point) override;

private: // Methods

	
	double phi1(const double &_xi);	// Mathes to the left vertex																			// Базисная функция линейного элемента phi(x) = 1 -x
	double phi2(const double &_xi);	// Mathes to the right vertex														// Базнсная фукнция линейного элемента phi(x) = x

private: // Properties

	double Length; 																									// Длинна элемента
	double StartPoint;		                                  														// Точка начала отрезка

public:

	std::vector<std::size_t> GlobalIndices;	
	std::size_t NBasisFunction = GlobalIndices.size(); // Тут поменять когда будет массив базисных функций																		// Массив глобальных идексов объекта
	std::vector<std::vector<double>> MassMatrix {{Length / 3, Length / 6}, {Length / 6, Length / 3}};			// Локальная матрица масс
	std::vector<std::vector<double>> StiffnessMatrix {{Length / 3, Length / 6}, {Length / 6, Length / 3}}; 	// Локальная матрица Жесткости
	std::vector<double> LumpedMassMatrix {Length / 2, Length / 2}; 	// Локальная lumped mass матрица
};

#endif