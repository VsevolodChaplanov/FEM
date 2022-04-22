## Установка
В командной строке необходимо выполнить следующие команды:
``` sh
mkdir build
cd build
cmake ..
make
```
## Запуск тестов
Для запуска тестов необходимо выполнить следующие команды:
``` sh
cd bin
./fem_pde_lib-test
``` 

## Использование
### Конечно-элементная сетка
Формат хранения сеток *.vtk*, с использованием этого же формата осущетсвляется сохранение решения

Считываение конечно-элементной сетки происходит по команде:
``` c++
FemGrid YourGrid = Builder::BuildFromFile("NameoftheFile.vtk")
```
### Задание параметров решателя
Параметры решателя задаются следующией командой:
``` c++
MatrixSolverParams* YourParams = new MatrixSolverParams(
	Method,
	PreconditionMethod,
	MaximumNumberOfIterations,
	Accurancy,
	ResidualsSaveInterval,
	RelaxationParamforSOR/SSOR,
	RelaxationParamforSSORPreconditioner
);
```
Дефолтные параметры установленные для решателя:
``` с++
MatrixSolverParams* YourParams = new MatrixSolverParams(
	Method = MatrixSolverParams::Methods::Thomas
	PreconditionMethod
	MatrixSolverParams::Preconditioners::None, 
	MaximumNumberOfIteration = 1000,
	Accurancy = 1.e-5, 
	ResidualsSaveInterval = 10,
	RelaxationParamforSOR/SSOR = 1.95,
	RelaxationParamforSSORPreconditioner = 1.95
);
```
### Создание объекта ДУ
Решается уравнение типа:
$\nabla({k \nabla u}) + u = f$
Объект ДУ создается следующим образом:
``` c++
FemPDE YourPDE(&YourGrid, function f, function k, YourParams);
```
Для задания граничных условий Дирихле используется селектирующая функция:
``` c++
YourPDE.apply_boundary_condition_dirichlet(value, YourGrid.boundary_element_indices(function YourSelectorFunction));
```

### Получение решения
Получить решение и сохранить его можно следующим образом:
``` c++
std::vector<double> solution = YourPDE.solve();
YourGrid.savevtk(solution, "TargetFileName.vtk");
```

#### Пример можно найти в файле *main.cpp*

