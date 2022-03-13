set(FIN_ELEMS_LIB 	Builder.cpp
				 	LinElem.cpp
				  	RectElem.cpp
					CubeElem.cpp)

add_library(FINELEM STATIC ${FIN_ELEMS_LIB})