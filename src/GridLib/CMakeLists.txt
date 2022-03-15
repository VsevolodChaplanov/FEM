set(SRC IFiniteElem.cpp
	LinElem.cpp
	RectElem.cpp
	CubeElem.cpp
	MeshBuilder.cpp)

set(HEADERS MeshBuilder.h
	IFiniteElem.h
	CubeElem.h
	LinElem.h
	RectElem.h
	IFiniteElem.h)

add_library(FINELEM STATIC ${SRC} ${HEADERS})