cmake_minimum_required(VERSION 3.5)

Include(FetchContent)

FetchContent_Declare(
	Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG        v3.0.0-preview3
)

FetchContent_MakeAvailable(Catch2)

set(SRC catch-tests.cpp
	main-tests.cpp)
set(HEADERS catch-tests.h)

set(LIB_TESTS fem_pde_lib-test)

add_executable(${LIB_TESTS} ${SRC} ${HEADERS})

target_link_libraries(${LIB_TESTS} PRIVATE ${FEM_PDE_lib})
target_link_libraries(${LIB_TESTS} PRIVATE Catch2::Catch2)
target_include_directories(${LIB_TESTS} PRIVATE ${CMAKE_SOURCE_DIR}/src/lib/include)
link_directories(${catch2_SOURCE_DIR})

include(CTest)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
