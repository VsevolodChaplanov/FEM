cmake_minimum_required(VERSION 3.5)

Include(FetchContent)

FetchContent_Declare(
	Catch2
	GIT_REPOSITORY https://github.com/catchorg/Catch2.git
	GIT_TAG        v3.0.0-preview3
)

FetchContent_MakeAvailable(Catch2)

add_executable(catch_tests catch-tests.cpp)
target_link_libraries(catch_tests PRIVATE ${FEM_PDE_lib})
target_link_libraries(catch_tests PRIVATE Catch2::Catch2WithMain)

set_target_properties(catch_tests
	PROPERTIES 
		RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests
)

include(CTest)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
