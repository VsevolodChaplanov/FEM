cmake_minimum_required(VERSION 3.0.0)
project(FEM VERSION 0.1.0)

set(TEST_EXE fem)

message(STATUS ${TEST_EXE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory("src/")
