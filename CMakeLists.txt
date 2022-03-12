cmake_minimum_required(VERSION 3.0.0)
project(FEM VERSION 0.1.0)

set(TEST_EXE fem)

message(STATUS ${TEST_EXE})


#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#set(BOOST_DIR "/home/vsevolod/Programs/BOOST/boost_1_61_0/boost")
#set(Boost_NO_SYSTEM_PATHS ON)
#find_package(Boost COMPONENTS geometry REQUIRED)
#find_package(Threads REQUIRED)
#find_package(OpenMP REQUIRED)


add_subdirectory("src/")
