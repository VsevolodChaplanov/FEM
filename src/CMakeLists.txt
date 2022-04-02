set(SRC main.cpp)

add_executable(${TEST_EXE} ${SRC})
target_include_directories(${TEST_EXE} PRIVATE lib/include)

add_subdirectory(lib)

target_link_libraries(${TEST_EXE} PRIVATE ${FEM_PDE_lib})