set(SRC main.cpp)

add_subdirectory(lib)

add_executable(${TEST_EXE} ${SRC})

target_link_libraries(${TEST_EXE} PRIVATE ${FEM_PDE_lib})