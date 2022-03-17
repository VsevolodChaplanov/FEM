set(SRC main.cpp)

add_executable(${TEST_EXE} ${SRC})

target_link_libraries(${TEST_EXE} PRIVATE FEM_lib HEADERS)