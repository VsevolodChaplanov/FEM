include_directories(GridLib)

set(SRC main.cpp)

add_executable(${TEST_EXE} ${SRC} )		  
add_subdirectory(GridLib)

target_link_libraries(${TEST_EXE} GridLib)

