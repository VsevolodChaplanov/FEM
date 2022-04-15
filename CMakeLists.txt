cmake_minimum_required(VERSION 3.0.0)

project(FEM VERSION 0.1.0)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory("src/")

# add_custom_target(test_resources)
# add_custom_command(
#         TARGET test_resources POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E copy
#                 ${CMAKE_CURRENT_BINARY_DIR}/test_resources
#                 ${CMAKE_SOURCE_DIR})

file(INSTALL ./test_resources DESTINATION ${CMAKE_BINARY_DIR}/bin/test_resources)