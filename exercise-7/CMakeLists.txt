cmake_minimum_required(VERSION 3.10)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(Exercise7)

add_executable(exercise-7 ./main.c)
target_link_libraries(exercise-7 m)
