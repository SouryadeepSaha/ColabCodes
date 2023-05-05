cmake_minimum_required(VERSION 3.25)
set( CMAKE_CXX_COMPILER "C:/MinGW/bin/g++.exe" )
set( CMAKE_C_COMPILER "C:/MinGW/bin/gcc.exe" )
set( CMAKE_GENERATOR "MinGW Makefiles")

project(ClionProject)

set(CMAKE_CXX_STANDARD 20)

add_executable(${PROJECT_NAME} main.cpp dependencies.h)
add_subdirectory(source)

target_link_libraries(${PROJECT_NAME} Source )
