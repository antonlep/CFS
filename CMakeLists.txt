cmake_minimum_required(VERSION 3.18)

project(cfs LANGUAGES CXX)

option(CFS_BUILD_TESTS "Build CFS tests" OFF)

include(CTest)
enable_testing()

add_subdirectory(src/solver)

