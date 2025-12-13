cmake_minimum_required(VERSION 3.18)
project(solver)

find_package(pybind11 CONFIG REQUIRED)
add_library(Eigen3::Eigen INTERFACE IMPORTED)
target_include_directories(Eigen3::Eigen INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen
)

pybind11_add_module(solver
    bindings.cpp
    fem.cpp
    mesh.cpp
)

target_link_libraries(solver PRIVATE Eigen3::Eigen)

install(TARGETS solver DESTINATION .)