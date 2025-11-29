cmake_minimum_required(VERSION 3.18)
project(solver)

find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)

pybind11_add_module(solver
    bindings.cpp
    fem.cpp
    mesh.cpp
)

target_link_libraries(solver PRIVATE Eigen3::Eigen)