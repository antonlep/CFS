cmake_minimum_required(VERSION 3.18)
project(solver LANGUAGES CXX)

find_package(pybind11 REQUIRED)

add_library(solver_core
    core/fem_element.cpp
    core/fem_solve.cpp
    io/read_file.cpp
    api/solver_api.cpp
)

set_target_properties(solver_core PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

target_include_directories(solver_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen
)

pybind11_add_module(_solver
    bindings.cpp
)

target_link_libraries(_solver PRIVATE solver_core)

install(
  TARGETS _solver
  LIBRARY DESTINATION solver
)

install(
  FILES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
  DESTINATION solver
)

if (CFS_BUILD_TESTS)
    find_package(Catch2 3 REQUIRED)

    add_executable(solver_tests
        tests/test_fem_basic.cpp
    )

    target_link_libraries(solver_tests
        PRIVATE
            solver_core
            Catch2::Catch2WithMain
    )

    add_test(NAME solver_tests COMMAND solver_tests)
endif()


if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(solver_core PRIVATE
        -fsanitize=address,undefined
        -fno-omit-frame-pointer
        -g
    )

    target_link_options(solver_core PRIVATE
        -fsanitize=address,undefined
    )

    target_link_options(solver_tests PRIVATE
        -fsanitize=address,undefined
    )

    target_compile_options(_solver PRIVATE
    -fsanitize=address,undefined
    -fno-omit-frame-pointer
    )

    target_link_options(_solver PRIVATE
        -fsanitize=address,undefined
    )
endif()