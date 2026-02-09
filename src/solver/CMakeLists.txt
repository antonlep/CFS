cmake_minimum_required(VERSION 3.19)
project(solver LANGUAGES CXX)

# -------------------------
# pybind11 (always fetched, never manual)
# -------------------------
find_package(pybind11 CONFIG REQUIRED)

include(FetchContent)

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v3.0.1
)

option(CFS_FETCH_PYBIND11 "Fetch pybind11 if not found" ON)

if (CFS_FETCH_PYBIND11)
  FetchContent_MakeAvailable(pybind11)
endif()

# -------------------------
# Core C++ library
# -------------------------
add_library(solver_core
    core/fem_element.cpp
    core/fem_solve.cpp
    core/householder.cpp
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

target_compile_features(solver_core PUBLIC cxx_std_20)

# -------------------------
# Python extension module
# -------------------------
pybind11_add_module(_solver
    bindings.cpp
)

target_link_libraries(_solver PRIVATE solver_core)

# -------------------------
# Install for Python
# -------------------------
install(
  TARGETS _solver
  LIBRARY DESTINATION solver
)

install(
  FILES __init__.py
  DESTINATION solver
)

# -------------------------
# Tests (optional)
# -------------------------
if (CFS_BUILD_TESTS)
    include(FetchContent)

    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.5.2
    )
    FetchContent_MakeAvailable(Catch2)

    add_executable(solver_tests tests/test_fem_basic.cpp)

    target_link_libraries(solver_tests
        PRIVATE
            solver_core
            Catch2::Catch2WithMain
    )

    include(Catch)
    catch_discover_tests(solver_tests)
endif()

add_executable(solver_cli
    solver_cli.cpp
)

target_link_libraries(solver_cli
    PRIVATE solver_core
)

target_include_directories(solver_cli PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# -------------------------
# Debug sanitizers
# -------------------------
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (UNIX AND NOT APPLE)
        # Linux / WSL: sanitizers OK
        foreach(tgt solver_core _solver)
            target_compile_options(${tgt} PRIVATE
                -g
                -fno-omit-frame-pointer
                -fsanitize=address,undefined
            )
            target_link_options(${tgt} PRIVATE
                -fsanitize=address,undefined
            )
        endforeach()
    else()
        # Windows: NO sanitizers
        message(STATUS "Debug build on Windows: sanitizers disabled")
    endif()
endif()
