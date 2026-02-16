#define CATCH_CONFIG_MAIN
#include "solver_api.h"
#include <catch2/catch_test_macros.hpp>
#include <iostream>

TEST_CASE("Single element tension") {
  SolverInput in;

  in.nodes = {{0, 0}, {0, 5}, {0, 10}, {10, 0}, {10, 5}, {20, 0}};

  in.elements = {{0, 5, 2, 3, 4, 1}};

  in.fixed_dofs = {0, 1, 2, 4};
  in.fixed_values = {0, 0, 0, 0};

  in.forces = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10000, 0};

  SolverOutput out = solve_from_data(in);

  REQUIRE(out.disp[5][0] > 0.0);
}

TEST_CASE("Single element compression") {
  SolverInput in;

  in.nodes = {{0, 0}, {0, 5}, {0, 10}, {10, 0}, {10, 5}, {20, 0}};

  in.elements = {{0, 5, 2, 3, 4, 1}};

  in.fixed_dofs = {0, 1, 2, 4};
  in.fixed_values = {0, 0, 0, 0};

  in.forces = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -10000, 0};

  SolverOutput out = solve_from_data(in);

  REQUIRE(out.disp[5][0] < 0.0);
}