#include "solver_api.h"
#include "fem_solve.h"
#include "read_file.h"

SolverOutput solve_from_data(const SolverInput &input) { return solve(input); }

SolverOutput solve_from_file(const char *path) {
  SolverInput input = read_file(path);

  return solve(input);
}