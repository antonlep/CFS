#include "solver_api.h"
#include <iostream>

int main(int argc, char **argv) {

  if (argc < 2) {
    std::cerr << "usage: solver_cli input.inp\n";
    return 1;
  }

  SolverOutput result = solve_from_file(argv[1]);

  std::cout << "Stress" << "\n";
  std::cout << "\n";
  for (std::array<double, 3> r : result.stress) {
    std::cout << r[0] << " " << r[1] << " " << r[2] << std::endl;
  }
  return 0;
}