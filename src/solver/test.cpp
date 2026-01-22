#include "fem.h"
#include <iostream>

int main() {

  std::vector<std::vector<double>> result = solve_from_file("fem.inp");

  for (std::vector<double> r : result) {
    std::cout << r[0] << " " << r[1] << " " << r[2] << std::endl;
  }

  return 0;
}