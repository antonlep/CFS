#ifndef FEM_H
#define FEM_H
#include "read_file.h"

int add(int a, int b);

std::vector<std::vector<double>>
solve_from_data(const std::vector<std::vector<double>> &node_coords,
                const std::vector<std::vector<int>> &elements,
                const std::vector<int> &u_indices, const std::vector<double> &u,
                const std::vector<double> &f);

std::vector<std::vector<double>> solve_from_file(const std::string &path);

#endif