#pragma once
#include <array>
#include <vector>

struct Node {
  double x;
  double y;
};

using Element = std::array<int, 6>;
using Nodes = std::vector<Node>;
using Elements = std::vector<Element>;

struct SolverInput {
  std::vector<Node> nodes;
  std::vector<Element> elements;
  std::vector<size_t> fixed_dofs;
  std::vector<double> fixed_values;
  std::vector<double> forces;
};

struct SolverOutput {
  std::vector<std::array<double, 3>> stress;
  std::vector<std::array<double, 2>> disp;
};
