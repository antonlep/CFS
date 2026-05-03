#pragma once
#include <Eigen/Core>
#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

// Index cast helpers
inline Eigen::Index eidx(std::size_t i) { return static_cast<Eigen::Index>(i); }
inline std::size_t uidx(int i) {
  assert(i >= 0);
  return static_cast<std::size_t>(i);
}
inline int sidx(std::size_t i) { return static_cast<int>(i); }

struct Node {
  double x;
  double y;
};

using Element = std::array<std::size_t, 6>;
using Nodes = std::vector<Node>;
using Elements = std::vector<Element>;

struct BoundaryEdge {
  size_t n1, n2, n3; // quadratic edge nodes
  double h;
  double Tinf;
};

struct TractionEdge {
  size_t n1, n2, n3; // quadratic edge nodes
  double tx, ty;     // traction components [force/length]
};

struct MaterialProperties {
  double E = 210e6;
  double nu = 0.3;
  double t = 1.0;
  double k = 45.0;
  double alpha = 12e-6;
  double T0 = 273.0;
};

struct SolverInput {
  std::vector<Node> nodes;
  std::vector<Element> elements;
  std::vector<size_t> fixed_dofs;
  std::vector<double> fixed_values;
  std::vector<double> forces;
  std::vector<BoundaryEdge> boundary_edges;
  std::vector<TractionEdge> traction_edges;
  MaterialProperties material;
};

struct SolverOutput {
  std::vector<std::array<double, 3>> stress;
  std::vector<std::array<double, 2>> disp;
  std::vector<double> temperature;
};
