#include "Eigen/Core"
#include "read_file.h"
#include <cstddef>
#include <iostream>
#include <unordered_set>
using Eigen::MatrixXd;
using Eigen::VectorXd;

double area(Node i, Node j, Node k) {
  return 0.5 * (i.x * (j.y - k.y) + j.x * (k.y - i.y) + k.x * (i.y - j.y));
}

Nodes filter_by_index(const Nodes &data, const Elem &indices) {
  Nodes result;

  for (Index i : indices) {
    result.push_back(data[i]);
  }

  return result;
}

MatrixXd stiffness_matrix(double E, double v, double t, const Nodes &all_nodes,
                          const Elem &elem) {

  auto nodes = filter_by_index(all_nodes, elem);
  Node i = nodes[0];
  Node j = nodes[1];
  Node k = nodes[2];
  double A = area(i, j, k);

  MatrixXd B(3, 6);
  MatrixXd D(3, 3);
  MatrixXd K(6, 6);

  B << j.y - k.y, 0, k.y - i.y, 0, i.y - j.y, 0, 0, k.x - j.x, 0, i.x - k.x, 0,
      j.x - i.x, k.x - j.x, j.y - k.y, i.x - k.x, k.y - i.y, j.x - i.x,
      i.y - j.y;
  B *= 1 / (2 * A);

  D << 1, v, 0, v, 1, 0, 0, 0, (1 - v) / 2;
  D *= E / (1 - v * v);

  K = t * A * B.transpose() * D * B;
  return K;
}

void update_global_stiffness(MatrixXd &K, const MatrixXd &k,
                             const Elem &order) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      K(order[i] * 2, order[j] * 2) += k(i * 2, j * 2);
      K(order[i] * 2, order[j] * 2 + 1) += k(i * 2, j * 2 + 1);
      K(order[i] * 2 + 1, order[j] * 2) += k(i * 2 + 1, j * 2);
      K(order[i] * 2 + 1, order[j] * 2 + 1) += k(i * 2 + 1, j * 2 + 1);
    }
  }
}

VectorXd element_stress(double E, double v, double t, const Nodes &all_nodes,
                        const Elem &elem, const VectorXd &u_all) {
  Eigen::VectorXi ind_order(6);
  for (int i = 0; i < 3; ++i) {
    ind_order(2 * i) = static_cast<int>(elem[i] * 2);
    ind_order(2 * i + 1) = static_cast<int>(elem[i] * 2 + 1);
  }
  VectorXd u = u_all(ind_order);

  Node i = all_nodes[elem[0]];
  Node j = all_nodes[elem[1]];
  Node k = all_nodes[elem[2]];
  double A = area(i, j, k);

  MatrixXd B(3, 6);
  MatrixXd D(3, 3);
  MatrixXd K(6, 6);

  B << j.y - k.y, 0, k.y - i.y, 0, i.y - j.y, 0, 0, k.x - j.x, 0, i.x - k.x, 0,
      j.x - i.x, k.x - j.x, j.y - k.y, i.x - k.x, k.y - i.y, j.x - i.x,
      i.y - j.y;
  B *= 1 / (2 * A);

  D << 1, v, 0, v, 1, 0, 0, 0, (1 - v) / 2;
  D *= E / (1 - v * v);

  VectorXd s = D * B * u;
  return s;
}

int add(int a, int b) { return a + b; }

SolverInput make_input_from_file(const std::string &path) {
  return read_file(path);
};

SolverInput
make_input_from_parameters(const std::vector<std::vector<double>> &node_coords,
                           const std::vector<std::vector<int>> &elements,
                           const std::vector<int> &u_indices,
                           const std::vector<double> &u,
                           const std::vector<double> &f) {
  Nodes nodes;
  Elems elems;
  for (std::vector<double> n : node_coords) {
    nodes.push_back(Node{n[0], n[1]});
  }
  for (std::vector<int> e : elements) {
    elems.push_back(Elem{e[0], e[1], e[2]});
  }
  std::vector<size_t> u_indices_t(u_indices.begin(), u_indices.end());
  SolverInput input{nodes, elems, u_indices_t, u, f};
  return input;
};

double solve(const SolverInput &input) {
  double E = 210E6;
  double v = 0.3;
  double t = 1;
  // Node node1 = {0, 0};
  // Node node2 = {0, 10};
  // Node node3 = {20, 10};
  // Node node4 = {20, 0};
  // Node node5 = {30, 0};
  // Node node6 = {30, 10};
  // Nodes nodes = {node1, node2, node3, node4, node5, node6};
  // int dofs = nodes.size() * 2;
  // Elem elem1 = {0, 2, 1};
  // Elem elem2 = {0, 3, 2};
  // Elem elem3 = {2, 3, 5};
  // Elem elem4 = {3, 4, 5};
  // Elems elems = {elem1, elem2, elem3, elem4};
  // // std::vector<size_t> ind_zero = {0, 1, 2, 3, 8, 10};
  // // VectorXd F{{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  // // VectorXd u_initial{{0, 0, 0, 0, 0.014, 0.014}};
  // std::vector<size_t> ind_zero = {0, 1, 2, 3};
  // VectorXd F{{0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 100, 0}};
  // VectorXd u_initial{{0, 0, 0, 0}};
  Nodes nodes = input.nodes;
  Elems elems = input.elems;
  std::vector<size_t> u_indices = input.u_indices;
  std::vector<double> u_initial = input.u;
  VectorXd F = Eigen::Map<const VectorXd>(input.F.data(), input.F.size());
  // for (Node n : nodes) {
  //   std::cout << n.x << n.y << std::endl;
  // }
  // for (Elem e : elems) {
  //   std::cout << e[0] << e[1] << e[2] << std::endl;
  // }
  // for (size_t t : u_indices) {
  //   std::cout << t << std::endl;
  // }
  // for (double t : u_initial) {
  //   std::cout << t << std::endl;
  // }
  // std::cout << F << std::endl;
  int dofs = nodes.size() * 2;

  MatrixXd K = MatrixXd::Zero(dofs, dofs);
  for (Elem e : elems) {
    MatrixXd K1 = stiffness_matrix(E, v, t, nodes, e);
    update_global_stiffness(K, K1, e);
  }

  for (size_t k = 0; k < u_indices.size(); ++k) {
    int j = u_indices[k];
    F -= K.col(j) * u_initial[k];
  }

  std::unordered_set<size_t> remove(u_indices.begin(), u_indices.end());
  std::vector<size_t> ind_nonzero;
  for (int i = 0; i < dofs; ++i) {
    if (!remove.count(i)) {
      ind_nonzero.push_back(i);
    }
  }

  Eigen::VectorXi ind_nonzero_eig(ind_nonzero.size());
  for (int i = 0; i < ind_nonzero.size(); ++i) {
    ind_nonzero_eig(i) = static_cast<int>(ind_nonzero[i]);
  }

  MatrixXd K_nonzero = K(ind_nonzero_eig, ind_nonzero_eig);

  VectorXd F_nonzero = F(ind_nonzero_eig);

  VectorXd u_nonzero = K_nonzero.colPivHouseholderQr().solve(F_nonzero);

  VectorXd u = VectorXd::Zero(dofs);
  for (int i = 0; i < u_nonzero.size(); ++i)
    u(ind_nonzero_eig[i]) = u_nonzero(i);

  for (size_t k = 0; k < u_indices.size(); ++k)
    u(u_indices[k]) = u_initial[k];

  std::cout << "u\n";
  std::cout << u << "\n";
  std::cout << "\n";

  for (int i = 0; i < elems.size(); i++) {
    VectorXd S = element_stress(E, v, t, nodes, elems[i], u);
    std::cout << "stress element " << i << "\n";
    std::cout << S << "\n";
    std::cout << "\n";
  }

  return u_nonzero(0);
}

double solve_from_data(const std::vector<std::vector<double>> &node_coords,
                       const std::vector<std::vector<int>> &elements,
                       const std::vector<int> &u_indices,
                       const std::vector<double> &u,
                       const std::vector<double> &f) {
  SolverInput input =
      make_input_from_parameters(node_coords, elements, u_indices, u, f);
  double res = solve(input);
  return res;
};

double solve_from_file(const std::string &path) {
  SolverInput input = make_input_from_file(path);
  double res = solve(input);
  return res;
};
