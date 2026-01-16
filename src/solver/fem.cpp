#include "Eigen/Core"
#include <Eigen/Dense>
#include <iostream>
#include <unordered_set>
#include <vector>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Index = Eigen::Index;
using Dir = std::size_t;
using NodeId = Index;
struct Node {
  double x;
  double y;
};
struct Disp {
  Index dir;
  double mag;
};
using Nodes = std::vector<Node>;
using Elem = std::vector<NodeId>;
using Elems = std::vector<Elem>;
using Dirs = std::vector<Dir>;

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

double plate(double h, double w, double t, double f) {
  double E = 210E6;
  double v = 0.3;
  f = 0.5 * f;
  Nodes nodes = {Node{0, 0},  Node{0, 10}, Node{20, 10},
                 Node{20, 0}, Node{30, 0}, Node{30, 10}};
  int dofs = nodes.size() * 2;
  Elem elem1 = {0, 2, 1};
  Elem elem2 = {0, 3, 2};
  Elem elem3 = {2, 3, 5};
  Elem elem4 = {3, 4, 5};
  Elems elems = {elem1, elem2, elem3, elem4};
  std::vector<size_t> ind_zero = {0, 1, 2, 3};
  VectorXd F{{0, 0, 0, 0, 0, 0, 0, 0, f, 0, f, 0}};

  MatrixXd K = MatrixXd::Zero(dofs, dofs);
  for (Elem e : elems) {
    MatrixXd K1 = stiffness_matrix(E, v, t, nodes, e);
    update_global_stiffness(K, K1, e);
  }

  std::unordered_set<size_t> remove(ind_zero.begin(), ind_zero.end());
  std::vector<size_t> ind_nonzero;
  for (int i = 0; i < dofs; ++i) {
    if (!remove.count(i)) {
      ind_nonzero.push_back(i);
    };
  };

  Eigen::VectorXi ind_nonzero_eig(ind_nonzero.size());
  for (int i = 0; i < ind_nonzero.size(); ++i) {
    ind_nonzero_eig(i) = static_cast<int>(ind_nonzero[i]);
  }

  MatrixXd K_nonzero = K(ind_nonzero_eig, ind_nonzero_eig);

  VectorXd F_nonzero = F(ind_nonzero_eig);

  VectorXd u_nonzero = K_nonzero.colPivHouseholderQr().solve(F_nonzero);

  VectorXd u = VectorXd::Zero(dofs);
  for (int i = 0; i < u_nonzero.size(); i++) {
    u(ind_nonzero_eig[i]) = u_nonzero(i);
  };
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
